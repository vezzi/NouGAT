from __future__ import absolute_import
import click
import argparse
import os
from collections import OrderedDict
import yaml


class State(object):
    def __init__(self):        
        self.params = OrderedDict() 
        self.params = {}


pass_state = click.make_pass_decorator(State, ensure=True)
CONTEXT_SETTINGS = dict(default_map=None)


def callback(ctx, param, value):
    state = ctx.ensure_object(State)
    if value is not None:
        state.params[param.name] = value
    return value


def _to_namespace(args_dict, **kwargs):
    args = argparse.Namespace()
    d = vars(args)
    for key, val in args_dict.items():
        d[key] = val
    for key, val in kwargs.items():
        d[key] = val
    return args


def common_options(f):
    f = click.option('--dry-run', is_flag=True, default=False, expose_value=False, callback=callback,
            help='A dry-run (of sorts). Will create .yaml and .slurm but not run the analysis')(f)
    f = click.option('--global-config', expose_value=False, type=click.Path(exists=True), required=True, callback=callback,
            help='Path to NouGAT config file')(f)
    return f


def slurm_options(f):
    f = click.option('--time', expose_value=False, callback=callback, default="1-00:00:00",
            help="SLURM: Wall time for the slurm job (eg. 1-00:00:00 for 1 day.)", 
            show_default=True)(f)
    f = click.option('--project', expose_value=False, callback=callback, default="a2012043",
            help="SLURM: Uppnex project", show_default=True)(f)
    f = click.option('--qos', expose_value=False, callback=callback, default="seqver",
            help="SLURM: Quality of service", show_default=True)(f)
    f = click.option('--env', expose_value=False, callback=callback, default="DeNovoPipeline",
            help="SLURM: Anaconda environment to be used", show_default=True)(f)
    f = click.option('--threads', expose_value=False, callback=callback, default=16,
            help="SLURM: Number of cores to request", type=int, show_default=True)(f)
    f = click.option('--email', expose_value=False, callback=callback, default="",
            help="SLURM: Email address for notifications", show_default=True)(f)
    return f


def lib_options(f):
    f = click.option('--sample-data-dir', expose_value=False, required=True, type=click.Path(exists=True), callback=callback,
            help=("LIBRARY: Path to directory (usually INBOX) containing the project "
            "(one dir per sample, scilife structure project/sample/flowcell/)"))(f)
    f = click.option('--insert', expose_value=False, required=True, callback=callback,
            help="LIBRARY: insert size", type=int)(f) 
    f = click.option('--std', expose_value=False, required=True, callback=callback,
            help="LIBRARY: Insert size standard deviation", type=int)(f)
    f = click.option('--orientation', type=click.Choice(['innie', 'outtie']), expose_value=False, required=True, callback=callback,
            help="LIBRARY: orientation, ie. 'innie' for paired-end and 'outtie' for mate-pairs")(f)
    return f


@click.group(context_settings=CONTEXT_SETTINGS)
@click.pass_context
def cli(ctx):
    """A set of convenience functions for the bioinformaticians at NGI Stockholm for using
    NouGAT. This assumes that data is in a specific format and that is is run on a SLURM cluster"""   
    pass


@cli.command()
@common_options
@slurm_options
@lib_options
@click.option('--adapter', required=True, type=click.Path(exists=True), help="Adapter sequences to remove")
@click.option('--reference', type=click.Path(exists=True), help="Fasta reference file. Will run BWA if specified")
@pass_state
def qc_analysis(state, **kwargs):
    """Perform QC for de novo / mate-pair libraries"""
    from . import run_QC_analysis
    run_QC_analysis.main(_to_namespace(state.params, **kwargs))


@cli.command()
@common_options
@click.option('--QC-folder', type=click.Path(exists=True), required=True, help="Path to where denovo QC pipeline was run. Where you find, "
    "the sample folders P????_??? eg. in, /foo/B.Baz_15_01/01-QC/")
@pass_state
def qc_report(state, **kwargs):
    """Generate a PDF report detailing the QC analysis results. 
    This will also create a results folder for each sample, which 
    contains the QC data to be delivered to the platform user."""
    from . import run_QC_report
    run_QC_report.main(_to_namespace(state.params, **kwargs))


@cli.command()
@common_options
@slurm_options
@lib_options
@click.argument('assemblers', required=True, nargs=-1, type=str)
@click.option('--kmer', required=True, type=int,
        help="K-mer size. Set to a reasonabe value or estimate on using Kmergenie")
@click.option('--genomesize', required=True, type=int ,help="The size of the genome (in bp)")
@click.option('--afterqc', is_flag=True, default=False, 
        help="The data sample data is the output of the QC pipeline. --sample-data-dir should be set accordingly")
@click.option('--keep-tmp-files', is_flag=True, default=False,
        help="This will override the default behaviour of deleting temporary files, eg. assembly graphs and error corrections.")
@pass_state
def assembly(state, **kwargs):
    """Start assembly using with sample sequenced with a single
    library. For samples sequenced with multiple libraries the user
    needs to prepare the sample_config file and run the
    deNovoPipeline manually. If multiple samples are present in the
    specified folder then one assembly for each sample will be
    performed. If a sample is splitted across multiple runs all the
    data willl be used"""
    from . import run_assemblies
    run_assemblies.main(_to_namespace(state.params, **kwargs))


@cli.command()
@common_options
@slurm_options
@click.option('--assembly-dir', type=click.Path(exists=True), required=True, help="Path to directory containg assemblies")
@click.option('--multiple-lib-project', is_flag=True, default=False, help="To be specified if we are running a "
            "mulitple library assembly")
@pass_state
def validation(state, **kwargs):
    """Runs assembly validation with sample sequenced with a single
    library. For samples sequenced with multiple libraries the user
    needs to prepare the sample_config file and run the
    deNovoPipeline manually. If multiple samples are present in the
    specified folder then one assembly for each sample will be
    performed. If a sample is split across multiple runs all the
    data will be used"""
    from . import run_validation
    run_validation.main(_to_namespace(state.params, **kwargs))


@cli.command()
@common_options
@click.option('--validation-dirs', type=click.Path(exists=True), required=True, help="Directory where validation are stored for each sample")
@click.option('--assemblies-dirs', type=click.Path(exists=True), required=True, help="Directory where assemblies are stored for each sample")
@click.option('--min-contig-length', default=1000, type=int, help="Minimum contig length to be contig length for contiguity stats", show_default=True)
@click.option('--sample-name', help="In case of multi-project sample (ie. PE+MP) specify a sample name. Must be specified when --no-uppmax option is used.")
@click.option('--no-uppmax', is_flag=True, default=False, help="if specified the validation-dir and the assemblies-dir "
                        "is assumed to contain the assemblies (and not samples)")
@pass_state
def report(state, **kwargs):
    """Write a pdf report and generate folder structure for delivery
    It is assumed that assemble part and evalaution part have been run with this
    pipeline, otherwise the assumptions done on the file names and on the results
    """
    from . import run_assembly_report
    run_assembly_report.main(_to_namespace(state.params, **kwargs))


def main():
    try:
        conf_file = os.path.join(os.environ.get('HOME'), '.nougat','scilifelab.conf')
        with open(conf_file) as f:
            config = yaml.load(f)
    except IOError:
        click.secho("Could not open the config file {}, will use hardcoded defaults".format(conf_file), fg="red")
    else:   
        CONTEXT_SETTINGS["default_map"] = config
    cli()

