from __future__ import absolute_import
import sys, os, yaml, glob
import subprocess
import pandas as pd
import argparse
import re
import string
import shutil
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from nougat import common
from nougat import pdf
from nougat.pdf.theme import colors, DefaultTheme


def main(args):
    workingDir = os.getcwd()
    validation_dirs = args.validation_dirs
    assemblies_dirs = args.assemblies_dirs
    min_contig_length = args.min_contig_length
    # store all samples folders
    assemblies_samples_dirs = [sample for sample in \
            os.listdir(assemblies_dirs) \
            if os.path.isdir(os.path.join(assemblies_dirs,sample))]
    # store all samples folders
    validation_samples_dirs = [sample for sample in \
            os.listdir(validation_dirs) \
            if os.path.isdir(os.path.join(validation_dirs,sample))]
    if args.no_uppmax:
        collect_results_and_report(validation_dirs, assemblies_dirs, "",
                args.sample_name, min_contig_length, args.no_uppmax)
    else:
        for sample in assemblies_samples_dirs:
            if sample not in validation_samples_dirs:
                system.exit("ATTENTION: sample {} is present in dir {} but \
                        absent in dir {}".format(sample, assemblies_dirs,
                        validation_dirs))
            else: #otherwise I have one or more sample folders
                 #It means that I have assembled with a set of assemblers and
                 #all of them have gone through validation... I can proceed
                os.chdir(workingDir)
                validation_sample_dir = os.path.join(validation_dirs, sample)
                assemblies_sample_dir = os.path.join(assemblies_dirs, sample)
                sample_folder = os.path.join(workingDir, sample)
                if not os.path.exists(sample_folder):
                    os.makedirs(sample_folder)
                os.chdir(sample_folder)
                collect_results_and_report(validation_sample_dir,
                        assemblies_sample_dir, sample_folder, sample,
                        min_contig_length, args.no_uppmax)


def collect_results_and_report(validation_sample_dir, assemblies_sample_dir,
        sample_folder, sample, min_contig_length, no_uppmax):
    assemblers_validation = [assembler for assembler in \
            os.listdir(validation_sample_dir) \
            if os.path.isdir(os.path.join(validation_sample_dir,assembler))]
    assemblers_assemblies = [assembler  for assembler in \
            os.listdir(assemblies_sample_dir) \
            if os.path.isdir(os.path.join(assemblies_sample_dir,assembler))]
    if set(assemblers_validation) != set(assemblers_assemblies):
        sys.exit("Error: different assemblies in assemblies and validation "
                "folder: {} {}".format(
                assemblies_sample_dir,validation_sample_dir))
    assemblers_assemblies = sorted(assemblers_assemblies)
    assemblers_validation = sorted(assemblers_validation)
    #now I am in the folder that will contain all the results
    if not os.path.exists("assemblies"):
        os.makedirs("assemblies")
        for assembler in assemblers_assemblies:
            shutil.copytree(os.path.join(assemblies_sample_dir, assembler),
                    os.path.join("assemblies", assembler))
    ##copied original assemblies
    if not os.path.exists("evaluation"):
        os.makedirs("evaluation")

    # let us start with QAcompute results
    QA_pictures = os.path.join(sample_folder, "evaluation", "QA_pictures")
    if not os.path.exists(QA_pictures):
        os.makedirs(QA_pictures)
    picturesQA = {}
    for assembler in assemblers_assemblies:
        original_QA_dir = os.path.join(validation_sample_dir,
                    assembler, "QAstats")
        if not os.path.exists(original_QA_dir):
            continue
        cur_ass_dir = os.path.join(QA_pictures, assembler)
        if not os.path.exists(cur_ass_dir):
            os.makedirs(cur_ass_dir)
        ##QA pictures
        Original_CoverageDistribution200 = os.path.join(validation_sample_dir,
                assembler, "QAstats", "Coverage_distribution_noOutliers.png")
        Original_GC_vs_Coverage = os.path.join(validation_sample_dir,
                assembler, "QAstats", "GC_vs_Coverage_noOutliers.png")
        Original_GC_vs_CtgLength = os.path.join(validation_sample_dir,
                assembler, "QAstats", "GC_vs_CtgLength.png")
        Original_MedianCov_vs_CtgLength = os.path.join(validation_sample_dir,
                assembler, "QAstats", "MedianCov_vs_CtgLength_noOutliers.png")
        Original_QAstats_gc_result = ""
        QAstats_gc_result_file_name = ""
        if no_uppmax:
            QAstats_gc_result_file_name = [name for name in \
                    os.listdir(os.path.join(validation_sample_dir,
                    assembler, "QAstats")) \
                    if name.endswith(".bam.cov.gc")][0]
            Original_QAstats_gc_result = os.path.join(validation_sample_dir,
                    assembler, "QAstats",
                    "{}".format(QAstats_gc_result_file_name))
        else:
            Original_QAstats_gc_result = os.path.join(validation_sample_dir,
                    assembler, "QAstats", "{}.bam.cov.gc".format(sample))

        Copied_CoverageDistribution200 = os.path.join(cur_ass_dir,
                "Coverage_distribution_noOutliers.png")
        Copied_GC_vs_Coverage = os.path.join(cur_ass_dir,
                "GC_vs_Coverage_noOutliers.png")
        Copied_GC_vs_CtgLength = os.path.join(cur_ass_dir,
                "GC_vs_CtgLength.png")
        Copied_MedianCov_vs_CtgLength = os.path.join(cur_ass_dir,
                "MedianCov_vs_CtgLength_noOutliers.png")
        if no_uppmax:
            Copied_QAstats_gc_result = os.path.join(cur_ass_dir,
                    QAstats_gc_result_file_name)
        else:
            Copied_QAstats_gc_result = os.path.join(cur_ass_dir,
                    "{}.bam.cov.gc".format(sample))

        shutil.copy(Original_CoverageDistribution200,
                Copied_CoverageDistribution200)
        shutil.copy(Original_GC_vs_Coverage,
                Copied_GC_vs_Coverage)
        shutil.copy(Original_GC_vs_CtgLength,
                Copied_GC_vs_CtgLength)
        shutil.copy(Original_MedianCov_vs_CtgLength,
                Copied_MedianCov_vs_CtgLength)
        shutil.copy(Original_QAstats_gc_result,
                Copied_QAstats_gc_result)

        picturesQA[assembler] = [[Copied_CoverageDistribution200,
            "Contig coverage distribtion" ],
            [Copied_GC_vs_Coverage, "GC-content versus contig-coverage"],
            [Copied_GC_vs_CtgLength, "GC-content versus contig-Length"],
            [Copied_MedianCov_vs_CtgLength,
            "Median-coverage vs Contig-Length"]]

    # now FRCurve results
    FRC_folder = os.path.join(sample_folder, "evaluation", "FRCurves")
    if not os.path.exists(FRC_folder):
        os.makedirs(FRC_folder)
    Features = ["_FRC", "COMPR_MP_FRC", "COMPR_PE_FRC", "HIGH_COV_PE_FRC",
            "HIGH_NORM_COV_PE_FRC", "HIGH_OUTIE_MP_FRC", "HIGH_OUTIE_PE_FRC",
            "HIGH_SINGLE_MP_FRC", "HIGH_SINGLE_PE_FRC", "HIGH_SPAN_MP_FRC",
            "HIGH_SPAN_PE_FRC", "LOW_COV_PE_FRC", "LOW_NORM_COV_PE_FRC",
            "STRECH_MP_FRC", "STRECH_PE_FRC"]

    FRC_to_print = ""
    for feature in Features:
        FRCurves = []
        for assembler in assemblers_assemblies:
            FRCurve_Orig_name = os.path.join(validation_sample_dir, assembler,
                    "FRCurve", "{}{}.txt".format(sample, feature))
            FRCurves.append([assembler, FRCurve_Orig_name])
        FRCname = _plotFRCurve(os.path.join(FRC_folder, "{}_{}".format(sample,
            feature)), FRCurves)
        if feature == "_FRC":
            FRC_to_print = FRCname
        FRCurves = []

    #Contiguity stats
    contig_stats = []
    source_stats = []
    for assembler in assemblers_assemblies:
        asm_stats = [assembler]
        stat_file = os.path.join(validation_sample_dir, assembler,
                "contig_stats", "contiguity.out")
        source_stats.append((stat_file, assembler))
        with open(stat_file, "r") as sf:
            for line in sf:
                if line.startswith("scaffolds"):
                    sl = line.strip().split("\t")
                    asm_stats.extend([sl[1], sl[4], sl[8], sl[9], sl[10], sl[2], sl[6]])
        contig_stats.append(asm_stats)
    contig_stats = sorted(contig_stats, key=lambda x: x[0])

    # Copy the contig stats
    stat_folder = os.path.join(sample_folder, "evaluation", "contig_stats")
    if not os.path.exists(stat_folder):
        os.makedirs(stat_folder)
    for src_stat in source_stats:
        shutil.copy(src_stat[0], os.path.join(stat_folder, "{}.contiguity.out".format(src_stat[1])))


    #### now I can produce the report
    write_report(sample_folder, sample, assemblies_sample_dir,
            assemblers_assemblies,  picturesQA, FRC_to_print,
            min_contig_length, contig_stats)
    return


def write_report(sample_folder, sample, assemblies_sample_dir, assemblers,
        picturesQA, FRCname ,min_contig_length, contig_stats):
    """This function produces a pdf report """

    with open(os.path.join(assemblies_sample_dir,
        "{}_assemble.yaml").format(sample)) as sample_config_handle:
        sample_config = yaml.load(sample_config_handle)

    reportDir   = os.path.join(sample_folder, "report")
    if not os.path.exists(reportDir):
        os.makedirs(reportDir)
    PDFtitle = os.path.join(sample_folder, "report",
            "{}_assembly_report.pdf".format(sample))

    # this you cannot do in rLab which is why I wrote the helper initially
    TABLE_WIDTH = 540
    class MyTheme(DefaultTheme):
        doc = {
            'leftMargin': 25,
            'rightMargin': 25,
            'topMargin': 20,
            'bottomMargin': 25,
            'allowSplitting': False
            }
    # let's create the doc and specify title and author
    doc = pdf.Pdf('{}'.format(sample),
            'NGI-Stockholm, Science for Life Laboratory')
    # now we apply our theme
    doc.set_theme(MyTheme)
    # give me some space
    doc.add_spacer()
    # this header defaults to H1
    scriptDirectory = os.path.split(os.path.abspath(__file__))[0]
    logo_path = os.path.join(scriptDirectory,
            '../pictures/ngi_scilife.png')
    doc.add_image(logo_path, 540, 50, pdf.CENTER)
    # give me some space
    doc.add_spacer()

    doc.add_header('NGI-Stockholm -- Science For Life Laboratory')
    doc.add_header('De Novo Assembly Best-Practice Analysis Report')
    doc.add_header('{}'.format(sample))

    doc.add_spacer()
    doc.add_paragraph("For sample {}  NGI-Stockholm best-practice analysis "
            "for de novo assembly and assembly evaluation have been "
            "performed.".format(sample))
    doc.add_paragraph("Sample has been assembled using the following "
            "assemblers:")
    bollet_list = []
    for assembler in assemblers:
        bollet_list.append("{}".format(assembler)) # TODO: add version
    doc.add_list(bollet_list)
    doc.add_spacer()
    doc.add_paragraph("Each assembler has been evaluated by aligning a subset "
            "of Illumina reads back to the assembled sequence. Consequently "
            "statistics on coverage, GC-content, contig length distribution "
            "have been computed. A global ranking with FRCurve is also "
            "performed.")

    doc.add_spacer()
    doc.add_paragraph("De novo assembly and de novo assembly evaluation are "
            "two difficult computational exercises. Currently there is no tool "
            "(i.e., de novo assembler) that is guaranteed to always outperform "
            "the others. Many recent publications (e.g., GAGE, GAGE-B, "
            "Assemblathon 1 and 2) showed how the same assembler can have "
            "totally different performances on slightly different datasets. "
            "For these reasons, at NGI-Stockholm we do not limit our de novo "
            "analysis to a single tool, instead we employ several assemblers "
            "and we provide our users with a semi-automated evaluation in "
            "order to allow them to choose the best assembler based on their "
            "specific needs. The assembly or assemblies judged to be the best "
            "can be directly employed to answer important biological "
            "questions, or they can be used as a backbone for a specific user "
            "defined assembly pipeline (i.e., use of extra data, use of non "
            "supported tools, variation of parameters).")

    doc.add_spacer()
    doc.add_paragraph("For each assembly the following information is "
            "provided:")
    doc.add_list(["Table with Standard Assembly Statistics: number of "
            "scaffolds, number of scaffold longer than {}bp, N50, N80, length "
            "of the longest scaffold, total assembly length, and sum of "
            "scaffolds scaffolds longer than "
            "{}bp.".format(min_contig_length,min_contig_length),
            "For each individual assembler four plots are automatically "
            "generated: Contig-coverage distribution, GC-content versus "
            "Contig-Coverage, GC-content vs Contig-Length, and Contig "
            "Coverage vs Contig Length.",
            "FRCurve plot: ROC-curve inspired method for assembly "
            "validation."])

    doc.add_spacer()
    doc.add_paragraph("Only contigs longer than {}bp are used during "
            "validation. This is done in order to avoid problems with outlier "
            "points and to partially circumvent the fact that some assemblers "
            "output small contigs while others perform automatic trimming. "
            "Statistics like N50, N80, etc. are computed on the expected "
            "genome length in order to normalise the numbers and allow a "
            "fair comparison among various tools.".format(min_contig_length))

    doc.add_paragraph("Coverage information and FRCurve features are obtained "
            "by aligning the same reads used in the assembling phase against "
            "the assembled sequences using BWA-MEM algorithm.")

    doc.add_paragraph("This report is delivered both via e-mail and via "
            "Uppmax. In particular on Uppmax the following files are available "
            "for further result inspection:")
    doc.add_list([
            "The report saved in the folder report.",
            "All the assemblies (ie., contigs and scaffolds).",
            "All the evaluations. For QA pictures, the same pictures included "
            "in this report plus the original cov.gc table are present. For "
            "FRC, the FRCurve for each invidual feature is plotted"])

    doc.add_paragraph("Please note that the pipeline generates all the plots "
            "automatically. The pipeline tries to eliminate outliers in order "
            "to visualize data in a meaningful and useful way (e.g., a single "
            "contig with extremely high coverage can jeopardize the "
            "visualization of all the other contigs). However, there might be "
            "situations where interesting points are discarded. We recommend "
            "to always inspect the original tables that are delivered on "
            "Uppmax altogether with this report.")
    doc.add_pagebreak()
    doc.add_header("Standard Contiguity Metrics", pdf.H2)
    doc.add_paragraph("Contiguity measures give an idea of the connectivity "
            "of the assembly. For all the assemblies that have been generated "
            "we report:")
    doc.add_list([
            "n. scaff: number of scaffolds produced by the assembler.",
            "n. scaff>{}: number of scaffolds produced by the assembler "
            "longer or equal to {}bp.".format(
            min_contig_length, min_contig_length),
            "NG50: length of the longest scaffold such that the sum of all the "
            "scaffolds longer than it is at least 50% of the estimated genome "
            "length.",
            "NG80: length of the longest scaffold such that the sum of all the "
            "scaffolds longer than it is at least 80% of the estimated genome "
            "length.",
            "max_scf_lgth: maximum scaffold length.",
            "Ass_length: total assembly length.",
            "Ass_lgth_scfs>{}: total assembly length considering only "
            "scaffolds longer or equal to {}bp.".format(
            min_contig_length,min_contig_length)])
    doc.add_spacer()

    assemblyStats = [['assembler', 'n. scaff',
        'n. scaff>{}'.format(min_contig_length), 'NG50', 'NG80', 'max_scf_lgth',
        'Ass_lgth', 'Ass_lgth_scfs>{}'.format(min_contig_length)]]
    assemblyStats.extend(contig_stats)

    doc.add_spacer()
    doc.add_table(assemblyStats, TABLE_WIDTH)
    doc.add_spacer()
    doc.add_paragraph("The table can be used to give an idea of the "
            "conectivity of the assemblies. N.B. a high connectivity "
            "(i.e., a long N50) does not necessarily implies a high quality "
            "assembly as this assembly might be the result of a too "
            "aggressive strategy.")

    #Now QC pictures
    doc.add_header("QC plots", pdf.H2)
    doc.add_paragraph("QC plots try to represent the relations between "
            "coverage, GC content, and scaffolds length in a visual way. "
            "The idea is to use the four plots to see if the assembly "
            "coincides with the expected results and to check the presence of "
            "biases. For each assembler we aligned the same reads used for de "
            "novo assembly back to the assembly itself and we compute for each "
            "scaffold its coverage, its GC content, and its length. In this "
            "way the following plots can be generated:")
    doc.add_list(["Contig Coverage Distribution: this plot shows the scaffold "
            "coverage distribution. Ideally the picture should look like a "
            "gaussian distribution with the maximum around the expected "
            "coverage. If the assembly is highly connected (i.e., formed by "
            "only tens of scaffolds) this shape might be not visible. ",
            "GC-Content versus Contig-Coverage: this plot shows for each "
            "scaffold its GCs content on the x-axis and its coverage on the "
            "y-axis. Typically scaffolds should cluster forming a cloud. "
            "The presence of two distict clouds might suggest the presence of "
            "contamination.",
            "GC-Content versus Contig-Length: this plot shows for each "
            "scaffold its GCs content on the x-axis and its length on the "
            "y-axis. The main purpose is to identify biases towards long or "
            "short scaffolds.",
            "Median-Coverage vs Contig-Length: this plot shows for each "
            "scaffold its coverage on the x-axis and its length on the y-axis. "
            "The main purpose is to identify biases towards long or short "
            "scaffolds."])

    doc.add_paragraph("N.B.: pictures are produced automatically thus it "
            "might be possible that some outliers are not printed. The "
            "original cov.gc table is present under the "
            "evaluation/QA_pictures folder")

    for assembler, assembler_QC_pictures in picturesQA.iteritems():
        doc.add_pagebreak() 
        doc.add_header("QC plots for {}".format(assembler) , pdf.H3)
        doc.add_image(assembler_QC_pictures[0][0], 280, 200, pdf.CENTER,
                "Contig Coverage Distribution")
        doc.add_image(assembler_QC_pictures[1][0], 280, 200, pdf.CENTER,
                "GC-Content versus Contig-Coverage")
        doc.add_image(assembler_QC_pictures[2][0], 280, 200, pdf.CENTER,
                "GC-Content versus Contig-Length")
        doc.add_image(assembler_QC_pictures[3][0], 280, 200, pdf.CENTER,
                "Median-Coverage vs Contig-Length")

    #FRCurves
    doc.add_pagebreak()
    doc.add_header("FRCurves", pdf.H2)
    doc.add_paragraph("Inspired by the standard receiver operating "
            "characteristic (ROC) curve, the Feature-Response curve (FRCurve) "
            "characterizes the sensitivity (coverage) of the sequence "
            "assembler output (contigs) as a function of its discrimination "
            "threshold (number of features/errors). Generally speaking, "
            "FRCurve can be used to rank different assemblies: the sharper "
            "the curve is the better the assembly is (i.e., given a certain "
            "feature threshold, we prefer the assembler that reconstructs a "
            "higher portion of the genome with fewer features). FRCurve is one "
            "of the few tools able to evaluate de novo assemblies in absence "
            "of a reference sequence. Results are not always straightforward "
            "to interpret and must be always used in conjunction with other "
            "sources (e.g., quality plots and standard assembly statistics).")
    doc.add_image(FRCname, 336, 220, pdf.CENTER, "Feature-Response curve. "
            "On x-axis is the number of features in total, on y-axis is "
            "coverage (based on estimated genome size).")
    doc.render(PDFtitle)
    return 0


def _plotFRCurve(outputName, FRCurves):
    marks = [None,'o', 's', 'v', '>', 'd', 'x', '<', '^', '+']

    # These are the "Tableau 20" colors as RGB.
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)
    plt.rcParams['lines.linewidth'] = 2.0

    FRCurveName = "{}_FRCurve.png".format(outputName)
    maxXvalues   = []
    for i, FRCurveData in enumerate(FRCurves):
        assembler = FRCurveData[0]
        FRCurve = FRCurveData[1]
        FRC_data = pd.io.parsers.read_csv(FRCurve, sep=' ', header=None)
        FRC_features = FRC_data[FRC_data.columns[0]].tolist()
        FRC_coverage = FRC_data[FRC_data.columns[1]].tolist()
        plt.plot(FRC_features, FRC_coverage, label="{}".format(assembler),
                 color=tableau20[i], marker=marks[i], markersize=5, markeredgecolor='k', markevery=2)
        maxXvalues.append(max(FRC_features))
    maxXvalues.sort()
    maxXvalues.reverse()
    maxXvalue = maxXvalues[0]
    for i in range(1, len(maxXvalues)-1):
        if maxXvalue > maxXvalues[i]*2 \
                and (maxXvalues[i-1] - maxXvalues[i] > 1000):
            maxXvalue = maxXvalues[i] + int(maxXvalues[i]*0.10)

    plt.ylim((-5,140))
    plt.xlim((-1,maxXvalue))
    plt.legend(loc=4, ncol=1, borderaxespad=0.)
    plt.savefig(FRCurveName)
    plt.clf()
    return FRCurveName


if __name__ == '__main__':
    parser = argparse.ArgumentParser("This utility scripts will generate the "
            "report files for each assembled sample. It is assumed that "
            "assemble part and evalaution part have been run with this "
            "pipeline, otherwise the assumptions done on the file names "
            "and on the results")
    parser.add_argument('--validation-dirs', type=str, required=True,
            help=("Directory where validation are stored  for each sample "
            "(one assembler per folder)"))
    parser.add_argument('--assemblies-dirs', type=str, required=True,
            help=("Directory where assemblies are stored for each sample "
            "(one assembler per folder)"))
    parser.add_argument('--min-contig-length', type=int, default=1000,
            help=("minimum length that a contig must have to be considered "
            "long"))
    parser.add_argument('--global-config', type=str,
            help="global configuration file")
    parser.add_argument('--no-uppmax', action='store_true',
            default = False, help=("if specified the validation-dir and the "
            "assemblies-dir is assumed to contain the assemblies "
            "(and not samples) -- this is useful for large multi-library "
            "projects"))
    parser.add_argument('--sample-name', type=str,
            help=("It must be specifed when --no-uppmax is present, in this "
            "case you need to tell the porgram under which iouput name the "
            "validation has been saved (in the validation yaml file)"))
    args = parser.parse_args()
    main(args)

