import sys, os, yaml, glob
import subprocess
import argparse



def main(args):
    workingDir = os.getcwd()
    samples_data_dir = os.path.abspath(args.sample_data_dir)
    latex_document = _latexHeader(args.project_name)
    
    for sample_dir_name in [dir for dir in os.listdir(samples_data_dir) if os.path.isdir(os.path.join(samples_data_dir, dir))]:
        SampleQC_folder = os.path.join(samples_data_dir, sample_dir_name)
        Sample_fastqc   = os.path.join(SampleQC_folder, "fastqc")
        PE1_BasicsStatistics    = ""
        PE1_per_base_quality    = ""
        PE1_length_distribution = ""
        PE2_per_base_quality    = ""
        PE2_length_distribution = ""
        PE2_BasicsStatistics    = ""
        SE_per_base_quality     = ""
        SE_length_distribution  = ""
        SE_BasicsStatistics    = ""
        kmer_plot               = ""
        
        latex_document = _new_section(latex_document, sample_dir_name)
        
        #process the data
        
        
        for Sample_fastqc_dir in [fastqc_dir for fastqc_dir in os.listdir(Sample_fastqc) if os.path.isdir(os.path.join(Sample_fastqc, fastqc_dir))]:
            if "PE_1" in Sample_fastqc_dir:
                PE1_BasicsStatistics    = _getBasicStatsFromFastqc(os.path.join(Sample_fastqc, Sample_fastqc_dir, "fastqc_data.txt"), "pair 1")
                PE1_per_base_quality    = os.path.join(Sample_fastqc, Sample_fastqc_dir, "Images", "per_base_quality.png")
                PE1_length_distribution = os.path.join(Sample_fastqc, Sample_fastqc_dir, "Images", "sequence_length_distribution.png")
            elif "PE_2" in Sample_fastqc_dir:
                PE2_BasicsStatistics    = _getBasicStatsFromFastqc(os.path.join(Sample_fastqc, Sample_fastqc_dir, "fastqc_data.txt"), "pair 2")
                PE2_per_base_quality    = os.path.join(Sample_fastqc, Sample_fastqc_dir, "Images", "per_base_quality.png")
                PE2_length_distribution = os.path.join(Sample_fastqc, Sample_fastqc_dir, "Images", "sequence_length_distribution.png")
            elif "SE" in Sample_fastqc_dir:
                SE_BasicsStatistics    = _getBasicStatsFromFastqc(os.path.join(Sample_fastqc, Sample_fastqc_dir, "fastqc_data.txt"), "single")
                SE_per_base_quality = os.path.join(Sample_fastqc, Sample_fastqc_dir, "Images", "per_base_quality.png")
                SE_length_distribution = os.path.join(Sample_fastqc, Sample_fastqc_dir, "Images", "sequence_length_distribution.png")
            else:
                print "folder {} is not supposed to be here!!!!".format(Sample_fastqc_dir)
    
    
    
        latex_document = _insert_stat_table(latex_document, "\n".join([PE1_BasicsStatistics,PE2_BasicsStatistics, SE_BasicsStatistics]))
    
    
        Sample_kmer_dir = os.path.join(SampleQC_folder, "abyss_kmer")
        if os.path.isfile(os.path.join(Sample_kmer_dir, "kmer_coverage.png")):
            kmer_plot = os.path.join(Sample_kmer_dir, "kmer_coverage.png")
        else:
            print "folder {} does not contain kmer_coverage.png file: have you run the quality control pipeline?".format(Sample_kmer_dir)

        if SE_per_base_quality:
            latex_document = _new_three_figure(latex_document, PE2_per_base_quality, PE2_per_base_quality, SE_per_base_quality, "fastqc: per base quality distribution for read1, read2, and single ended reads")
        else:
            latex_document = _new_double_figure(latex_document, PE2_per_base_quality, PE2_per_base_quality, "fastqc: per base quality distribution for read1 and read2")

        if SE_per_base_quality:
            latex_document = _new_three_figure(latex_document, PE1_length_distribution, PE2_length_distribution, SE_length_distribution, "fastqc:  Sequence length distribution for read1, read2, and  single ended reads:")
        else:
            latex_document = _new_double_figure(latex_document, PE1_length_distribution, PE2_length_distribution, "fastqc: Sequence length distribution for read1 and read2")
    
        if kmer_plot:
            latex_document = _new_figure(latex_document, kmer_plot, "kmer plot showing the kmer distribution (computed with ABYSS-P)")
    
    
    latex_document = _latexFooter(latex_document)
    with open("{0}.tex".format(args.output),'w') as f:
        f.write(latex_document)
    command = ["pdflatex",  "{0}.tex".format(args.output)]
    latex_stdOut = open("latex.stdOut", "a")
    latex_stdErr = open("latex.stdErr", "a")
    subprocess.call(command, stdout=latex_stdOut , stderr=latex_stdErr)


def _latexHeader(projectName):
    LaTeX_header  = "\\documentclass{article}\n"
    LaTeX_header += "\\usepackage{graphicx}\n"
    LaTeX_header += "\\usepackage{caption}\n"
    LaTeX_header += "\\begin{document}\n"
    projectName = projectName.replace("_", "-")
    
    LaTeX_header += "\\title{{Quality report for {0} de novo project}}\n".format(projectName)
    LaTeX_header += "\\author{Francesco Vezzi}\n"
    LaTeX_header += "\\maketitle\n"
    LaTeX_header += "\\begin{abstract}\n"
    LaTeX_header += "This document shows general informations about the sequenced libraries that need to be used for de novo assembly\n"
    LaTeX_header += "\\end{abstract}\n"

    LaTeX_header += "\\section{Introduction}\n"
    LaTeX_header += "For each sample the following information is provided: \
    \\begin{itemize} \
    \\item Table showing number of reads, read length, and \%GC content (N.B. in case of trimming and/or merging single-ended reads will be present) \
    \\item Per base quality plot (generated by fastqc) \
    \\item Read length distribution (generated by fastqc) \
    \\item k-mer distribution (generated by abyss-pe) \
    \\end{itemize}\n"
    return LaTeX_header

def _new_section(latex_document, title):
    if "_" in title:
        title = "-".join(title.split("_"))
    latex_document += "\\subsection{"+ title + "}\n\n"
    return latex_document


def _new_figure(latex_document, picture, caption):
    figure  = "\n"
    figure += "\\begin{center}\n"
#    figure += "\\begin{figure}[H]\n"
#    figure += " \\centering\n"
    figure += " \\includegraphics[width=0.8\\textwidth]{" + picture + "}\n"
    figure += " \\captionof{figure}{" + caption + "}\n"
#    figure += " \\caption{" + caption + "}\n"
#    figure += "\end{figure}\n\n"
    figure += "\\end{center}\n"
    latex_document += figure
    return latex_document



def _new_double_figure(latex_document, picture1, picture2, caption):
    figure  = "\n"
    figure += "\\begin{center}\n"
#    figure += "\\begin{figure}[H]\n"
#    figure += " \\centering\n"
    figure += " \\includegraphics[width=0.48\\textwidth]{" + picture1 + "}\n"
    figure += " \\includegraphics[width=0.48\\textwidth]{" + picture2 + "}\n"
    figure += " \\captionof{figure}{" + caption + "}\n"
#    figure += " \\caption{" + caption + "}\n"
#    figure += "\end{figure}\n\n"
    figure += "\\end{center}\n"
    latex_document += figure
    return latex_document


def _new_three_figure(latex_document, picture1, picture2, picture3, caption):
    figure  = "\n"
    figure += "\\begin{center}\n"
    figure += " \\includegraphics[width=0.48\\textwidth]{" + picture1 + "}\n"
    figure += " \\includegraphics[width=0.48\\textwidth]{" + picture2 + "}\n"
    figure += " \\includegraphics[width=0.48\\textwidth]{" + picture3 + "}\n"
    figure += " \\captionof{figure}{" + caption + "}\n"
    figure += "\\end{center}\n"
    latex_document += figure
    return latex_document


def _latexFooter(latex_document):
    latex_document += "\\end{document}\n"
    return latex_document

def _getBasicStatsFromFastqc(fastqc_dataFile, line_header):
    fastqc_dataFileHandle = open(fastqc_dataFile, 'r')
    currentline = fastqc_dataFileHandle.readline().rstrip()
    while ">>Basic Statistics" in currentline:
        currentline = fastqc_dataFileHandle.readline().rstrip()
    currentline = fastqc_dataFileHandle.readline().rstrip() # >>Basic Statistics
    currentline = fastqc_dataFileHandle.readline().rstrip() # #Measure        Value
    currentline = fastqc_dataFileHandle.readline().rstrip() # Filename        lib2_PE_1.fastq.gz
    currentline = fastqc_dataFileHandle.readline().rstrip() # File type       Conventional base calls
    currentline = fastqc_dataFileHandle.readline().rstrip() # Encoding        Sanger / Illumina 1.9
    seqNumber   = fastqc_dataFileHandle.readline().rstrip().split("\t")[1]   # Total Sequences
    currentline = fastqc_dataFileHandle.readline().rstrip() # Filtered Sequences
    seqLength   = fastqc_dataFileHandle.readline().rstrip().split("\t")[1]   # Sequence length
    GCcontent   = fastqc_dataFileHandle.readline().rstrip().split("\t")[1]  # %GC
    return "{0} & {1} & {2} & {3} \\\\".format(line_header, seqNumber, seqLength, GCcontent)


def _insert_stat_table(latex_document, content):
    tabular = "\n"
    tabular += "\\begin{center}\n"
    tabular += "\\begin{tabular}{c|lll}\n"
    tabular += " type & \\#reads & length & \\%GC \\\\\\hline \n"
    tabular += content
    tabular += "\n\\end{tabular}\n"
    tabular += "\\end{center}\n"
    latex_document += tabular
    #print tabular
    return latex_document

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-data-dir', help="full path to directory containing one folder per sample. Each sample contaains only one library (i.e., one PE lib)", type=str)
    parser.add_argument('--project-name', help="project name", type=str)
    parser.add_argument('--output', help="name of the report", type=str)
    parser.add_argument('--global-config', help="global configuration with")
    
    args = parser.parse_args()

    main(args)
