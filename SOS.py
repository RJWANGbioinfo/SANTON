#!/usr/bin/env python3

import os, sys
import glob
import argparse
import pandas as pd
import ray
from ray.util.multiprocessing import Pool
from modules.Oligo import *
from modules.Visualization import manhattanPlotOTscoresLoop

    
### Default constants and sequences  ###
constant1="GACTAGACACTGCTATCACACTCTCTCA"
constant2="AGACGTTCTCACAGCAATTCGTACAGTCGACGTCGATTCGTGT"
constant3="TTGACATTCTGCAATTA"
constant4="AGTATGTATGCTTCGCGCAGTGCGACTTCGCAGCGCATCACTTCA"
constant5="AGAGCTGCGAGTCTTACAGCATTGCA"
illumina_adapter="AGATCGGAAGAGC"
umi_length=11
wkdir=os.getcwd()


def main(job_number):
    # List merged fastQ files
    fastq_files = sorted(glob.glob(os.path.join(input_folder, '*.fastq.gzip')))
    if not fastq_files:
        print("No fastq.gzip files found in the input folder.")
        sys.exit(1)
    
    spike_name = os.path.basename(fastq_files[0]).replace(".extendedFrags.fastq.gzip", "").replace(".fastq.gzip", "")

    # Extract sample name from analysis file
    df_analysis = pd.read_csv(analysis_file, sep='\t')
    df_analysis = df_analysis.assign(Merge_Fq=fastq_files[0])
    df_sample = df_analysis[df_analysis['Fastq_R1'].str.contains(spike_name, regex=False)]
    
    if df_sample.empty:
        print("No matching samples found in the analysis file.")
        sys.exit(1)
    
    sample_name = df_sample['Sample_name'].unique()[0]
    print(f"Processing sample: {sample_name}")

    temp_dir = os.path.join(output_folder, f"Tmp_{sample_name}")
    out_dir = os.path.join(output_folder, sample_name)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    temp_analysis_file = os.path.join(temp_dir, os.path.basename(analysis_file))
    df_analysis.to_csv(temp_analysis_file, index=False, sep='\t')

    # Adjust job_numbers when number of off-targets is lower than number of CPUs
    job_number = min(job_number, len(df_sample))

    # Prepare arguments for parallel processing
    Reverse, Type, runnerlist = MPIwarpAnalysisFile(
        temp_analysis_file, sample_name, illumina_adapter,
        constant1, constant2, constant3, constant4, constant5, 
        umi_length, temp_dir, job_number
    )
    
    print(f"{len(runnerlist)} jobs initialized...")

    # Parallel processing
    with Pool(processes=job_number) as pool:
        data_bulk = pool.map(quantifyReadsUMIs, runnerlist)
    ray.shutdown()

    # Combine extracted data and summarize
    n = data_bulk[0][0]
    c1, r1 = 0, 0
    df_constant123_combined, df_constant456_combined = pd.DataFrame(), pd.DataFrame()

    for data in data_bulk:
        df_constant123, df_constant456, count123, count456 = data[3], data[4], data[1], data[2]
        c1 += count123
        r1 += count456
        df_constant123_combined = pd.concat([df_constant123_combined, df_constant123])
        df_constant456_combined = pd.concat([df_constant456_combined, df_constant456])

    # Calculate unique UMI counts
    df_constant123_umi_info, df_constant123_umi_count = UMIunique(df_constant123_combined)
    df_constant456_umi_info, df_constant456_umi_count = UMIunique(df_constant456_combined)

    # Summarize and export results
    summarySOS(n, c1, r1, Reverse, out_dir, sample_name, 
                  df_constant123_combined, df_constant456_combined, 
                  df_constant123_umi_info, df_constant456_umi_info, 
                  df_constant456_umi_count)

    # Annotate sequences using parallel processing
    runner_list_constant123 = preRunnerAnnotation(df_constant123_umi_count, 100, out_dir, Reverse, Type, 
                                            temp_analysis_file, "c1umicount.tsv", sample_name, temp_dir)
    runner_list_constant456 = preRunnerAnnotation(df_constant456_umi_count, 100, out_dir, Reverse, Type, 
                                            temp_analysis_file, "r1umicount.tsv", sample_name, temp_dir)
    
    with Pool(processes=job_number) as pool:
        annotation_constant123 = pool.map(MPIAnnotateEditStatusC123, runner_list_constant123)
        annotation_constant456 = pool.map(MPIAnnotateEditStatusC456, runner_list_constant456)

    ray.shutdown()

    # Combine annotated data
    df_constant123_umi_count = pd.concat(annotation_constant123)
    df_constant456_umi_count = pd.concat(annotation_constant456)

    df_constant456_umi_count=calculateUMIPercentage(df_constant456_umi_count)
    df_constant123_umi_count=calculateUMIPercentage(df_constant123_umi_count)
    
    # Save dataframe to local file 
    df_constant123_umi_count = df_constant123_umi_count.iloc[:, 1:]
    df_constant456_umi_count = df_constant456_umi_count.iloc[:, 1:]
    df_constant123_umi_count.columns = ['Template_ID', 'Sequence',  'UMI_count', 'Index', 'Chunk', 'Editing_status', 'UMI_percentage']
    df_constant456_umi_count.columns = ['Template_ID', 'Sequence',  'UMI_count', 'Index', 'Chunk', 'Editing_status', 'UMI_percentage']
    pam_side_umi_file = os.path.join(out_dir, f"{sample_name}.allele.PAMSide.UMI_unique.UMIcount.txt") 
    proto_side_umi_file = os.path.join(out_dir, f"{sample_name}.allele.ProtoSide.UMI_unique.UMIcount.txt") 
    if Reverse=="Yes":
        df_constant123_umi_count.to_csv(pam_side_umi_file, index=False, sep='\t')
        df_constant456_umi_count.to_csv(proto_side_umi_file, index=False, sep='\t')
    else:
        df_constant123_umi_count.to_csv(proto_side_umi_file, index=False, sep='\t')
        df_constant456_umi_count.to_csv(pam_side_umi_file, index=False, sep='\t')

    # Calulate score and plot
    calculateOTScore(df_analysis, out_dir, sample_name, Type)
    manhattanPlotOTscoresLoop(analysis_file, out_dir, score_cutoff)

    # Clean and create temp_dir
    os.system("rm -rf " + temp_dir)


if __name__ == "__main__":

    # compile arguments
    parser = argparse.ArgumentParser(description="Analysis merged NGS reads from synthetic oligonucleotide-based sequencing.")
    parser.add_argument("-a", "--analysis_file", 
                        type=str, action="store", default="",
                        help="Required, define local path to an analysis file")
    parser.add_argument("-i", "--input_folder", 
                        type=str, action="store", default=wkdir,
                        help="Required, define local path to a folder containing input fastQ files.")
    parser.add_argument("-o", "--output_folder", 
                        type=str, action="store", default=wkdir, 
                        help="Optional, define local path to a folder holding output files. Default is current directory.")
    parser.add_argument('-n', "--processor_number",
                        type=int, action="store", default="1",
                        help="Optional, number of processor to be used for parallel computing")
    parser.add_argument('-s', "--score_cutoff",
                        type=float, action="store", default="0.1",
                        help="Optional, mininal score for an off-target region to be included in manhattan plot")

    # Parse user arguments
    args = parser.parse_args()
    analysis_file = args.analysis_file
    processor_number = args.processor_number
    input_folder = args.input_folder
    output_folder = args.output_folder
    score_cutoff = args.score_cutoff

    # Check user inputs
    if len(sys.argv)==1: 
        parser.print_help(sys.stderr); 
        sys.exit(1)

    # Data analysis
    main(processor_number)

