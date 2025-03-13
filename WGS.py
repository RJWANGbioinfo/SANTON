#!/usr/bin/env python3

import os, sys
import ast
import argparse
import numpy as np
import pandas as pd
from IPython.display import display
import multiprocessing

from modules.BamFile import wgs_analysis, isRedFlag
from modules.Alignment import *
from modules.Visualization import *
from modules.Bed import *
from modules.Utility import *



def checkARead(read):
    """
    Filter pysam read based on mapping quality and read flags.

    This function is to be passed to pysam function 'count_coverage' as callback 
    and accepts cutoff_mapping_quality variable as global variable. Because 
    global variable change has no effect on function executed in Ray thread, this
    function has to stay in main script so that user specified cutoff_mapping_quality value
    can be used in multitreading process.

    Returns 
        boolean: True or False
    """
    # cutoff_mapping_quality has to be defined in outer closure as global variable
    # because 
    if read.mapping_quality < cutoff_mapping_quality or isRedFlag(read.flag):
        return False
    else:
        return True
    

def processBamFiles():
    """
    Process one or more input BAM file(s)
    """

    control_sample_outfile = getOutputFilePaths(output_file)

    # scanning cut sites in input bam file(s)
    site_pair_collections = []
    for input_bam_file in [bam_file_edited, bam_file_control]:
        site_pairs_from_one_sample = pd.DataFrame()
        if input_bam_file:
            site_pairs_from_one_sample = wgs_analysis(
                input_bam_file,
                cutoff_supporting_reads_forward, 
                cutoff_supporting_reads_reverse, 
                cutoff_mapping_quality, 
                gamma,
                cutoff_ratio, 
                cutoff_coverage, 
                cutoff_score,
                output_file,
                guide_rna, 
                pam_location,
                enzyme,
                reference_genome,
                checkARead,
                processor_number,
                )
        site_pair_collections.append(site_pairs_from_one_sample)
        
    site_pairs_in_edited, site_pairs_in_control = site_pair_collections
    print(f'{site_pairs_in_edited.shape[0]} sites are found in bam file {bam_file_edited}')
    if bam_file_control:
        print(f'{site_pairs_in_control.shape[0]} sites are found in bam file {bam_file_control}')
    print('\n')

    if site_pairs_in_control.shape[0]:
        site_pairs_in_edited = removeSitesFoundInControl(site_pairs_in_edited, site_pairs_in_control)
    


    # report and save final list of site pairs
    print(f'{site_pairs_in_edited.shape[0]} coordinate pairs (forward & reverse) meeting all requirements')
    if output_file:
        if enzyme.upper() in ['CPF1', 'ABE', 'CBE']:
            # write raw predictions and cluster
            site_pairs_in_edited_cluster = add_cluster_annotations(site_pairs_in_edited)
            save_and_plot_result(site_pairs_in_edited_cluster, f'{remove_file_suffix(output_file)}.raw_predictions.csv')
            # write representative site per cluster
            site_pairs_in_edited_nonredundant = get_highest_score_per_cluster(site_pairs_in_edited_cluster)
            save_and_plot_result(site_pairs_in_edited_nonredundant, output_file)
        else:
            save_and_plot_result(site_pairs_in_edited, output_file)
    if site_pairs_in_edited.shape[0]:
        display(site_pairs_in_edited.head())

    if bam_file_control:
        save_and_plot_result(site_pairs_in_control, control_sample_outfile)


if __name__ == "__main__":

    # set up argument parser
    text= "This is a tool to analyzed WGD NGS data derived from " \
          "samples treated by enzyme such as Cas9, Cas12a and ABE/CBE for off-target nomination."
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument('-be', '--bam_file_edited',
                        type=str, default='',
                        help=": a bam file derived from edited sample")
    parser.add_argument("-q", "--cutoff_mapping_quality",
                        type=int, default=1,
                        help=": minimal mapping quality for a read to be considered in analysis")
    parser.add_argument("-f", "--cutoff_supporting_reads_forward",
                        type=int, default=5,
                        help=": minimal count of supporting reads aligned forwardly to reference")
    parser.add_argument("-r", "--cutoff_supporting_reads_reverse",
                        type=int, default=5,
                        help=": minimal count of supporting reads aligned reverse-complementarily to reference")
    parser.add_argument("-d", "--cutoff_coverage",
                        type=int, default=10,
                        help=": minimal count of reads that mapped to a position in reference regardless of mapping orientation")
    parser.add_argument("-R", "--cutoff_ratio",
                        type=float, default=0.02,
                        help=": minimal fraction of reads supporting edit at a position in reference")
    parser.add_argument("-G", "--gamma",
                        default=0,
                        help=": number of base pair between edit: 0 for blunt-end cut; >0 for 5'overhang cut; <0 for 3'overhang")
    parser.add_argument("-s", "--cutoff_score",
                        type=float, default=2.5,
                        help=": minimal score for a cut site to be included in final report")
    parser.add_argument("-e", "--enzyme",
                        type=str, default='',
                        help=": enzyme used genome editing. Supported options include 'Cas9', 'ABE', 'CBE', and 'Cpf1'.")
    parser.add_argument("-g", "--reference_genome",
                        type=str, default='',
                        help=": path to reference genome (optional)")
    parser.add_argument("-grna", "--guide_rna",
                        type=str, default='',
                        help=": guide RNA sequence including spacer and PAM (optional)")
    parser.add_argument("-pam", "--pam_location",
                        type=str, default='',
                        help=": location of pam relative to spacer, either 'left' or 'right' (optional)")
    parser.add_argument("-n", "--processor_number",
                        type=int, default= multiprocessing.cpu_count() - 1,
                        help=": number of processor to be used for parallel computing (optional)")
    parser.add_argument("-o", "--output_file",
                        type=str, default='',
                        help=": path to a *.csv output file (optional)")
    parser.add_argument('-bc', '--bam_file_control',
                        type=str, default='',
                        help=": a bam file derived from control sample")
    

    # get parser ready
    args=parser.parse_args()

    if len(sys.argv)==1: 
        parser.print_help(sys.stderr); 
        sys.exit(1)
    elif not args.bam_file_edited:
        print(f"bam_file_edited is mandatory\n")
        parser.print_help(sys.stderr); 
        sys.exit(1)


    # population user input arguments
    bam_file_edited = args.bam_file_edited
    bam_file_control = args.bam_file_control
    cutoff_supporting_reads_forward = args.cutoff_supporting_reads_forward
    cutoff_supporting_reads_reverse = args.cutoff_supporting_reads_reverse
    cutoff_coverage = args.cutoff_coverage
    cutoff_mapping_quality = args.cutoff_mapping_quality
    cutoff_ratio = args.cutoff_ratio
    cutoff_score = args.cutoff_score
    gamma = args.gamma
    reference_genome = args.reference_genome
    guide_rna = args.guide_rna
    processor_number = args.processor_number
    output_file = args.output_file
    enzyme = args.enzyme
    pam_location = args.pam_location


    # check number of processor 
    if processor_number < 1: processor_number = 1


    # check gamma variable
    enzymes = ['CAS9', 'ABE', 'CBE', 'CPF1']
    if enzyme and enzyme.upper() not in enzymes:
        print(f'Invalid enzyme type. Supported enzyme list is {enzymes}.')
        exit(1)
    if enzyme.upper() == 'CAS9':
        gamma = [0]
    elif enzyme.upper() in ['ABE', 'CBE']:
        gamma = list(np.arange(1,20))
    elif enzyme.upper() in ['CPF1']:
        gamma = list(np.arange(-10,11))

    if isinstance(gamma, str):
        try:
            gamma = ast.literal_eval(gamma)
        except ValueError:
            print(f'Invalid gamma "{gamma}". Gamma must be an integer, a list of integer, or a string expressing a list of integers')
            sys.exit(1)
    if not (isinstance(gamma, int) or isinstance(gamma, list)):
        print(f'Invalid gamma "{gamma}". Gamma must be an integer, a list of integer, or a string expressing a list of integers')
        sys.exit(1)
    if isinstance(gamma, int):
        gamma = [gamma]



    # print to console all paremeters
    user_argument_string = f"input paremeters:\n-q {cutoff_mapping_quality} -f {cutoff_supporting_reads_forward} -r {cutoff_supporting_reads_reverse} -d {cutoff_coverage} -R {cutoff_ratio} -s {cutoff_score} -n {processor_number}"
    if enzyme:
        user_argument_string = f"{user_argument_string}\n-e {enzyme}"
    user_argument_string = f"{user_argument_string} -G {gamma}"
    if guide_rna:
        user_argument_string = f"{user_argument_string}\n-grna {guide_rna}"
    if pam_location:
        user_argument_string = f"{user_argument_string} -pam_location {pam_location}"
    if bam_file_edited:
        user_argument_string = f"{user_argument_string}\n-be {bam_file_edited}"
    if bam_file_control:
        user_argument_string = f"{user_argument_string}\n-bc {bam_file_control}"
    if reference_genome:
        user_argument_string = f"{user_argument_string}\n-g {reference_genome}"
    if output_file:
        user_argument_string = f"{user_argument_string}\n-o {output_file}"
    print(user_argument_string, '\n')
    

    # sanity check input paths/files
    input_file_errors = []
    input_file_errors = input_file_errors + sanityCheckInputFile(bam_file_edited, 'bam')
    input_file_errors = input_file_errors + sanityCheckInputFile(bam_file_control, 'bam')
    input_file_errors = input_file_errors + sanityCheckInputFile(reference_genome, 'fasta')
    if len(input_file_errors):
        print('\n'.join(input_file_errors))
        exit(1)


    # search cut sites in input bam files, and go through downstream steps
    processBamFiles()
