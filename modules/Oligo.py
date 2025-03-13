import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import gzip

from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    from Bio import pairwise2



def UMIunique(df):
    """
    Consolidate unique UMIs and calculate read counts and sequence quality.

    Args:
        df (DataFrame): Input DataFrame with UMI, template, read, and quality information.

    Returns:
        Tuple[DataFrame, DataFrame]: 
            - df_consolidated_umi: DataFrame with unique UMI sequences, sequence lengths, quality scores, and read counts.
            - df_umi_count_per_template: DataFrame with UMI counts aggregated per template and sequence.
    """
    # Calculate read count per UMI and assign sequence length
    df_read_count_per_umi = df.groupby(['Template_ID', 'UMI', "Reads"]).size().reset_index()
    df_read_count_per_umi.columns = ['Template_ID', 'UMI', 'Sequence', 'Count']
    df_read_count_per_umi = df_read_count_per_umi.assign(seqlen=df_read_count_per_umi['Sequence'].str.len())
    df_quality = df.groupby(['Template_ID', 'UMI', "Reads"])['AvgQual'].max().reset_index()
    df_quality.columns = ['Template_ID', 'UMI', 'Sequence', 'AvgQual']
    df_read_count_per_umi = pd.merge(df_read_count_per_umi, df_quality, on=['Template_ID', 'UMI', "Sequence"], how="left")
    df_read_count_per_umi = df_read_count_per_umi.sort_values(by=['Template_ID', 'UMI', 'Sequence', 'Count', 'seqlen', 'AvgQual'], ascending=False)

    # Consolidate by retaining the highest-ranked sequence for each UMI
    df_consolidated_umi = df_read_count_per_umi.drop_duplicates(subset=['Template_ID', 'UMI'])
    df_consolidated_umi = df_consolidated_umi[['Template_ID', 'UMI', 'Sequence', 'seqlen', 'AvgQual']]
    df_reads_count = df.groupby(['Template_ID', 'UMI']).size().reset_index()
    df_reads_count.columns = ['Template_ID', 'UMI', 'ReadsCount']
    df_consolidated_umi = pd.merge(df_consolidated_umi, df_reads_count, on=['Template_ID', 'UMI'], how='left')
    df_consolidated_umi = df_consolidated_umi.sort_values(by=['ReadsCount'], ascending=False)

    # Calculate UMI counts per template
    df_umi_count_per_tempalte = df_consolidated_umi.groupby(['Template_ID', "Sequence"]).size().reset_index()
    df_umi_count_per_tempalte.columns = ['Template_ID', "Sequence", "UMI_Count"]
    df_umi_count_per_tempalte = df_umi_count_per_tempalte.sort_values(by=['UMI_Count'], ascending=False)

    return df_consolidated_umi, df_umi_count_per_tempalte



def read_align_reference(read, reference):
    """
    Perform local alignment between a read and a reference sequence.

    Args:
        read (str): The off-target read sequence.
        reference (str): The template or reference sequence.

    Returns:
        list: Alignment object(s) representing the optimal local alignment.
    """
    return pairwise2.align.localms(read, reference, 2, -1, -3, -3)



def checkEditedStatusC123(off_target_seq, template, expcut=30, type='Cas9'):
    """
    Align off-target sequence adjacent to constant regions (1, 2, 3) against its template 
    and annotate the edit status.

    Args:
        off_target_seq (str): Sequence adjacent to constant regions to be aligned.
        template (str): Reference template sequence.
        exp_cut (int): Expected cut site (1-based, adjusted to 0-based internally).

    Returns:
        tuple: 
            - edite_status (str): 'Edited' or 'NotValidEdited' based on alignment.
            - query (str): Aligned off-target sequence.
            - reference (str): Aligned reference sequence.
    """
    # Adjust expected cut site to 0-based indexing
    expcut -= 1  
    
    # Default edit status
    edit_status = "NotValidEdited"
    type = type.upper()
    
    # Perform sequence alignment
    alignments = read_align_reference(off_target_seq, template)
    
    # Return default if no alignment is found
    if not alignments:
        return edit_status, "", ""
    
    # Extract the best alignment
    aln = alignments[0]
    query, reference = aln[0], aln[1]
    
    # Initialize tracking variables
    change_detected = False
    difference = 0
    cutsite = -1
    
    # Iterate in reverse to find the cut site and mismatches
    for i in reversed(range(len(query))):
        qs = query[i]  # Character from query
        rs = reference[i]  # Character from reference
        
        if qs != "-" or change_detected:
            if not change_detected:
                cutsite = i  # Record cut site
            change_detected = True
            
            # Count mismatches
            if qs != rs:
                difference += 1
    
    # Set cutoffs based on enzyme type
    cut_window_protoside = cut_window_pamside = 3 
    max_difference = 2
    if type == 'CPF1':
        cut_window_protoside = 6
        cut_window_pamside = 10
    if type == 'ABE':
        max_difference = 4

    # Determine edit status based on cutsite proximity and mismatch tolerance,
    # and A>G conversion for ABE case
    if (expcut - cut_window_pamside) <= cutsite <= (expcut + cut_window_protoside) and difference <= max_difference:
        if type == 'ABE':
            for p in range(cutsite):
                if query[p] == "G" and reference[p] == "A":
                    edit_status="Edited"
                    break
        else:
            edit_status = "Edited"
    
    return edit_status, query, reference



def checkEditedStatusC456(off_target_seq, template, expcut=30, type='Cas9'):
    """
    Align off-target sequence adjacent to constant regions (4, 5, 6) against its template
    and annotate the edit status.

    Args:
        off_target_seq (str): PAM-side sequence to be aligned.
        template (str): Reference template sequence.
        expcut (int): Expected cut site (1-based, adjusted to 0-based internally).

    Returns:
        tuple:
            - edit_status (str): 'Edited' or 'NotValidEdited' based on alignment.
            - query (str): Aligned PAM-side sequence.
            - reference (str): Aligned reference sequence.
    """

    # Adjust expected cut site to 0-based indexing
    expcut -= 1  
    
    # Default edit status
    edit_status = "NotValidEdited"
    type = type.upper()
    
    # Perform local alignment
    alignments = read_align_reference(off_target_seq, template)
    
    # Return default if no alignment is found
    if not alignments:
        return edit_status, "", ""
    
    # Extract best alignment
    aln = alignments[0]
    query, reference = aln[0], aln[1]
    
    # Initialize variables to track changes and mismatches
    change_detected = False
    difference = 0
    cut_site = -1
    
    # Analyze alignment to detect mismatches and cut sites
    for i, (qs, rs) in enumerate(zip(query, reference)):
        if qs != "-" or change_detected:
            if not change_detected:
                cut_site = i  # Record cut site at the first non-gap character
            change_detected = True
            
            # Count mismatches
            if qs != rs:
                difference += 1
    
    # Set cutoffs based on enzyme type
    cut_window_protoside = cut_window_pamside = 3 
    max_difference = 2
    if type == 'CPF1':
        cut_window_protoside = 6
        cut_window_pamside = 10
    elif type == 'CBE':
        cut_window_protoside = cut_window_pamside = 6

    # Determine edit status based on cut site proximity and mismatch tolerance,
    # and 'C' presence around alignment edge for CBE case
    if (expcut - cut_window_pamside) <= cut_site <= (expcut + cut_window_protoside) and difference <= max_difference:
        if type == 'CBE':
            if off_target_seq[0]=="C" or template[cut_site-1]=="C":
                edit_status="Edited"
        else:
            edit_status = "Edited"
    
    return edit_status, query, reference



def trimAdaptor(illumina_adapter, sequence):
    """
    Trim Illumina adapter sequences from the input DNA sequence.

    Args:
        illumina_adapter (str): The Illumina adapter sequence to be trimmed.
        sequence (str): The input DNA sequence from which the adapter should be removed.

    Returns:
        str: The trimmed sequence without the adapter or reverse complement.
    """
    if illumina_adapter in sequence:
        sequence = sequence.split(illumina_adapter, 1)[0]        
    elif str(Seq(illumina_adapter).reverse_complement()) in sequence:
        sequence = sequence.split(str(Seq(illumina_adapter).reverse_complement()), 1)[0]
    return sequence



def calculateSingleSideScore(suffix_name, df_analysis, out_dir, sample_name):
    """
    Calculate proto- or pam-side scores based on UMI counts and generate a DataFrame with scores.

    Args:
        suffix_name (str): Filename suffix for allele data.
        df_analysis (DataFrame): DataFrame containing analysis and target information.
        out_dir (str): Directory containing allele files.
        sample_name (str): Name of the sample to process.

    Returns:
        DataFrame: DataFrame with sample name, TemplateID, and calculated side-specific scores.
    """
    on = df_analysis[(df_analysis.Sample_name == sample_name) & (df_analysis.On_target == "Yes")].Template_ID.values[0]
    print("Reading allele files")
    df_allele = pd.read_csv(f"{out_dir}/{sample_name}.{suffix_name}", sep='\t')
    df_allele.columns = ['Template_ID', 'Sequence', 'UMI_Count', 'inx', 'chunk', 'EditingStatus', 'UMI_Pct']

    # Sum On-target Edit UMI
    print("Calculating scores")
    on_sum = df_allele[(df_allele.Template_ID == on) & (df_allele.EditingStatus == "Edited")].UMI_Count.sum()
    df_allele = df_allele[(df_allele.Template_ID != on) & (df_allele.EditingStatus == "Edited")]
    df_allele_summary = df_allele.groupby(['Template_ID'])['UMI_Count'].sum().reset_index()
    df_allele_summary = df_allele_summary.assign(Sample_name=sample_name)
    df_allele_summary = df_allele_summary.assign(score=df_allele_summary.UMI_Count / on_sum)
    df_score = df_allele_summary[['Sample_name', 'Template_ID', 'score']]
    key_name = suffix_name.replace("allele.", "").replace(".UMI_unique.UMIcount.txt", "")
    df_score.columns = ['Sample_name', 'Template_ID', key_name + "_score"]
    return df_score



def calculateOTScore(dfAnalysis, out_dir, sample_name, Type):
    """
    Calculate the OT score by combining proto-side and PAM-side scores.

    Args:
        dfAnalysis (DataFrame): DataFrame containing analysis and target information.
        out_dir (str): Directory to save the output file.
        sample_name (str): Name of the sample to process.
        Type (str): Editing type, either "Cas9", "CBE", or "ABE".

    Returns:
        None: Saves the OT score as a tab-delimited file in the specified output directory.
    
    Notes:
        - For "Cas9" or "CBE", the OT score is the average of PAM-side and proto-side scores.
        - For "ABE", the proto-side score is directly used as the OT score.
    """
    df_pamside_score = calculateSingleSideScore("allele.PAMSide.UMI_unique.UMIcount.txt", dfAnalysis, out_dir, sample_name)
    df_protoside_score = calculateSingleSideScore("allele.ProtoSide.UMI_unique.UMIcount.txt", dfAnalysis, out_dir, sample_name)
    if Type=="Cas9" or Type=="CBE": 
        df_ot_score = pd.merge(df_pamside_score, df_protoside_score, on=["Sample_name","Template_ID"])
        df_ot_score["OT_score"] = df_ot_score[["PAMSide_score", "ProtoSide_score"]].mean(axis = 1)
    elif Type=="ABE":
        df_ot_score = df_protoside_score.copy()
        df_ot_score["OT_score"] = df_ot_score["ProtoSide_score"].copy()

    df_ot_score.columns = ['Sample_name', 'Template_ID', 'PAMside_score', 'Protoside_score', 'OT_score']    
    df_ot_score.to_csv(out_dir+"/"+sample_name+".OT_Score.txt", index=False, sep='\t')



def MPIwarpAnalysisFile(analysis_file, sample_name, illumina_adapter, constant1, constant2, constant3, constant4, constant5, umi_length, temp_dir, job):
    """ 
    Prepare array data for parallel fastQ data processing.

    Args:
        analysis_file (str): Path to the input analysis file (tab-delimited).
        sample_name (str): Sample name to filter analysis data.
        illumina_adapter (str): Illumina adapter sequence for trimming.
        constant1-5 (various): Constant regions used in downstream processing.
        umi_length (int): Length of UMI sequences.
        temp_dir (str): Directory for storing temporary files.
        job (int): Number of parallel jobs to divide the workload.

    Returns:
        tuple: 
            - reverse (str): Reverse sequence indicator for the sample.
            - sample_type (str): Editing type (e.g., Cas9, ABE).
            - runnerlist (list): List of parameters for parallel processing tasks.
    """
    # Read analysis file and filter for the specified sample
    df_analysis = pd.read_csv(analysis_file, sep='\t')
    df_sample = df_analysis[df_analysis['Sample_name'] == sample_name]
    
    # Extract metadata from the sample
    extension = df_sample['Extension'].iloc[0]
    reverse = df_sample['Reverse'].iloc[0]
    sample_type = df_sample['Type'].iloc[0]
    fastq_file = df_sample['Merge_Fq'].iloc[0]
    
    # Determine chunk size based on job count
    chunk_size = len(df_sample) / (job - 1) if job > 1 else len(df_sample) * 2
    
    # Assign chunk indices for parallel processing
    df_sample=df_sample.assign(inx=df_sample.index.map(int).copy())
    df_sample=df_sample.assign(chunk=round(df_sample['inx'].map(int)/int(chunk_size),0))
    df_sample['chunk']=df_sample['chunk'].map(int)

    # Save filtered sample to a temporary file
    sample_file = f"{temp_dir}/{sample_name}.SAMPLFILE.tsv"
    df_sample.to_csv(sample_file)
    
    # Generate runner list for each unique chunk
    runnerlist = [
        [
            chunk, fastq_file, constant1, constant2, constant3, constant4, constant5,
            umi_length, extension, illumina_adapter, sample_file, temp_dir, sample_name
        ]
        for chunk in df_sample['chunk'].unique()
    ]

    return reverse, sample_type, runnerlist



def findBetween(string, first, last):
    """
    Extract and return the substring between two specified substrings.

    Args:
        string (str): The input string to search within.
        first (str): The starting delimiter substring.
        last (str): The ending delimiter substring.

    Returns:
        str: The substring found between 'first' and 'last'. 
             Returns an empty string if the delimiters are not found.
    """
    try:
        start = string.index( first ) + len( first )
        end = string.index( last, start )
        return string[start:end]
    except ValueError:
        return ""



def quantifyReadsUMIs(runnerlist):
    """
    Process fastQ reads to extract UMI, template ID, sequence, and quality metrics.

    Args:
        runnerlist (list): A list of parameters necessary for processing, including:
            - chunk (int): Current chunk of data being processed.
            - fastq_file (str): Path to the input fastQ file.
            - constant1-constant5 (str): Constant regions used for sequence identification.
            - umi_length (int): Length of the UMI sequence.
            - extension (int): Length of the sequence extension to trim.
            - illumina_adapter (str): Adapter sequence for trimming.
            - sample_file (str): Path to the sample metadata file.
            - temp_dir (str): Directory to save output files.
            - sample_name (str): Name of the sample being processed.

    Returns:
        list: A list containing:
            - n (int): Total number of reads processed.
            - count123 (int): Number of proto-side reads passing UMI extraction.
            - count456 (int): Number of PAM-side reads passing UMI extraction.
            - df_umi123 (DataFrame): DataFrame of extracted proto-side UMI reads.
            - df_umi456 (DataFrame): DataFrame of extracted PAM-side UMI reads.
    """
    (chunk, fastq_file, constant1, constant2, constant3, constant4,
     constant5, umi_length, extension, illumina_adapter,
     sample_file, temp_dir, sample_name) = runnerlist
    
    # Load sample metadata
    df_sample = pd.read_csv(sample_file)
    df_sample = df_sample[df_sample.chunk == chunk]

    barcode_list = df_sample['Barcode'].tolist()
    df_sample['common_seq'] = constant2 + df_sample['Barcode'] + constant3
    df_sample['commonseq2'] = constant4 + df_sample['Barcode'] + constant5
    
    target_dict_23 = dict(zip(df_sample['common_seq'], df_sample['Template_ID']))
    target_dict_45 = dict(zip(df_sample['commonseq2'], df_sample['Template_ID']))

    # Initialize counters and lists
    n = 1
    count123, count456 = 0, 0
    read_list123, umi_list123, quality_list123, target_list123 = [], [], [], []
    read_list456, umi_list456, quality_list456, target_list456 = [], [], [], []

    # Read fastQ file
    read1file = gzip.open(fastq_file,"rt")
    fastq_reads = SeqIO.parse(read1file, "fastq")

    for read in fastq_reads:
        n += 1
        read_sequence = str(read.seq)
        read_sequence_rc = str(read.seq.reverse_complement())

        # Check if all constant regions are intact, skip if not
        have_intact_constants = any([
            all(c in read_sequence for c in [constant1, constant2, constant3]),
            all(c in read_sequence for c in [constant4, constant5]),
            all(c in read_sequence_rc for c in [constant1, constant2, constant3]),
            all(c in read_sequence_rc for c in [constant4, constant5])
        ])
        if not have_intact_constants:
            continue

        # Reverse complement if constants are detected in reverse
        if not (constant1 in read_sequence or constant4 in read_sequence):
            read = read.reverse_complement()
            read_sequence = read_sequence_rc

        # Extract barcodes
        barcode23 = findBetween(read_sequence, constant2, constant3)
        barcode45 = findBetween(read_sequence, constant4, constant5)

        # Extract annotation from reads encoding constant regions1,2,3
        if barcode23 and barcode23 in barcode_list and (read_sequence.startswith(constant1)):
            common_seq = constant2 + barcode23 + constant3
            index = read_sequence.find(common_seq)
            c1_umi = read_sequence[index - umi_length:index]
            
            if len(c1_umi)==umi_length:
                count123 += 1
                quality_off_target = read.letter_annotations["phred_quality"][index + len(common_seq):]
                c1_seq = read_sequence[index + len(common_seq):]
                c1_seq = trimAdaptor(illumina_adapter, c1_seq)
                c1_seq = c1_seq[:len(c1_seq) - extension] # added double clean up of seq
                read_list123.append(c1_seq)
                umi_list123.append(c1_umi)
                target_list123.append(target_dict_23[common_seq])
                
                if (len(c1_seq) > 0) and (type(quality_off_target) == list) and (len(quality_off_target) > 0):
                    quality_list123.append(sum(quality_off_target)/len(quality_off_target))
                else:
                    quality_list123.append(0)
        # Extract annotation from reads encoding constant regions4,5,6
        elif barcode45 and barcode45 in barcode_list:
            common_seq = constant4 + str(barcode45)+constant5
            index = read_sequence.find(common_seq)
            r1_umi = read_sequence[index + len(common_seq):index + len(common_seq) + umi_length]
            if len(r1_umi) == umi_length:
                count456 = count456 + 1
                quality_off_target = read.letter_annotations["phred_quality"][:index]
                r1_seq = read_sequence[:index]
                r1_seq = trimAdaptor(illumina_adapter, r1_seq)
                r1_seq = r1_seq[extension:] # added double clean up of seq
                read_list456.append(r1_seq)
                umi_list456.append(r1_umi)
                target_list456.append(target_dict_45[common_seq])
                    
                if (len(r1_seq) > 0) and (type(quality_off_target)==list) and (len(quality_off_target)>0):
                    avgqr1=sum(quality_off_target)/len(quality_off_target)
                    quality_list456.append(avgqr1)
                else:
                    quality_list456.append(0)

        if n%10000==0:
            print("-----Processing ",str(chunk)+" :"+str(n)+ f" reads-----")

    # Save extracted data
    def save_to_file(read_list, umi_list, quality_list, target_list, suffix):
        df_umi = pd.DataFrame({
            'Template_ID': target_list,
            'UMI': umi_list,
            'Reads': read_list,
            'AvgQual': quality_list
        })
        df_umi = df_umi[(df_umi['UMI'] != "") & (df_umi['Reads'] != "")]
        df_umi.to_csv(os.path.join(temp_dir, f'{sample_name}__Batch{chunk}.UMI_Reads.{suffix}.raw.txt'),
                      index=False, sep='\t')
        return df_umi

    df_umi123 = save_to_file(read_list123, umi_list123, quality_list123, target_list123, "PAMSide")
    df_umi456 = save_to_file(read_list456, umi_list456, quality_list456, target_list456, "ProtoSide")

    return [n, count123, count456, df_umi123, df_umi456]



def summarySOS(n, count_constant123, count_constant456, reverse, out_dir, sample_name, 
                  df_umi_constant123, df_umi_constant456, df_unique_umi_constant123, 
                  df_unique_umi_constant456, df_valid_umi):
    """
    Generate a summary of valid read and UMI statistics and save the results to files.

    Args:
        n (int): Total number of reads processed.
        count_constant123 (int): Count of valid reads for constant region 123 (proto-side).
        count_constant456 (int): Count of valid reads for constant region 456 (PAM-side).
        reverse (str): Indicator if reverse orientation is applied ("Yes" or "No").
        out_dir (str): Directory to save output files.
        sample_name (str): Name of the sample being processed.
        df_umi_constant123 (DataFrame): DataFrame of proto-side UMI reads.
        df_umi_constant456 (DataFrame): DataFrame of PAM-side UMI reads.
        df_unique_umi_constant123 (DataFrame): Proto-side UMI unique read counts.
        df_unique_umi_constant456 (DataFrame): PAM-side UMI unique read counts.
        df_valid_umi (DataFrame): DataFrame containing valid UMI counts.

    Returns:
        None: Saves summary statistics and UMI data to the specified output directory.
    """
    
    print("---- Generating Summary Statistics ----")
    
    # Calculate summary statistics
    total_valid_reads = count_constant123 + count_constant456
    total_valid_pct = (total_valid_reads / (n - 1)) * 100
    total_valid_umi = df_valid_umi['UMI_Count'].sum()
    umi_reads_ratio = total_valid_reads / total_valid_umi
    

    summary_headers=["total merged reads", "total_valid_reads", 
                    "Valid_reads_pct", "total_valid_umi", "ReadsUmiRatio"]
    
    summary_values = [
        n - 1,
        total_valid_reads, total_valid_pct, 
        total_valid_umi, umi_reads_ratio
    ]
    
    # Create summary DataFrame
    df_summary = pd.DataFrame({
        'variable': summary_headers,
        'values': summary_values
    }).assign(Sample_name=sample_name)
    
    # Pivot and save summary
    df_summary_pivot = df_summary.pivot(
        index='Sample_name', 
        columns='variable', 
        values='values'
    ).reset_index()
    
    summary_output_path = f"{out_dir}/{sample_name}.summary.txt"
    df_summary_pivot.columns = ['Sample_name', 'Read_UMI_ratio', 'Valid_read_percentage', 'Total_merged_reads', 'Total_valid_reads', 'Total_valid_UMIs']
    df_summary_pivot.to_csv(summary_output_path, index=False, sep='\t')
    
    # Sort UMI read data
    df_umi_constant123.sort_values(['Template_ID', 'UMI'], inplace=True)
    df_umi_constant456.sort_values(['Template_ID', 'UMI'], inplace=True)
    
    # Define file paths based on reverse flag
    if reverse == "Yes":
        side123_suffix = "PAMSide"
        side456_suffix = "ProtoSide"
    else:
        side123_suffix = "ProtoSide"
        side456_suffix = "PAMSide"
    
    # Save UMI read files
    df_umi_constant123.columns = ['Template_ID', 'UMI', 'Reads', 'Average_quality']
    df_umi_constant456.columns = ['Template_ID', 'UMI', 'Reads', 'Average_quality']
    df_umi_constant123.to_csv(
        f"{out_dir}/{sample_name}.UMI_Reads.{side123_suffix}.raw.txt", index=False, sep='\t'
    )
    df_umi_constant456.to_csv(
        f"{out_dir}/{sample_name}.UMI_Reads.{side456_suffix}.raw.txt", index=False, sep='\t'
    )

    df_unique_umi_constant123.columns = ['Template_ID', 'UMI', 'Sequence', 'Sequence_length', 'Average_quality', 'Read_count']
    df_unique_umi_constant456.columns = ['Template_ID', 'UMI', 'Sequence', 'Sequence_length', 'Average_quality', 'Read_count']
    df_unique_umi_constant123.to_csv(
        f"{out_dir}/{sample_name}.UMI_Reads.{side123_suffix}.UMI_unique.Readscount.txt", index=False, sep='\t'
    )
    df_unique_umi_constant456.to_csv(
        f"{out_dir}/{sample_name}.UMI_Reads.{side456_suffix}.UMI_unique.Readscount.txt", index=False, sep='\t'
    )
    
    print("---- Summary Generation Complete ----")


def preRunnerAnnotation(df_umi_count, job, out_dir, reverse, type, analysis_file, suffix_name, sample_name, temp_dir):
    """ 
    Prepare data (files and variable list) for parallel proto/PAM-side sequence annotation.

    Args:
        df_umi_count (DataFrame): DataFrame containing UMI counts per template.
        job (int): Number of parallel jobs to distribute the workload.
        out_dir (str): Directory to save output annotation results.
        reverse (str): Reverse sequence indicator for the sample.
        type (str): Editing type (e.g., Cas9, ABE).
        analysis_file (str): Path to the analysis file used for annotation.
        suffix_name (str): Suffix to append to the output file name.
        sample_name (str): Sample name to use for file labeling.
        temp_dir (str): Directory for storing temporary files.

    Returns:
        list: A list of lists, where each sublist contains parameters needed for 
              parallel annotation tasks (chunk ID, file paths, sample metadata).
    """

    # Assing chunk identifier to each template for parallele processing, 
    if len(df_umi_count)<job:
        job=len(df_umi_count)
    chunk_size=len(df_umi_count)/(job-1)
    df_umi_count=df_umi_count.assign(inx=df_umi_count.index.map(int).copy())
    df_umi_count=df_umi_count.assign(chunk=round(df_umi_count['inx'].map(int)/int(chunk_size),0))
    df_umi_count['chunk']=df_umi_count['chunk'].map(int)

    # Save dataframe with chunk identifier to temporary folder
    sample_file=os.path.join(temp_dir, f'{sample_name}.{suffix_name}')
    df_umi_count.to_csv(sample_file)

    # Build a list of inputs for parallele pocessing at annotation step
    runner_list_annotation=[]
    for chunk in df_umi_count['chunk'].unique():
        runner_list_annotation.append([chunk, sample_file, out_dir, reverse, analysis_file, type, sample_name])

    return runner_list_annotation



def calculateUMIPercentage(df):
    """
    Calculate the percentage of UMI counts relative to the total UMI counts per TemplateID.

    Args:
        df (DataFrame): Input DataFrame containing 'Template_ID' and 'UMI_Count' columns.

    Returns:
        DataFrame: Updated DataFrame with an additional 'UMI_Pct' column representing 
                   the percentage of UMI counts for each row within the same TemplateID.
    """
    df_umi_sum = df.groupby(['Template_ID'])['UMI_Count'].sum().reset_index()
    df_umi_sum.columns=['Template_ID', 'SumCount']
    df = pd.merge(df,df_umi_sum, on='Template_ID', how="left" )
    df = df.assign(UMI_Pct=df.UMI_Count/df.SumCount*100)
    del df['SumCount']
    return df



def MPIAnnotateEditStatusC123(runner):
    """
    Align and annotate proto-side sequences using parallel processing.

    Args:
        runner (list): A list containing parameters for processing, including:
            - chunk (int): Current chunk identifier.
            - sample_file (str): Path to the sample UMI file.
            - out_dir (str): Output directory for results.
            - reverse (str): Reverse sequence indicator.
            - analysis_file (str): Path to the analysis file.
            - type (str): Editing type (ABE, Cpf1, etc.).
            - sample_name (str): Name of the sample being processed.

    Returns:
        DataFrame: Annotated DataFrame with editing status for each sequence.
    """
    # Unpack runner parameters
    chunk, sample_file, out_dir, reverse, analysis_file, type, sample_name = runner

    # Load analysis and sample data
    dfAnalysis = pd.read_csv(analysis_file, sep='\t')
    df_umi_count = pd.read_csv(sample_file)
    df_umi_count = df_umi_count[df_umi_count['chunk'] == chunk]

    print(f"---- Annotation Edit on c1 side: {chunk} ----")
    
    # Initialize EditingStatus column
    df_umi_count['EditingStatus'] = "NA"
    
    # Extract unique TemplateID and Sequence combinations
    df_umi_count_unique = df_umi_count[['Template_ID', 'Sequence']].drop_duplicates()

    # Annotate sequences
    for i, (templateID, seq) in enumerate(zip(df_umi_count_unique['Template_ID'], df_umi_count_unique['Sequence']), 1):
        
        if i % 1000 == 0:
            print(f"----- Processing {chunk} : {i} sequences -----")
        
        # Retrieve template and expected cut/user information
        template_row = dfAnalysis[dfAnalysis['Template_ID'] == templateID]
        template = template_row['Template'].values[0]
        expcut = template_row['Expected_nick'].values[0]
        
        # Perform editing annotation based on editing type
        if type == "ABE":
            seq = seq[:-1] if seq.endswith("A") else seq
        editestatus, inseq, refseq = checkEditedStatusC123(seq, template, expcut=expcut, type=type)
        
        # Update editing status in the DataFrame
        df_umi_count.loc[
            (df_umi_count['Sequence'] == seq) & 
            (df_umi_count['Template_ID'] == templateID), 
            'EditingStatus'
        ] = editestatus
        
    return df_umi_count



def MPIAnnotateEditStatusC456(runner): 
    """
    Align and annotate PAM-side sequences using parallel processing.

    Args:
        runner (list): A list containing parameters for processing, including:
            - chunk (int): Current chunk identifier.
            - sample_file (str): Path to the sample UMI file.
            - out_dir (str): Output directory for results.
            - reverse (str): reverse sequence indicator.
            - analysis_file (str): Path to the analysis file.
            - type (str): Editing type (Cas9, Cpf1, CBE, ABE).
            - SampleName (str): Name of the sample being processed.

    Returns:
        DataFrame: Annotated DataFrame with editing status for each sequence.
    """
    
    chunk, sample_file, out_dir, reverse, analysis_file, type, SampleName = runner
    
    # Load analysis and sample data
    df_analysis = pd.read_csv(analysis_file, sep='\t')
    df_umi_count = pd.read_csv(sample_file)
    df_umi_count = df_umi_count[df_umi_count.chunk == chunk]
    
    print(f"---- Annotation Edit on r1 side: {chunk} ----")
    
    # Initialize EditingStatus column
    df_umi_count['EditingStatus'] = "NA"
    
    # Extract unique TemplateID and Sequences
    dfr1umicountinfo = df_umi_count[['Template_ID', 'Sequence']].drop_duplicates()
    
    # Annotate sequences
    for i, (templateID, seq) in enumerate(zip(dfr1umicountinfo['Template_ID'], dfr1umicountinfo['Sequence']), 1):
        
        if i % 1000 == 0:
            print(f"----- Processing {chunk} : {i} sequences -----")
        
        # Retrieve template and expected cut/user information
        template_row = df_analysis[df_analysis.Template_ID == templateID]
        template = template_row['Template'].values[0]
        expcut = template_row['Expected_nick'].values[0]
        expuser = template_row['Expected_USER'].values[0]
        
        # Perform edited UMI annotation
        if type == "CBE":
            expcut = expuser
        elif type == "ABE":
            seq = seq[1:] if seq.startswith("T") else seq
        editestatus, inseq, refseq = checkEditedStatusC456(seq, template, expcut, type)

        # Update editing status in the DataFrame
        df_umi_count.loc[
            (df_umi_count['Sequence'] == seq) & 
            (df_umi_count['Template_ID'] == templateID), 
            'EditingStatus'
        ] = editestatus
        
    return df_umi_count