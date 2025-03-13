import os
import copy
import math
import pysam
import pybedtools
import numpy as np
import pandas as pd
import ray
from ray.util.multiprocessing import Pool
from IPython.display import display

from modules.Alignment import *



def parseCoordinate(x):
    """
    Parse coordiante string into 3 variables

    Parameters:
        x (string): such as "chr1:123-678"

    Returns:
        chr (string): chromosome or contig name 
        start (integer): start position
        end (integer): end position
    """
    chr = ''
    start = end = -1

    # assign value to chr, start, end
    if ':' in x:
        array = x.replace(':', '-').split('-')
        if len(array) == 3:
            [chr, start, end] = array
            start = int(start)
            end = int(end)
        elif len(array) == 2:
            [chr, start] = array
            start = int(start)
    else:
        chr = x

    return chr, start, end



def bitsContainsRedFlag(bits):
    """
    Test if input binary code encode red flags

    Parameters:
        bits (string): a string of 12 digits
    Returns 
        boolean: True or False
    """
    red_flag_index = [ 3,   #BAM_FUNMAP(4),
                       9,   #BAM_FSECONDARY(256)
                      10,   #BAM_FQCFAIL(512)
                      11]   #BAM_FDUP(1024)
    for i in red_flag_index:
        flag = bits[0 - i]
        if flag == '1':
            return True
    return False



def isRedFlag(x):
    """
    Test if flag associated with a read contains unwanted scenarions

    Parameters:
        x (integer): an integer from pysam read.flag

    Returns 
        a boolean, True or False
    """

    def _integerTo12BitBinary(x):
        binary_str = bin(x).replace('0b', '')
        binary_str = '0' * 12 + binary_str
        binary_str = binary_str[-12:]
        return binary_str

    return bitsContainsRedFlag(_integerTo12BitBinary(x))



def coordinatesOnEdgesInOneRegion(sam_object, region, read_count_cutoff=3, mapping_quality=1, strand='both'):
    """
    Search for coordinates in a single region of bam file that match the beginning of multiple reads using 
    a single processor

    Parameters:
        sam_object (AlignmentFile object): an object returned from class pysam.AlignmentFile
        region (string): such as 'chr1:123-678'
        read_count_cutoff (integer): minimal number of reads to be considered
        mapping_quality (integer): minimal mapping quality for a mapped read to considered
        strand (string): one of the 3 options: 'both', 'forward' and 'reverse'
        outfile (string): path to a local file to write output dataframe to

    Returns:
        dataframe: dataframe with index for coordinates and a column 'count' for supporting reads
    """
    site_counts = {}

    chr, start, end = parseCoordinate(region)
    end = end if end >=0 else start + 1

    # get sam fetch object
    if chr and start >=0:
        fetch_object = sam_object.fetch(chr, start, end)
    elif chr:
        fetch_object = sam_object.fetch(chr)

    # iterate through fetch object and summaize position align with begining of reads
    for read in fetch_object:
        if read.mapping_quality < mapping_quality or isRedFlag(read.flag):
            continue    # skip reads red-flagged or of low mapping quality

        try:
            all_positions = read.get_aligned_pairs(matches_only=False)   ####

            # skip if no alignment
            if len(all_positions) < 1:
                continue
            # skip if strand requirement does'nt meet
            is_forward = read.is_forward
            if (is_forward and strand=='reverse') or (not is_forward and strand=='forward'):
                continue

            # extract alignment start position in reference genome
            start_position = ''
            if is_forward:
                start_position = all_positions[0][1]
            else:
                start_position = all_positions[-1][1]
            if not start_position:
                continue

            if type(start_position) is int:
                site = f"{chr}:{start_position+1}:{'+' if is_forward else '-'}"
                if site in site_counts:
                    site_counts[site] += 1
                else:
                    site_counts[site] = 1                
        except:
            print(read)
            print(all_positions)
            pass

    # remove position with low supporting reads
    df = pd.DataFrame.from_dict(site_counts, orient='index', columns=['count'])
    df2 = df[df['count'] >= read_count_cutoff]
    df2_with_flanking_sites = df2.copy(deep=True)

    """Get a second dataframe with additional coordinates (regardless of counting 
    of reads aligned to it) within 3bp-flanking region of those included in df2.
    this dataframe is to be used down the road in the scoring step to calculate 
    supporting ratio for flanking positions"""
    if df2.shape[0] > 0:
        df2_expanded = df2.copy(deep=True).reset_index()
        df2_expanded[['chr', 'start', 'strand']] = df2_expanded['index'].str.split(':', expand=True, regex=False)
        df2_expanded['start'] = df2_expanded['start'].astype(int)

        dfs_for_flanking_sites = [df2_expanded]
        for i in [-3, -2, -1, 1, 2, 3]:
            dfx = df2_expanded.copy(deep=True)
            dfx['start'] = dfx['start'] + i
            dfs_for_flanking_sites.append(dfx)
        df_with_flanking_sites = pd.concat(dfs_for_flanking_sites)
        df_with_flanking_sites['coordinates'] = df_with_flanking_sites['chr'] + ':' + df_with_flanking_sites['start'].astype(str) + ':' + df_with_flanking_sites['strand']
        df_with_flanking_sites = df_with_flanking_sites.drop_duplicates('coordinates')
        df2_with_flanking_sites = df[df.index.isin(df_with_flanking_sites['coordinates'])]

    return df2, df2_with_flanking_sites



def coordinatesOnEdgesInMultipleRegions(argument_dict):
    """
    Search for coordinates in a given bam file for postions matching the beginning of multiple reads using 
    a single processor

    Parameters:
        argument_dict (dict): a dict with the following key-value pairs:
            bam_path: a path to an bam file
            chromosomes: a list of chromosomes or contigs
            read_count_cutoff: minimal number of reads to be considered
            mapping_quality: minimal mapping quality for a mapped read to considered
            strand: one of the 3 options: 'both', 'forward' and 'reverse'

    Returns:
        dataframe: a dataframe with columns "coordinate" (example value: 'chr19:51167431:+') and 'count' (example value '20')
    """
    # parse input arguments
    samfile = argument_dict['bam_path']
    chromosomes = argument_dict['chromosomes']
    read_count_cutoff = argument_dict['read_count_cutoff'] if 'read_count_cutoff' in argument_dict else 5
    mapping_quality = argument_dict['mapping_quality'] if 'mapping_quality' in argument_dict else 1
    strand = argument_dict['strand'] if 'strand' in argument_dict else 'both'

    # read input bam file
    if isinstance(samfile, str) and os.path.exists(samfile):
        samfile = pysam.AlignmentFile(samfile, "rb")

    # extract infomation from pysam object
    dfs = []
    dfs_with_flanking_sites = []
    for chr in list(set(chromosomes)):
        df, df_with_flanking_sites = coordinatesOnEdgesInOneRegion(samfile, chr, read_count_cutoff, mapping_quality, strand)
        dfs.append(df)
        dfs_with_flanking_sites.append(df_with_flanking_sites)
    
    # combine data from multiple chromosomes
    all = pd.concat(dfs).reset_index().rename(columns={'index':'coordinates'})
    all_with_flanking_sites = pd.concat(dfs_with_flanking_sites).reset_index().rename(columns={'index':'coordinates'})
    return [all, all_with_flanking_sites]



def getCoordinatesOnEdges(sam_path, read_count_cutoff=3, processor_count=5, mapping_quality=1, strand='both'):
    """
    Search for coordinates in a given bam file for postions matching the beginning of multiple reads using 
    multiple processor

    Parameters:
        sam_path (string) path to an bam file
        read_count_cutoff (integer): minimal number of supporing reads to be considered
        processor_count (integer): number of processors to be used for parellel computing
        mapping_quality (integer): minimal mapping quality for a mapped read to considered
        strand (string): one of the 3 options: 'both', 'forward' and 'reverse'

    Returns
        dataframe: a dataframe with columns "coordinate" (example value: 'chr19:51167431:+') and 'count' (example value '20')
    """
    # devide chr/contigs in bam into several groups that are to be processed in parellel
    reference_groups = groupReferences(sam_path, processor_count)

    # start ray instance 
    ray.init(ignore_reinit_error=True)
    print(f'Screen of sites aligned with begining of multiple reads using {processor_count} processors')

    x = []
    y = []
    pool = Pool(processor_count)
    for result in pool.map(
        coordinatesOnEdgesInMultipleRegions, 
        [{'bam_path': sam_path, 
          'chromosomes': x, 
          'read_count_cutoff': read_count_cutoff, 
          'mapping_quality': mapping_quality, 
          'strand':strand} for x in reference_groups]):
        x.append(result[0])
        y.append(result[1])
    ray.shutdown()

    return pd.concat(x), pd.concat(y)



def assignSupportingReadCountToPositions(df, df_breaks_flanking, coordinate_col):
    """
    Assign supporting read count to positions given a supporting read count database. 
    The positions centered around potential cut sites and is listed under column 'neighbor(f/r)' in input dataframe

    Parameters:
        df (dataframe): a dataframe with column 'neighbor(f)' or 'neighbor(r)'
        df_breaks_flanking (dataframe): a dataframe contain supporting read counts for a list of sites
        coordinate_col (string): either 'neighbor(f)' or 'neighbor(r)'

    Returns
        a copy of input dataframe with new column 'neighbor_count(f)' or 'neighbor_count(r)'
    """

    breaks = df_breaks_flanking.copy(deep=True)
    breaks = breaks.drop_duplicates()
    new_column_name = 'neighbor_count'
    if '(r)' in coordinate_col:
        breaks = breaks[breaks['coordinates'].str.contains('-')]
        new_column_name = new_column_name + '(r)'
    elif '(f)' in coordinate_col:
        breaks = breaks[breaks['coordinates'].str.contains('+', regex=False)]
        new_column_name = new_column_name + '(f)'

    # manipulate data for merged
    if breaks.shape[0]:
        breaks = breaks.rename(columns={'count': new_column_name})
        breaks['coordinates'] = breaks['coordinates']\
            .str.replace('-', '', regex=False)\
            .str.replace('+', '', regex=False)\
            .str.replace(':$', '', regex=True)
    else:
        return df
    
    # append supporting read cout to input dataframe
    new_df = df.merge(breaks, left_on=coordinate_col, right_on='coordinates', how='left')
    new_df = new_df.fillna(0)
    new_df[new_column_name] = new_df[new_column_name].astype(int)
    new_df.drop('coordinates', axis=1, inplace=True)

    return new_df



def groupReferences(samfile, processor_count):
    """
    Separate chromosomes/contigs into a number of groups so that the total length of each group are about the same as possible

    Parameters:
        samfile (string) path to an bam file
        processor_count (integer) number of processors to be used for a relevant job

    Returns 
        list: a list of chromosome/contig groups. Each group is a list.
    """
    # read input bam file
    if isinstance(samfile, str) and os.path.exists(samfile):
        samfile = pysam.AlignmentFile(samfile, "rb")
    else:
        raise Exception(f"Invalid alignment file or file does not exist {samfile}")

    # extract reference names and lengths
    references = samfile.references
    reference_lengths = {}
    total_reference_length = 0
    for ref in references:
        length = samfile.get_reference_length(ref)
        reference_lengths[ref] = length
        total_reference_length = total_reference_length + length
    
    reference_lengths = dict(sorted(reference_lengths.items(), key=lambda x:x[1], reverse=False))

    # group references
    reference_groups = []
    fragment_size = math.ceil(total_reference_length / processor_count)

    a_reference_group = []
    a_reference_group_total_length = 0
    for ref, ref_length in reference_lengths.items():
        a_reference_group_total_length += ref_length
        a_reference_group.append(ref)
        if a_reference_group_total_length > fragment_size:
            reference_groups.append(copy.deepcopy(a_reference_group))
            a_reference_group = []
            a_reference_group_total_length = 0

    if len(a_reference_group) > 0:
        reference_groups.append(copy.deepcopy(a_reference_group))

    return reference_groups



def coverageAtOnePosition(coordinate, sam_object, checkARead):
    """
    Extract read coverage at a position

    Parameters:
        coordinate (string): such as "chr1:12345"
        sam_object (AlignmentFile object): an object returned from class pysam.AlignmentFile

    Returns 
        read_count (integer): read count (coverage or depth) for the input coordinate
    """

    array = coordinate.split(':')
    chr = array[0]
    positition = int(array[1])
    # get coverage from BAM
    count_tuple = sam_object.count_coverage(chr, positition-1, positition, 
                                            quality_threshold = 0,
                                            read_callback=checkARead)
    # do the sum
    read_count = 0
    for array in count_tuple:
        for count in array:
            read_count += count

    return read_count



def getSequenceCoverages(argument_dict):
    """
    Add coverage column to input dataframe

    Parameters:
        argument_dict (dict): a dictionary with key-value pairs:
        'df': a dataframe with column 'coordinates'
        'bampath': a path leading to a bam file
        
    Returns:
        dataframe: a copy of input dataframe with additional column 'coverage'
    """
    # Note: passing a pysam object 1) does not work with ray multiprocessing, 2) might cause 
    # discrepancy in results (https://groups.google.com/g/pysam-user-group/c/bBdqn7DkVtE)

    df = argument_dict['df'] if 'df' in argument_dict else pd.DataFrame(columns=['coordinates', 'coverage'])
    bampath = argument_dict['bampath'] if 'bampath' in argument_dict else ''
    checkARead = argument_dict['func']


    if df.shape[0] > 0 and bampath and os.path.exists(bampath):
        samfile = pysam.AlignmentFile(bampath, "rb")
        df['coverage'] = df['coordinates'].apply(coverageAtOnePosition, args=[samfile, checkARead])
        samfile.close()

    return df



def getSequenceCoveragesMPI(bam_path, df, processor_count, checkARead):
    """
    Retrieve coverage for a list of coordinates using more than one processor

    Parameters:
        bam_path (string): path to local bam file
        df (dataframe): it must contain a column named 'coordinates' with values like "chr1:23456"
        processor_count (integer): number of processor to leverage for parallel

    Returns
        dataframe: a new dataframe with coordinates and their corresponding read coverage
    """
    # figure out batch size and batch number
    cut_site_count = df.shape[0]
    batch_size = math.ceil(cut_site_count / processor_count)
    batch_number = math.ceil(cut_site_count / batch_size)
    # split input df into subsets
    coordinate_subgroups = []
    for batch in range(0, batch_number):
        index_start = batch * batch_size
        index_end  = index_start + batch_size
        index_end = index_end if index_end <= cut_site_count else cut_site_count
        coordinate_subgroups.append(df[index_start:index_end])


    # retrieve coverage in multi-threading manner
    ray.init(ignore_reinit_error=True)
    x = []
    pool = Pool()
    for result in pool.map(getSequenceCoverages, [{'df':x, 'bampath':bam_path, 'func':checkARead} for x in coordinate_subgroups]):
        x.append(result)
    ray.shutdown()

    # combined output from different threadings
    df_new = pd.concat(x)
    df_new = df_new.sort_values(['coordinates'])
    return df_new



def addCoverageColumns(bam_path, df, processor_count, checkARead):
    """
    Retrieve coverage for forward/reverse coordinates using more than one processor

    Parameters:
        bam_path (string): path to local bam file
        df (dataframe): it must contain 2 columns named 'coordinate(f)' and 'coordinate(r)' with values like "chr1:23456"
        processor_count (integer): number of processor to leverage for parallel

    Returns
        dataframe: a new dataframe with coordinats and new column for "sub score" 
    """
    # get a set unique cooridantes
    coordiante_f = df['coordinate(f)'].unique()
    coordiante_r = df['coordinate(r)'].unique()
    coordiante_for_coverage = np.concatenate([coordiante_f, coordiante_r])
    df_coordinates_for_coverage = pd.DataFrame(coordiante_for_coverage, columns=['coordinates'])
    df_coordinates_for_coverage = df_coordinates_for_coverage.drop_duplicates()
    print(f'{df_coordinates_for_coverage.shape[0]} unique positions to retrive covrage for')

    # retrive coverage
    df_coverage = getSequenceCoveragesMPI(bam_path, df_coordinates_for_coverage, processor_count, checkARead)
    print(f'received coverage for {df_coverage.shape[0]} positions')
    display(df_coverage.head())

    # add coverage to input dataframe
    df_coverage_f = df_coverage.copy(deep=True)
    df_coverage_f = df_coverage_f.drop_duplicates()
    df_coverage_f = df_coverage_f.rename(columns={'coordinates':'coordinate(f)', 'coverage':'coverage(f)'})

    df_coverage_r = df_coverage.copy(deep=True)
    df_coverage_r = df_coverage_r.drop_duplicates()
    df_coverage_r = df_coverage_r.rename(columns={'coordinates':'coordinate(r)', 'coverage':'coverage(r)'})

    df = df.merge(df_coverage_r, on='coordinate(r)', how='left')
    df = df.merge(df_coverage_f, on='coordinate(f)', how='left')
    return df



def neighborIntegers(x, initial, moves):
    """Given a number, derive a list of numbers surround it"""
    distances = list(range(1, moves+1)) if moves > 0 else list(range(moves, 0))
    return [x + initial + distance for distance in distances]



def neighborPositionsRequiredForScoring(position, initial=-4, moves=5):
    """Given a coordinate, derive a list of coordinates immediantely flanking it"""
    chr, start, end = parseCoordinate(position)
    return [f'{chr}:{x}' for x in neighborIntegers(start, initial, moves)]



def CalculateSubscores(samfile, df_breaks_flanking, df_input, func, orientation, initial=-4, moves=5):
    """
    Calculate single side scores for each cut site

    Parameters:
        samfile (AlignmentFile object): an object returned from class pysam.AlignmentFile
        df_breaks_flanking (dataframe): a dataframe contain supporting read counts for a list of sites
        df_input (dataframe): with at least 4 columns: coordinate{orientation}, coordinate{orientation}', count{orientation}, coverage{orientation}
        orientation (string): either "f" (forward) or "r" reverse
        initial (integer): initial position relevant to the coordinate in question
        moves (integer): number of sequential basepairs from initial position to be included in score calculation

    Returns
        dataframe : a new dataframe with coordinate columns and new column for "sub-score"
    """
    # sanity check
    if orientation not in ['f', 'r']:
        print(f'Invalde orientation input: {orientation}')
        return df
    opposite_orientation = 'r' if orientation =='f' else 'f'

    # fetch neighbouring positions
    df = df_input.copy(deep=True)
    df[f'neighbor({opposite_orientation})'] = df[f'coordinate({opposite_orientation})'].apply(neighborPositionsRequiredForScoring, args=[initial, moves])
    df = df.explode(f'neighbor({opposite_orientation})')
    df = df[[f'coordinate({orientation})', f'coordinate({opposite_orientation})', f'count({orientation})', f'coverage({orientation})', f'neighbor({opposite_orientation})']]
    # get coverage for neighbouring positions
    df[f'neighbor_coverage({opposite_orientation})'] = df[f'neighbor({opposite_orientation})'].apply(coverageAtOnePosition, args=[samfile, func])
    # get supporting count for neighbouring positions
    df = assignSupportingReadCountToPositions(df, df_breaks_flanking, f'neighbor({opposite_orientation})')

    # calculate subscore
    df['subscore-ratio'] = ((df[f'count({orientation})'] -1) / df[f'coverage({orientation})']) * ((df[f'neighbor_count({opposite_orientation})'] - 1)/df[f'neighbor_coverage({opposite_orientation})'])
    df['subscore'] = df['subscore-ratio'] * (df[f'count({orientation})'] + df[f'neighbor_count({opposite_orientation})'] - 2)
    # do sum to get final score
    df = df.groupby([f'coordinate({orientation})', f'coordinate({opposite_orientation})']).sum(numeric_only=True)
    df = df.reset_index()
    df = df[[f'coordinate({orientation})', f'coordinate({opposite_orientation})', 'subscore']]

    # grab useful columns
    df.columns = [f'coordinate({orientation})', f'coordinate({opposite_orientation})', f'subscore({orientation})']
    return df



def calculateScores(samfile, df_breaks_flanking, df, func):
    """
    Calculate cleavage score for each pair of coordinates (forward and reverse)

    Parameters:
        samfile (AlignmentFile object): an object returned from class pysam.AlignmentFile
        df_breaks_flanking (dataframe): a dataframe contain supporting read counts for a list of sites
        df (dataframe): a pandas dataframe with at least 2 columns 'coordinate(f)' and 'coordinate(r)'

    Returns
        dataframe : a copy of input dataframe with new column "score"
    """
    subscore_forward = CalculateSubscores(samfile, df_breaks_flanking, df, func, 'f', -3, 5)
    subscore_reverse = CalculateSubscores(samfile, df_breaks_flanking, df, func, 'r', 3, -5)

    scores = df.merge(subscore_forward, on=['coordinate(f)', 'coordinate(r)'], how='left')
    scores = scores.merge(subscore_reverse, on=['coordinate(f)', 'coordinate(r)'], how='left')
    scores['score'] = scores[['subscore(f)', 'subscore(r)']].sum(axis=1)
    scores = scores.drop(['subscore(f)', 'subscore(r)'], axis=1)
    scores = scores.sort_values('score', ascending=False)
    return scores



def getCoordinatePairsByRelativeDistance(df1, df2, allowed_distances=[1,2,3]):
    """
    Find paris of coordiantes associated with each other given defined distance between them

    Parameters:
        df1 (dataframe): dataframe for the first set of coordinates. Must contain at least 2 columns 'coordinates' and 'count'
        df2 (dataframe): dataframe for the first set of coordinates. Must contain at least 2 columns 'coordinates' and 'count'
        allowed_distances (list): a list of allowed distance between 2 coordinates

    Returns:
        dataframe: dataframe with each row for a pair of coordinates
    """
    allowed_distances.sort()

    # prepare bed for forwardly aligned reads
    dff = df1.copy(deep=True)
    dff[['chr', 'start', 'strand']] = dff['coordinates'].str.split(":", expand=True)
    dff['start'] = dff['start'].astype(int)
    dff['end'] = dff['start']
    dffbed = dff[['chr', 'start', 'end', 'coordinates']]
    # display(dffbed.head())

    # prepare bed for reversely aligned reads
    dfr = df2.copy(deep=True)
    dfr[['chr', 'start', 'strand']] = dfr['coordinates'].str.split(":", expand=True)
    dfr['start'] = dfr['start'].astype(int)
    dfr['end'] = dfr['start'] + 1 - allowed_distances[0]
    dfr['start'] = dfr['start'] + 1 - allowed_distances[-1]   # consider overlapping distance
    dfrbed = dfr[['chr', 'start', 'end', 'coordinates']]

    # intersect two sets of bed regions
    bed_object_f = pybedtools.BedTool.from_dataframe(dffbed)
    bed_object_r = pybedtools.BedTool.from_dataframe(dfrbed)
    overlaps = bed_object_f.intersect(bed_object_r, wo=True).to_dataframe()

    # extract selected columns    
    overlaps = overlaps.rename(columns={'name':'coordinate(f)', 'thickEnd':'coordinate(r)'})
    overlaps = overlaps[['coordinate(f)', 'coordinate(r)']]

    # calculate overlapping length and refine row selection
    overlaps[['chr', 'position1', 'strand']] = overlaps['coordinate(f)'].str.split(":", expand=True)
    overlaps['position1'] = overlaps['position1'].astype(int)
    overlaps[['chr', 'position2', 'strand']] = overlaps['coordinate(r)'].str.split(":", expand=True)
    overlaps['position2'] = overlaps['position2'].astype(int)
    overlaps['gamma'] = overlaps['position2'] - overlaps['position1'] + 1
    overlaps.drop(['chr', 'position1', 'strand', 'position2'], axis=1, inplace=True)
    overlaps = overlaps[overlaps['gamma'].isin(allowed_distances)]

    # combind intersetion with read count data
    df_overlap = df1.merge(overlaps, left_on='coordinates', right_on='coordinate(f)', how='inner')
    df_overlap = df_overlap.merge(df2, left_on='coordinate(r)', right_on='coordinates', how='inner')
    df_overlap['count(f)'] = df_overlap['count_x']
    df_overlap['count(r)'] = df_overlap['count_y']
    df_overlap['coordinate(f)'] = df_overlap['coordinate(f)'].str.replace(r':[-\+]$', '', regex=True)
    df_overlap['coordinate(r)'] = df_overlap['coordinate(r)'].str.replace(r':[-\+]$', '', regex=True)

    # extract useful features
    df_overlap = df_overlap[['coordinate(f)', 'count(f)', 'coordinate(r)', 'count(r)', 'gamma']]
    return df_overlap



def wgs_analysis(bam_path, cutoff_supporting_reads_forward, cutoff_supporting_reads_reverse, cutoff_mapping_quality, gamma, cutoff_ratio, cutoff_coverage, cutoff_score, output_file, gRNA, pam_location, enzyme, reference_genome, checkARead, processor_number):
    """
    Do a full-blown wgs analysis given a input bam file

    Parameters:
        bam_path (string): path to local bam file
        cutoff_supporting_reads_forward (int): minimal number of supporing reads (forwardly aligned) to be considered
        cutoff_supporting_reads_reverse (int): minimal number of supporing reads (reversely aligned) to be considered
        cutoff_mapping_quality (integer): minimal mapping quality for a mapped read to considered
        gamma (int): number of base pair allowed between 2 coordiantes
        cutoff_ratio (float): minimal fraction of reads supporting edit at a position in reference
        cutoff_coverage (int): minimal count of reads that mapped to a position in reference
        cutoff_score (float): minimal score for a cut site to be considered
        gRNA (string): guide RNA sequence
        reference_genome (string): path to a reference genome
        processor_number (int): number of processor to be used for parallel computing

    Returns:
        dataframe: a dataframe contains edited cuts meeting requirements
    """


    print(f'Searching cut sites in input bam file {bam_path}')

    # read bam file 
    samfile = pysam.AlignmentFile(bam_path, "rb")
    dict_chromosomes_lengths = getChromosomesAndLengths(samfile)

    # initial screen of sites aligned with begining of multiple reads forwardly aligned to genome
    df_breaks, df_breaks_flanking = getCoordinatesOnEdges(bam_path, 
                                min([cutoff_supporting_reads_forward, cutoff_supporting_reads_reverse]),
                                processor_number, cutoff_mapping_quality)
    print(f'{df_breaks.shape[0]} positions align to {cutoff_supporting_reads_forward} or more forwardly/reversely mapped reads')
    display(df_breaks.head())
    df_breaks_forward = df_breaks[df_breaks['coordinates'].str.endswith('+')]
    df_breaks_forward = df_breaks_forward[df_breaks_forward['count'] >= cutoff_supporting_reads_forward]
    df_breaks_reverse = df_breaks[df_breaks['coordinates'].str.endswith('-')]
    df_breaks_reverse = df_breaks_reverse[df_breaks_reverse['count'] >= cutoff_supporting_reads_reverse]
    
    # search site-pairs that match specified overlapping length
    df_site_pairs = getCoordinatePairsByRelativeDistance(df_breaks_forward, df_breaks_reverse, gamma)
    print(f'{df_site_pairs.shape[0]} coordinate pairs (forward & reverse) meet overlapping requirements')

    # retrive coverages for all sites
    df_site_pairs =  addCoverageColumns(bam_path, df_site_pairs, processor_number, checkARead)
    # get ratios
    df_site_pairs = df_site_pairs.drop_duplicates()
    df_site_pairs['ratio(f)'] = df_site_pairs['count(f)'].astype(int) / df_site_pairs['coverage(f)'].astype(int)
    df_site_pairs['ratio(r)'] = df_site_pairs['count(r)'].astype(int) / df_site_pairs['coverage(r)'].astype(int)
    # filter sites by ratio
    df_site_pairs = df_site_pairs[(df_site_pairs['ratio(f)']>=cutoff_ratio) & (df_site_pairs['ratio(r)']>=cutoff_ratio)]

    # scoring
    df_site_pairs = calculateScores(samfile, df_breaks_flanking, df_site_pairs, checkARead)

    # add guide-reference alignemnt
    df_site_pairs = addSequenceAlignmentFeatures(df_site_pairs, reference_genome, gRNA, 40, pam_location, enzyme)

    # reformat coordinates
    df_site_pairs[['chr', 'coordinate(f)']] =df_site_pairs['coordinate(f)'].str.split(':', expand=True)
    df_site_pairs[['chr', 'coordinate(r)']] =df_site_pairs['coordinate(r)'].str.split(':', expand=True)
    leading_features =  ['chr', 'coordinate(f)', 'coordinate(r)', 
                        'count(f)', 'count(r)', 'coverage(f)', 'coverage(r)', 'ratio(f)', 'ratio(r)', 
                        'score']
    df_site_pairs = df_site_pairs[leading_features + [x for x in df_site_pairs.columns if x not in leading_features]]

    # more filtering 
    df_site_pairs = df_site_pairs[df_site_pairs['score']>=cutoff_score]
    df_site_pairs = df_site_pairs[df_site_pairs['coverage(f)']>=cutoff_coverage]
    df_site_pairs = df_site_pairs[df_site_pairs['coverage(r)']>=cutoff_coverage]

    # add chromosome length
    df_site_pairs['chr_length'] = df_site_pairs['chr'].apply(getChromosomeLength, args=[dict_chromosomes_lengths])
    display(df_site_pairs)


    # export
    print(f'A total of {df_site_pairs.shape[0]} coordinate pairs (forward & reverse) in all categories')
    if output_file:
        df_site_pairs.to_csv(output_file, index=False)
        display(df_site_pairs.head())
    else:
        display(df_site_pairs.sort_values(['chr', 'coordinate(f)']))

    samfile.close()

    return df_site_pairs


def getChromosomesAndLengths(sam_object):
    """
    """

    chromosomes = sam_object.references
    lengths = sam_object.lengths

    dict = {}
    if len(chromosomes)>0 and len(chromosomes) == len(lengths):
        for i in range(len(lengths)):
            dict[chromosomes[i]] = lengths[i]
    return dict
    
    

def getChromosomeLength(chr, chromosomes_lengths):
    """
    """
    
    length = 0
    if chr in chromosomes_lengths:
        length = chromosomes_lengths[chr]
    return int(length)
