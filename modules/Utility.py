import os
import pandas as pd
from modules.Bed import addCoordinateColumn, getOverlappingCoordinates
from modules.Visualization import manhattanPlot, formatedExcelSheet, splitFileNameBySuffix



def getOutputFilePaths(file_path):
    """
    Create output file path based on user input arguments

    Parameters:
        file_path (string): file path
        background_gamma (list<number>): a list gamma used to define background sites
        enzyme (string): name of DNA modifying enzyme

    Returns:
        string: file path to save all background sites
        string: file path to save all foreground sites that does't overlap with background sites
        string: file path to save all foreground sites that overlap with background (potential artifacts)
        string: file path to save all sites found in control bam file
    """

    file_name, suffix = splitFileNameBySuffix(file_path)

    if file_name and suffix:
        return f'{file_name}.control.{suffix}'
    elif suffix:
        f'{suffix}.control.csv'
    else:
        return 'control.csv'



def sanityCheckInputFile(file_path, file_type):
    """
    Check existance of input file and its companion file if any

    Parameters:
        file_path (string): file path to an input file
        file_type (sring): type of input file including 'bam' or 'fasta'

    Returns
        list<string>: a list of strings encoding error message
    """
    errors = []
    companion_file = ''
    
    # derive companion file to be checked
    if file_path and file_type.lower() == 'bam':
        companion_file = f'{file_path}.bai'
    if file_path and file_type.lower() == 'fasta':
        companion_file = f'{file_path}.fai'

    # check existance of file and its companion file
    if file_path and not os.path.exists(file_path):
        errors.append(f'{file_path} does not exists')
    if companion_file and not os.path.exists(companion_file):
        errors.append(f'{companion_file} does not exists')

    return errors



def get_coordinate_to_cluster_links(df_overlapping):
    """
    Perform single-linkage cluster using dataframe resulted from bedtools instersection

    Parameters:
        df_overlapping (dataframe): dataframe returned from function getOverlappingCoordinates

    Returns:
        dictionary: coordinate as key and cluster_id as values

    """

    df_linkages = df_overlapping.copy(deep=True)

    round = 0
    clusters = []   # An array of arrays. Each sub-array represents a cluster

    while(1):   # process bed_object_f.intersect dataframe until it become empty
        if df_linkages.shape[0] == 0 or round >= 10000:
            break
        else:
            round+=1
            seed = list(df_linkages.loc[df_linkages.index[0], ['coordinates_x','coordinates_y']])
            seed = list(set(seed))
            cycle = 0

            while(1):   # expand seed array until it stops growing or intersect df is exhausted
                if df_linkages.shape[0] == 0:
                    break
                else:
                    cycle += 1
                    query_cluster = df_linkages[df_linkages['coordinates_x'].isin(seed) | df_linkages['coordinates_y'].isin(seed)]
                    df_linkages = df_linkages[~(df_linkages['coordinates_x'].isin(seed) | df_linkages['coordinates_y'].isin(seed))]

                    if query_cluster.shape[0] == 0 or df_linkages.shape[0] == 0:   # no room of seed expansion or intersect df is exhausted
                        clusters.append(seed)
                        break    
                    else:
                        seed = seed + list(query_cluster['coordinates_x']) + list(query_cluster['coordinates_y'])
                        seed = list(set(seed))

    # build coordinate-to-cluster dict for further annotation use
    cluster_count = len(clusters)
    clusters_id_base_number = 10 ** (len(str(cluster_count)))
    coordinate_to_cluster = {}
    for i, array in enumerate(clusters): 
        clusters_id_code_number = clusters_id_base_number + i + 1
        clusters_id_code_number = str(clusters_id_code_number)[1:]
        cluster_id = f'cluster_{clusters_id_code_number}'
        for coordinate in array:
            coordinate_to_cluster[coordinate] = cluster_id

    return coordinate_to_cluster



def add_cluster_annotations(df, radious=5):
    """
    Add columns to input dataframe that annotate cluster related features per row

    Parameters:
        df (dataframe): dataframe with cleavage site predictions
        radious (int): distance to between 2 positions when considering overlapping

    Returns:
        dataframe: a copy of input dataframe with more columns for cluster-level annotation

    """
    def _find_max_min_in_cell(cell_value):
        if isinstance(cell_value, list):
            return max(cell_value), min(cell_value)
        return cell_value, cell_value  # Return original value if not a list

    df = addCoordinateColumn(df.copy(deep=True))

    # intersect using bedtools
    cluster = getOverlappingCoordinates(df, df, radious)

    # extract linkage info from intersection
    coordinate_to_cluster_id = get_coordinate_to_cluster_links(cluster)

    # append cluster id to input dataframe
    df_temp = pd.DataFrame(data=coordinate_to_cluster_id.items(), columns=['coordinates', 'cluster_id'])
    df = df.merge(df_temp, on='coordinates')

    # extract and append merged coordinate to each cluster
    df_temp = df.groupby('cluster_id').agg({'coordinate(f)':lambda x:list(set(x.tolist()))}).reset_index()
    df_temp_r = df.groupby('cluster_id').agg({'coordinate(r)':lambda x:list(set(x.tolist()))}).reset_index()
    df_temp = df_temp.merge(df_temp_r, on='cluster_id') 
    df_temp['breakpoints'] = df_temp['coordinate(f)'] + df_temp['coordinate(r)']
    df_temp['end'], df_temp['start'] = zip(*df_temp['breakpoints'].apply(_find_max_min_in_cell))
    df_temp['cluster_coordinate'] = df_temp['start'].astype(str) + '-' + df_temp['end'].astype(str)
    df_temp = df_temp.drop(['coordinate(f)', 'coordinate(r)', 'breakpoints', 'start', 'end'], axis=1)
    df = df.merge(df_temp, on='cluster_id')
    df['cluster_coordinate'] = df['chr'] + ':' + df['cluster_coordinate']

    # extract and append score to each cluster
    df_temp = df.groupby('cluster_id').agg({'score':lambda x:list(set(x.tolist()))}).reset_index()
    df_temp = df_temp.rename(columns={'score':'all_scores'})
    df_temp['cluster_score'], df_temp['min_score'] = zip(*df_temp['all_scores'].apply(_find_max_min_in_cell))
    df_temp = df_temp.drop(['all_scores', 'min_score'], axis=1)
    df = df.merge(df_temp, on='cluster_id')

    # extract and append gamma list to each cluster
    df_temp = df.groupby('cluster_id').agg({'gamma':lambda x:list(set(x.tolist()))}).reset_index()
    df_temp = df_temp.rename(columns={'gamma':'cluster_gammas'})
    df_temp['cluster_gammas'] = df_temp['cluster_gammas'].astype(str)
    df = df.merge(df_temp, on='cluster_id')

    # extract and append a list of editable nucleotide presensce to each cluster
    if 'with_editable_nucleotide' in df.columns:
        df_temp = df.groupby('cluster_id').agg({'with_editable_nucleotide':lambda x:list(set(x.tolist()))}).reset_index()
        df_temp = df_temp.rename(columns={'with_editable_nucleotide':'cluster_editable_list'})
        df_temp['cluster_editable_list'] = df_temp['cluster_editable_list'].astype(str)
        df_temp['cluster_editability'] = 'False'
        df_temp.loc[df_temp['cluster_editable_list'].str.contains('True'), 'cluster_editability'] = 'True'
        df_temp = df_temp.drop('cluster_editable_list', axis=1)
        df = df.merge(df_temp, on='cluster_id')

    # extract and append cluster size to each cluster
    df_temp = df.copy(deep=True)
    df_temp['cluster_size'] = 1
    df_temp = df_temp.groupby('cluster_id')['cluster_size'].sum()
    df = df.merge(df_temp, on='cluster_id')

    # remove cluttered columns
    df = df.drop(['start', 'end', 'coordinates'], axis=1)

    return df



def get_highest_score_per_cluster(df):
    """
    Extract representative coordinate (row) with the highest score from each cluster 

    Parameters:
        df (dataframe): dataframe of cleavage sites with cluster annotation

    Returns:
        dataframe: each row represent a pair of coordinates overlapping each other
    """

    if 'cluster_id' in df.columns and 'score' in df.columns:
        df2 = df.copy(deep=True)
        df2 = df2.sort_values(['cluster_id', 'score'], ascending=[True, False] )
        df2 = df2.drop_duplicates('cluster_id')
        return df2
    else:
        return df
    


def remove_file_suffix(name):
    """
    Remove recognisable suffix from file name

    Parameters:
        name (string): file name in string

    Returns:
        string: file name without suffix
    """

    known_suffixes = ['csv', 'tsv', 'txt']
    array = name.split('.')
    if array[-1] in known_suffixes:
        array.pop()
    new_name = '.'.join(array)
    return new_name



def save_and_plot_result(df, output_table):
    """
    Saves a DataFrame to a CSV file and generates a Manhattan plot from the data.

    Parameters:
        df (pandas.DataFrame): The DataFrame containing the data to be saved and plotted.
        output_table (str): The file path for the output CSV. The function will save the plot 
                            with the same name but with a '.pdf' suffix.

    Returns:
        None
    """
    df.to_csv(output_table, index=False)
    manhattanPlot(df, f'{remove_file_suffix(output_table)}.pdf')
    formatedExcelSheet(df, f'{remove_file_suffix(output_table)}.xlsx')
