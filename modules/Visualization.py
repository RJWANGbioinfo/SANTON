import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype


def splitFileNameBySuffix(file_path):
    """
    Extract suffix from file path or file path
    
    Parameters:
        file_path (string): file path

    Returns:
        input filepath with suffix being removed
        suffix contained in input file path
    """
    array = file_path.split('.')
    return '.'.join(array[:-1]), array[-1]


def walkWithAbsolutePaths(directory):
    """
    Recursively walks through a directory, yielding absolute paths for directories and files.

    Parameters:
    directory (str): Path to the directory to traverse.

    Yields:
    tuple: (absolute_dirpath, dirnames, filenames) for each directory in the tree.
    """

    for dirpath, dirnames, filenames in os.walk(directory):
        absolute_dirpath = os.path.abspath(dirpath)
        yield absolute_dirpath, dirnames, filenames


def sortChromosomes(chromosomes):
    """
    Sort chromosomes by their digital suffix

    Parameters:
        chromosomes (list<string>): a list of chromosome names

    Returns:
        list<string>: a list of re-ordered chromosome names
    """

    if len(chromosomes) < 1:
        return chromosomes

    common_prefix = 'chr'
    if 'chromosome' in chromosomes[0]:
        common_prefix = 'chromosome'
    elif 'contig' in chromosomes[0]:
        common_prefix = 'contig'
    elif 'scaffold' in chromosomes[0]:
        common_prefix = 'scaffold'
    
    numbers = []
    letters = []
    for x in chromosomes:
        suffix = x.replace(common_prefix, '')
        if suffix.isdigit():
            numbers.append(int(suffix))
        else:
            letters.append(suffix)

    numbers.sort()
    letters.sort()
    sorted_suffix = numbers + letters
    return [f'{common_prefix}{x}' for x in sorted_suffix]



def manhattanPlot(df_input, output_file, log_scale=False):
    """
    Create Manhattan plot given a list of coordinates and their corresponding score

    Parameters:
        df_input (dataframe): a dataframe with columns 'score', 'chr', 'start' and 'end',
                                or with columns 'score', 'chr', ''coordinate(f)' and ''coordinate(r)'
        output_file (string): a path to output image file to save the plot
        log_scale (boolean): True (to plot score in log scale) or False (to plot score in a normal scale)

    Returns:
        None
    """
    df = df_input.copy(deep=True)
    if df.shape[0] == 0:
        return None

    sorter = list(df['chr'].unique())
    sorter = sortChromosomes(sorter)

    cat_size_order = CategoricalDtype(sorter, ordered=True)
    df=df[df['chr'].isin(sorter)]

    # add variables for downstream calculation of x-axis positions
    if 'start' not in df:
        df['start'] = df[['coordinate(f)','coordinate(r)']].min(axis=1).astype(int)
    if 'end' not in df:
        df['end'] = df[['coordinate(f)','coordinate(r)']].max(axis=1).astype(int)
    df['start']=df['start'].map(int)
    df['chr'] = df['chr'].astype(cat_size_order)

    # calculation of x-axis positions for each point and chromosome borders
    running_pos = 0
    cumulative_pos = []
    borders = [0]
    for chr, df_chr in df.groupby('chr'): 
        # collect x-axis position for datapoint
        cumulative_pos.append(df_chr['start'] + running_pos)
        # collect x-axis position for chromosome border
        chr_length = 0
        max_start_position = df_chr['start'].max()
        df_chr = df_chr.rename(columns={'Chr_length':'chr_length'})
        if 'chr_length' in df_chr.columns:
            chr_length = df_chr['chr_length'].max()        
        running_pos += max([max_start_position, chr_length])
        borders.append(running_pos)

    df['cumulative_pos'] = pd.concat(cumulative_pos)

    # plot data
    plot = sns.relplot(data=df, x='cumulative_pos', y='score', aspect=3.7, 
                       hue='chr', hue_order=sorter, palette = 'bright', legend=None) 
    
    if log_scale:
        plot.set(yscale="log")
    # add border between chromosomes
    for y in borders:
        ax1, = plot.axes[0]
        ax1.axvline(y, ls='--', color='silver')
    # add legends
    plot.fig.suptitle('Manhattan plot of scores')
    plot.ax.set_xlabel('Chromosome')
    plot.ax.set_ylabel('Score')
    # work on xticks
    chrom_df=df.groupby('chr')['cumulative_pos'].median()
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(sorter)
    plt.xticks(rotation=45, ha='center')
    plot.tight_layout()
    plot.savefig(output_file)



def formatedExcelSheet(df_input, output_file):
    """
    Writes a Pandas DataFrame to an Excel file with format for alignment visualization.
    
    Parameters:
        df_input (dataFrame): The DataFrame to be written.
        output_file (string): The path of the output Excel file.
    """

    def _max_character_count_per_column(df):
        """Calculate the max character count for each column."""
        return [max([len(str(s)) for s in df[col].values] + [len(col)]) for col in df.columns]


    # Create a Pandas Excel writer using XlsxWriter as the engine
    writer = pd.ExcelWriter(output_file, engine="xlsxwriter")

    # Convert the dataframe to an XlsxWriter Excel object
    df_input.to_excel(writer, sheet_name="Sheet1", index=False)

    # Get the xlsxwriter workbook and worksheet objects
    workbook = writer.book
    worksheet = writer.sheets["Sheet1"]

    # Set cell format
    format1 = workbook.add_format({"font_name": "Courier New", 'text_wrap': True})
    format1.set_align('vcenter')

    # Calculate max character count for column width adjustment
    column_widths = _max_character_count_per_column(df_input)

    for i, width in enumerate(column_widths):
        # Adjust column width based on content
        width = width / 2.5 if df_input.columns[i] == 'alignment' else width * 1.2
        worksheet.set_column(i, i, width, format1)

    # Save and close the writer
    writer._save()  # Use writer.save() for older pandas versions



def manhattanPlotOTscores(analysis_file, score_file, output_file, cutoff=0.1):
    """
    Generates a Manhattan plot for OT scores by merging analysis and score files, 
    filtering based on a cutoff, and extracting genomic coordinates.

    Parameters:
    analysis_file (str): Path to the analysis file with genomic coordinates and other annotations.
    score_file (str): Path to the score file containing OT scores.
    output_file (str): Path to save the Manhattan plot.
    cutoff (float, optional): Minimum OT score required for inclusion (default is 0.1).

    Returns:
    None: Saves the Manhattan plot to the specified output file.
    """

    # Check existance of input files
    invalid_file_paths = []
    for local_file in [analysis_file, score_file]:
        if not os.path.exists(local_file):
            invalid_file_paths.append(local_file)

    if len(invalid_file_paths):
        print('Essential file is not found for manhattan plot:', ','.join(invalid_file_paths))
        exit(0)

    # Read analysis file that contains OR region genomic coordinates, and chr length annotation
    df_analysis = pd.read_csv(analysis_file, sep='\t')
    df_analysis = df_analysis.rename(columns={'TemplateID':'Template_ID'})
    df_scores = pd.read_csv(score_file, sep='\t')
    
    # Merge data frames and filter row using cutoff
    annotation_column_names = ['Template_ID', 'Region']
    if 'chr_length' in df_analysis.columns:
        annotation_column_names = annotation_column_names + ['chr_length']
    df_scores = df_scores.merge(df_analysis[annotation_column_names], on='Template_ID', how='left')
    if not df_scores.shape[0]:
        print(f'No annotated OT regions for manhattan plot. Check these files for troubleshooting:', ','.join(invalid_file_paths))
        exit(0)

    df_scores = df_scores[df_scores['OT_score'] >= cutoff]
    if not df_scores.shape[0]:
        print(f'No OT region with score greater then {cutoff} for manhattan plot. Check these files for troubleshooting:', ','.join(invalid_file_paths))
        exit(0)


    # Create columns essential for data plotting
    df_scores[['chr', 'start-end', 'strand']] = df_scores['Region'].str.split('\:', regex=True, expand=True)
    df_scores[['start', 'end']] = df_scores['start-end'].str.split('\-', expand=True)
    df_scores = df_scores.drop('start-end', axis=1)
    df_scores = df_scores.rename(columns={'OT_score':'score'})

    # Plot data
    manhattanPlot(df_scores, output_file)



def manhattanPlotOTscoresLoop(analysis_file, analysis_outputs, cutoff=0.1):
    """
    Generates Manhattan plot for each OT score file found in the input dierectory

    Parameters:
    analysis_file (str): Path to the analysis file with genomic coordinates and other annotations.
    score_file (str): Path to the input directory that constains OT score files
    cutoff (float, optional): Minimum OT score required for inclusion in the plot (default is 0.1).

    Returns:
    None: Saves the Manhattan plot to the specified output file.
    """

    def _ot_score_tables(start_dir):
        """Retrun a list OT score files included in input directory """
        ot_tables = []
        for abs_dirpath, dirnames, filenames in walkWithAbsolutePaths(start_dir):
            for filename in filenames:
                absolute_filepath = os.path.join(abs_dirpath, filename)
                if absolute_filepath.endswith('.OT_Score.txt'):
                    ot_tables.append(absolute_filepath)
        return ot_tables
    
    score_files = _ot_score_tables(analysis_outputs)
    for file in score_files:
        file_no_suffix, suffix = splitFileNameBySuffix(file)
        manhattanPlotOTscores(analysis_file, file, f'{file_no_suffix}.pdf', cutoff)
