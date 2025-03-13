import pybedtools
import pandas as pd
from IPython.display import display



def addCoordinateColumn(input_df):
    """
    Add a column to input dataframe for string encoding chromosome and two cut positions

    Parameters:
        input_df (dataframe): a dataframe with at least 3 columns 'chr', 'coordinate(f)', and 'coordinate(r)'

    Returns:
        dataframe: A copy of input dataframe with additional columns 'coordinates', 'start', and 'end'
    """
    df = input_df.copy(deep=True)
    df['start'] = df[['coordinate(f)','coordinate(r)']].min(axis=1).astype(int)
    df['end'] = df[['coordinate(f)','coordinate(r)']].max(axis=1).astype(int)
    df['coordinates'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    return df



def getOverlappingCoordinates(df_cooridante1, df_cooridante2, radius=0):
    """
    Find overlapping between two list of coordinates

    Parameters:
        df_cooridante1 (dataframe): dataframe with at least three columns: 'chr', 'start', and 'end'
        df_cooridante2 (dataframe): dataframe with at least three columns: 'chr', 'start', and 'end'
        radius (int): allowed distance between two coordinates. Overlapping is considered if distance 
                      is not greater than this cutoff

    Returns:
        dataframe: each row represent a pair of coordinates overlapping each other
    """
    essential_columns = ['chr', 'start', 'end', 'coordinates']
    df1 = addCoordinateColumn(df_cooridante1)
    df2 = addCoordinateColumn(df_cooridante2)
    df1_backup = df1.copy(deep=True)
    df2_backup = df2.copy(deep=True)

    # prepare 'bed' for intersection
    df1['end'] = df1['end'] + 1
    df2['end'] = df2['end'] + 1
    if radius >0:
        df2['start'] = df2['start'] - radius
        df2.loc[df2['start'] < 0, 'start'] = 0
        df2['end'] = df2['end'] + radius

    dfbed1 = df1[essential_columns]
    dfbed2 = df2[essential_columns]

    # intersect two sets of bed regions
    bed_object_f = pybedtools.BedTool.from_dataframe(dfbed1)
    bed_object_r = pybedtools.BedTool.from_dataframe(dfbed2)
    overlaps = bed_object_f.intersect(bed_object_r, wo=True).to_dataframe()

    if overlaps.shape[0] > 0:
        overlaps = overlaps.merge(dfbed1, left_on='name', right_on ='coordinates', how='inner')
        overlaps = overlaps.merge(dfbed2, left_on='thickEnd', right_on ='coordinates', how='inner')
        overlaps['overlap_with_background'] = overlaps['itemRgb']
        return overlaps.iloc[:, 9:]
    else:
        return pd.DataFrame()
    


def removeBackgroundSites(df_foreground, df_background, radius=5):
    """
    Find overlapping between two list coordinates and remove the overlapping part from the first

    Parameters:
        df_foreground (dataframe): a dataframe with overlapping coordinates to be removed from
        df_cooridante2 (dataframe): a dataframe to identify overlapping with the first dataframe
        radius (int): allowed distance between two coordinates. Overlapping is considered if distance 
                      is not greater than this cutoff

    Returns:
        dataframe: modified df_foreground with overlapping coordinates removed
        dataframe: a dataframe contains pairs of foreground-background coordinates overlapping each other
    """
    df_bonafide = pd.DataFrame()
    df_artifacts = pd.DataFrame()
    foreground_size = df_foreground.shape[0]
    background_size = df_background.shape[0]

    # return input dataframe in simple scenarios
    if foreground_size == 0: 
        return df_bonafide, df_artifacts
    if background_size == 0:
        return df_foreground, df_artifacts

    # excise interval intersection when neither foreground nor background is empty
    df_shared = getOverlappingCoordinates(df_foreground, df_background, radius)
    if df_shared.shape[0] == 0:
        return df_foreground, df_artifacts

    # post-invertval intersection manipulation
    df_foreground = addCoordinateColumn(df_foreground)
    df_bonafide = df_foreground[~df_foreground['coordinates'].isin(df_shared['coordinates_x'])]
    df_artifacts = df_foreground[df_foreground['coordinates'].isin(df_shared['coordinates_x'])]
    df_artifacts = df_artifacts.merge(df_shared[['coordinates_x', 'coordinates_y', 'overlap_with_background']], left_on='coordinates', right_on='coordinates_x', how='inner')
    df_artifacts = df_artifacts.drop('coordinates_x', axis=1)
    df_artifacts = df_artifacts.rename(columns={'coordinates_y':'overlapping_background_coordinate'})
    return df_bonafide, df_artifacts



def removeSitesFoundInControl(edited, control):
    """
    Find shared coordinates between two dataframe and remove the shared sites from the first

    Parameters:
        edited (dataframe): a dataframe with coordinates to be cleaned through
        control (dataframe): a dataframe with coordinates as background noise

    Returns:
        dataframe: modified edited with shared coordinates being removed
    """
    df_edited = edited.copy(deep=True)
    df_control = control.copy(deep=True)
    if df_control.shape[0] > 0:
        df_edited = addCoordinateColumn(df_edited)
        df_control = addCoordinateColumn(df_control)
        df_true_edited = df_edited[~df_edited['coordinates'].isin(df_control['coordinates'])]
        return df_true_edited
    else:
        return df_edited