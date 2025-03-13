import os
from Bio import SeqIO
from Bio.Seq import Seq
import pybedtools
from Bio import BiopythonWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    from Bio import pairwise2


def reverseComplement(sequence):
    """Reverse complement a given nucleotide sequence

    Parameters:
        sequence (str): The input nucleotide sequence (e.g., 'ATCG').

    Returns:
        str: The reverse complement of the input sequence.
    """

    return str(Seq(sequence).reverse_complement())



def alignmentAnnotationString(query, reference):
    """
    Bulid an annotation string for two aligned sequences (handle degenerative nucleotide N and R/Y)

    Parameters:
        qryseq (string): aligned query sequence
        refseq (string): aligned reference sequence

    Returns
        string: string denote (mis)matches and gapes in alginment, such as "||||~||~~|"
    """
    query = query.upper()
    reference = reference.upper()

    aln_string=""
    for c in range(len(reference)):
        refbase = reference[c]
        qrybase = query[c]
        if (refbase=="-") or (qrybase=="-"):
            aln_string=aln_string+"~"
        elif refbase == qrybase or 'N' in [refbase, qrybase]:
            aln_string=aln_string+"|"
        elif qrybase == 'R' and refbase in ["A", "G"]:
            aln_string=aln_string+"|"
        elif qrybase == 'Y' and refbase in ["C", "T"]:
            aln_string=aln_string+"|"
        else:
            aln_string=aln_string+"."

    return aln_string



def adjustPositionConsideringGap(str, position):
    """
    Convert a positon in ungapped sequence to its correspondence in gapped version of the same sequence

    Parameters:
        str (string): a sequence with (or without) gaps
        position (int): a positon in the input sequence

    Returns
        int: a matched poistion in ungapped version of input sequence
    """
    
    adjusted_pos = position
    while 1 > 0:
        substr = str[0:adjusted_pos]
        gap = substr.count('-')
        char_count = len(substr) - gap
        
        if char_count == position:
            break
        else:
            adjusted_pos = position + gap
            
        if adjusted_pos > len(str):
            adjusted_pos = len(str)
            break

    return adjusted_pos



def sliceSequenceFlankingCutSites(row, flanking_size_ref, pam_location):
    """
    Extract cut site flanking sequences from gRNA-ref alignment
        This function is to be used in dataframe.apply(func) format. The dataframe has at least
        two columns 'alignment' and 'gRNA_strand'
        The reference in aligment should be built as  5'flanking seq + sequ between 2 cut sites + 3'flanking seq 
        (flanking seq are of the same length as defined by input argument)
    
    Parameters:
        row (pandas series): a row in a dataframe
        flanking_size_ref (int): length of flanking sequence on both sides of cut site in reference sequence
        pam_location (string): either 'right' or 'left'

    Returns:
        string: pam-side  grna sequence flanking cut site
        string: 5'-partial reference sequence flanking cut site
        string: 3'-partial reference sequence flanking cut site
    """
    # sanity check
    if not pam_location or not flanking_size_ref:
        return '', '', ''
    if 'alignment' not in row.index or 'gRNA_strand' not in row.index:
        return '', '', ''
    
    # extract sequence
    alignment_str = row['alignment']
    gRNA_strand = row['gRNA_strand']
    [ref, aln, grna] = alignment_str.split("\n")

    # slice reference sequence flanking cut sites
    ungapped_ref = ref.replace('-', '')
    ref_5flanking_cut = ungapped_ref[:flanking_size_ref]
    ref_3flanking_cut = ungapped_ref[-flanking_size_ref:]

    # slice sequence
    cut1_index_ungapped = flanking_size_ref
    cut2_index_ungapped = len(ref.replace('-', '')) - flanking_size_ref
    cut1_index_gapped = adjustPositionConsideringGap(ref, cut1_index_ungapped)
    cut2_index_gapped = adjustPositionConsideringGap(ref, cut2_index_ungapped)
    guide_flanking5 = grna[:cut1_index_gapped].replace('-', '')
    guide_flanking3 = grna[cut2_index_gapped:len(grna)].replace('-', '')

    # pick the right flanking sequence dependent on pam location
    if pam_location == 'right':
        pam_side = guide_flanking3 if gRNA_strand == 'forward' else reverseComplement(guide_flanking5)
    elif pam_location == 'left':
        pam_side = guide_flanking5 if gRNA_strand == 'forward' else reverseComplement(guide_flanking3)

    return pam_side, ref_5flanking_cut, ref_3flanking_cut



def alignGuideAndReference(ref, qry):
    """
    Align two sequences in fasta format (genomic reference and gRNA)

    Parameters:
        qryseq (string): aligned query sequence
        refseq (string): aligned reference sequence

    Returns:
        full_alignment (string): string containing aligned query and reference
        mmcount (int): number of mismatches and gaps in the alignemnt
        strand (string): the orientation of query relative to reference sequence
    """
    qryr=reverseComplement(qry)
    

    aln=None
    mmcount=-99
    full_alignment = ''

    # do pairwise alignment
    aln_forward=pairwise2.align.localms(ref, qry , 2.5, -0.5, -3, -2.5)
    aln_reverse=pairwise2.align.localms(ref, qryr, 2.5, -0.5, -3, -2.5)

    fscore = aln_forward[0][2] if len(aln_forward) > 0 else 0
    rscore = aln_reverse[0][2] if len(aln_reverse) > 0 else 0
    score = fscore if fscore > rscore else rscore
    strand = 'forward' if fscore > rscore else 'reverse'
    aln = aln_forward if fscore > rscore else aln_reverse

    if score>=0:
        refstr=aln[0][0]
        alnstr=aln[0][1]
        astr=alignmentAnnotationString(alnstr, refstr)
        full_alignment=refstr+"\n"+astr+"\n"+alnstr
        string=alnstr
        string2=astr
        starti=next((i for i, x in enumerate(string) if x!='-'), None)
        endi=max(i for i, x in enumerate(string) if x!='-')
        onstring2=string2[int(starti):int(endi)+1]
        mmcount=sum(x != "|" for x in onstring2)
    
    return full_alignment, mmcount, strand



def retrieveReferenceFlankingSequences(df, genome_reference, flanking_length=25):
    """
    Given a list coordinate pairs (forward & reverse) and reference sequence, retrieve sequence centered on 
    the coordinates with specified length on both sides

    Parameters:
        df (dataframe): a dataframe with at two columns 'coordinate(f)' and 'coordinate(r)' 
        genome_reference (string): path to a local file for reference genome in fasta format
        flanking_length (integer): length of flanking sequence to retrieve on both sides

    Returns
        dataframe: dataframe with columns 'coordinate(f)' ,'coordinate(r)' and 'region_sequence'
    """

    if not (genome_reference and os.path.exists(genome_reference)):
        return pd.DataFrame()

    # create BED from coordinate columns in input dataframe
    dfbed = df[['coordinate(f)' ,'coordinate(r)']].copy(deep=True)
    dfbed['chr'] = dfbed['coordinate(r)'].replace(r':.+$', '', regex=True)
    dfbed['start'] = dfbed['coordinate(f)'].replace(r'^[^\:]*:', '', regex=True).astype(int) - (1 + flanking_length)   # have to minus to get the first basepair 
    dfbed.loc[dfbed['start']<0, 'start'] = 0 
    dfbed['end'] = dfbed['coordinate(r)'].replace(r'^[^\:]*:', '', regex=True).astype(int) +  flanking_length
    dfbed['region_coordinate'] = dfbed['chr'] + ':' + dfbed['start'].astype(str) + '-' + dfbed['end'].astype(str)
    dfbed = dfbed[['chr', 'start', 'end', 'region_coordinate'] + ['coordinate(f)', 'coordinate(r)']]
    dfbed = dfbed.set_index('region_coordinate')

    # retrieve sequence of BED from input genome reference
    dfbed['region_sequence'] = ''
    surrounding_cuts_bedtool = pybedtools.BedTool.from_dataframe(dfbed)
    surrounding_cuts_seq = surrounding_cuts_bedtool.sequence(fi=genome_reference)
    surrounding_cuts_seq_blob = SeqIO.parse(open(surrounding_cuts_seq.seqfn), 'fasta')
    for fastablob in surrounding_cuts_seq_blob:
        surrounding_cuts_coordinate, surrounding_cuts_sequence = fastablob.id, str(fastablob.seq)
        dfbed.loc[surrounding_cuts_coordinate, 'region_sequence'] = surrounding_cuts_sequence

    # add more feature columns
    dfbed = dfbed.reset_index()
    dfbed = dfbed.drop(['chr','start','end'], axis=1)
    # dfbed['flanking_5'] = dfbed['region_sequence'].str.slice(0, flanking_length)
    # dfbed['overlap_nt'] = dfbed['region_sequence'].str.slice(flanking_length, 0-flanking_length)
    # dfbed['flanking_3'] = dfbed['region_sequence'].str.slice(0-flanking_length)
    # dfbed['edit_site_5'] = dfbed['flanking_5'].str.slice(-2, -1)
    # dfbed['edit_site_3'] = dfbed['flanking_3'].str.slice(1, 2)

    return dfbed



def addSequenceAlignmentFeatures(df, genome_reference, gRNA, flanking_length=25, pam_location='', enzyme=''):
    """
    Add alignment feature columns to input dataframe

    Parameters:
        df (dataframe): a dataframe with at two columns 'coordinate(f)' and 'coordinate(r)' 
        genome_reference (string): path to a local file for reference genome in fasta format
        gRNA (string): guide RNA sequence
        flanking_length (integer): length of flanking sequence to retrieve on both sides

    Returns
        dataframe: a copy of input dataframe with additional columns such as 'alignment, 'mm+gap' and 'gRNA_strand'
    """

    if genome_reference and os.path.exists(genome_reference) and len(gRNA) > 0:
        enzyme = enzyme.upper()
        # retrive genomic sequence for cuts
        df_features = retrieveReferenceFlankingSequences(df, genome_reference, flanking_length)
        if df_features.shape[0] == 0: return df
        
        # alignment genomic sequence with gRNA and add relevant features
        df_features['alignment'], df_features['mm+gap'], df_features['gRNA_strand'] = zip(*df_features['region_sequence'].apply(alignGuideAndReference, args=[gRNA]))
        df_features['pamside_flanking_grna_sequence'], df_features['ref_seq_5flanking_cut'], df_features['ref_seq_3flanking_cut'] = zip(*df_features.apply(sliceSequenceFlankingCutSites, args=[flanking_length, pam_location], axis=1))
        if enzyme in ['ABE', 'CBE']:
            df_features['edit_window'] = df_features['ref_seq_5flanking_cut']
            df_features.loc[df_features['gRNA_strand']=='reverse', 'edit_window'] = df_features.loc[df_features['gRNA_strand']=='reverse', 'ref_seq_3flanking_cut'].apply(reverseComplement)
            df_features['edit_window'] = df_features['edit_window'].str.slice(-5)
            # test existence of editable nucleotide
            editable_nucleotide = 'A|a' if enzyme == 'ABE' else 'C|c'
            df_features['with_editable_nucleotide'] = df_features['edit_window'].str.contains(editable_nucleotide).astype(str)

        df_features = df_features.drop(['region_coordinate', 'region_sequence', 'pamside_flanking_grna_sequence', 'ref_seq_5flanking_cut', 'ref_seq_3flanking_cut'], axis=1)

        # append features to input dataframe
        return df.merge(df_features, on=['coordinate(f)', 'coordinate(r)'])
    else:
        return df