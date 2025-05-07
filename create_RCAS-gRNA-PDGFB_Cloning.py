import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import datetime

df = pd.DataFrame()
sgRNA_dic = {}


def convert_excel_to_pd(file_name: str):
    """
    Converts excel documents into a panda dataframe.
    """
    global df
    df = pd.read_excel(file_name)

    df.columns = df.columns.map(str)  # Convert column names to strings


def convert_pd_to_dic():
    """
    Converts pandas dataframes into a dictionary. 
    """
    global sgRNA_dic
    sgRNA_dic = pd.Series(
        df["sgRNA Sequence (5'-3')"].values, index=df["Gene"]).to_dict()
    # print(sgRNA_dic)


def modify_sequence(sequence: str, sgRNA: str, old_sgRNA: str):
    """
    Modify the string by deleting the existing 20 bp sgRNA sequence and replacing it
    with a new 20 bp sgRNA sequence.
    """
    index = sequence.index(old_sgRNA)
    new_sequence = sequence[:index] + sgRNA + \
        sequence[(index + len(old_sgRNA)):]
    return new_sequence, index


def check(RCAS_sequence: str, RCAS_features: list[SeqFeature], old_sgRNASequence: str):
    """
    Uses a previously made plasmid map to check if the algorithm is working correctly.
    """
    # Check the replacement sgRNA sequence is 20 bp
    if len(old_sgRNASequence) != 20:
        raise ValueError(
            "The inputted sgRNA sequence for the example map is the wrong length.")

    ARID1A_sgRNA = "TCAATCGATGATCTCCCCAT".lower()
    ARID1APDGFB_SequenceMap_filename = "RCASBP-Y DV_ARID1AgRNA_PGKpuro2APDGFB.fa"
    ARID1APDGFBSequenceMap_path = os.path.join(
        "supplementary_files", ARID1APDGFB_SequenceMap_filename)
    ARID1APDGFB_SequenceMap = SeqIO.parse(ARID1APDGFBSequenceMap_path, "fasta")
    ARID1APDGFB_SeqRecord = next(ARID1APDGFB_SequenceMap)
    ARID1APDGFB_sequence = str(ARID1APDGFB_SeqRecord.seq).lower()

    experimental_sequence, index = modify_sequence(
        RCAS_sequence, ARID1A_sgRNA, old_sgRNASequence)

    print(len(experimental_sequence))
    print(len(ARID1APDGFB_sequence))
    # if ARID1APDGFB_sequence == experimental_sequence:
    #     print("Algorithm passes check.")
    # else:
    #     raise ValueError("There is a problem with the algorithm.")


def main(gRNASequences_filename: str, RCAS_SequenceMap_filename: str, old_sgRNASequence: str):
    # convert excel to pd
    convert_excel_to_pd(gRNASequences_filename)
    convert_pd_to_dic()

    # import RCAS sequence
    RCASSequenceMap_path = os.path.join(
        "supplementary_files", RCAS_SequenceMap_filename)
    RCAS_SequenceMap = SeqIO.parse(RCASSequenceMap_path, "genbank")
    RCAS_SeqRecord = next(RCAS_SequenceMap)
    RCAS_sequence = str(RCAS_SeqRecord.seq).lower()
    RCAS_features = RCAS_SeqRecord.features

    # check with an existing plasmid map
    check(RCAS_sequence, RCAS_features, old_sgRNASequence)


if __name__ == "__main__":
    gRNASequences_filename = "gRNASequences.xlsx"
    RCAS_SequenceMap_filename = "0052 RCASBP-Y DV_ATMgRNA_PGKpuro2APDGFB.gb"
    old_sgRNASequence = "ataattcatgctgtcaccag".lower()
    main(gRNASequences_filename, RCAS_SequenceMap_filename,
         old_sgRNASequence)
