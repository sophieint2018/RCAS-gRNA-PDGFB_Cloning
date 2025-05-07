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


def check():
    """
    Uses a previously made plasmid map to check if the algorithm is working correctly.
    """
    ARID1A_sgRNA = "TCAATCGATGATCTCCCCAT"
    ARID1APDGFB_SequenceMap_filename = "pDONR-U6-ARID1AgRNA-PGKpuro2APDGFB.fa"
    ARID1APDGFBSequenceMap_path = os.path.join(
        "supplementary_files", ARID1APDGFB_SequenceMap_filename)
    ARID1APDGFB_SequenceMap = SeqIO.parse(ARID1APDGFBSequenceMap_path, "fasta")
    ARID1APDGFB_SeqRecord = next(ARID1APDGFB_SequenceMap)
    ARID1APDGFB_sequence = str(ARID1APDGFB_SeqRecord.seq)

    experimental_sequence, experimental_features = modify_sequence_with_features(
        pDONR_sequence, BbsIoverhang, ARID1A_sgRNA, pDONR_features, "ARID1A Check")

    if ARID1APDGFB_sequence.lower() == experimental_sequence.lower():
        print("Algorithm passes check.")
    else:
        raise ValueError("There is a problem with the algorithm.")


def main(gRNASequences_filename: str, RCAS_SequenceMap_filename: str):
    print("hello")


if __name__ == "__main__":
    gRNASequences_filename = "gRNASequences.xlsx"
    RCAS_SequenceMap_filename = "0052 RCASBP-Y DV_ATMgRNA_PGKpuro2APDGFB.fa"
    main(gRNASequences_filename, RCAS_SequenceMap_filename)
