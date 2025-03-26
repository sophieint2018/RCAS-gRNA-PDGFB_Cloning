import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import glob
import copy

df = pd.DataFrame()
sgRNA_dic = {}
BbsIoverhang = "ggtcttcgagaagac"


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
    sgRNA_dic = pd.Series(
        df["sgRNA Sequence (5'-3')"].values, index=df["Gene"]).to_dict()
    # print(sgRNA_dic)


def adjust_feature_location(feature, offset, index):
    """
    Adjust the feature location based on the modification (insertions or deletions) in the sequence.
    """
    # Get the original start and end positions of the feature
    feature_start = -1
    feature_end = -1
    if (feature.location.start > index):
        feature_start = feature.location.start + offset
    else:
        feature_start = feature.location.start
    if (feature.location.end > index):
        feature_end = feature.location.end + offset
    else:
        feature_end = feature.location.end

    # Create a new FeatureLocation object with the updated positions
    updated_location = FeatureLocation(
        feature_start, feature_end, strand=feature.location.strand)

    # Return a new SeqFeature with the updated location
    return SeqFeature(location=updated_location, type=feature.type, qualifiers=feature.qualifiers)


def modify_sequence(sequence: str, overhang: str, sgRNA: str):
    index = sequence.index(overhang)
    new_sequence = sequence[:index] + sgRNA + \
        sequence[(index + len(overhang)):]
    return new_sequence, index


def modify_sequence_with_features(pDONR_sequence, BbsIoverhang, sgRNA, pDONR_features, key):
    """
    Modify the sequence and adjust the features based on the sequence modification.
    """
    # Modify the sequence (e.g., add the sgRNA and BbsI overhang)
    modified_sequence, index = modify_sequence(
        pDONR_sequence, BbsIoverhang, sgRNA)
    print("Modified the sequence.")

    # Calculate the difference in length
    sequence_length_diff = len(modified_sequence) - len(pDONR_sequence)

    # Update the feature locations
    updated_features = []
    for feature in pDONR_features:
        updated_feature = adjust_feature_location(
            feature, sequence_length_diff, index)
        updated_features.append(updated_feature)

    print("Updated features.")

    # Add in the gRNA feature
    sgRNA_featureLoc = FeatureLocation(index, index + 20, strand=1)
    sgRNA_feature = SeqFeature(
        location=sgRNA_featureLoc, type="misc_feature", qualifiers={"label": "{} sgRNA".format(key), "color": "blue", "Description": "Cloned {} sgRNA, mouse genome".format(key)})
    updated_features.append(sgRNA_feature)
    print("Added the {} sgRNA".format(key))

    return modified_sequence, updated_features


def check(pDONR_sequence, pDONR_features):
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


def main(gRNASequences_filename: str, pDONRSequenceMap_filename: str):
    # convert excel to pd
    convert_excel_to_pd(gRNASequences_filename)
    convert_pd_to_dic()

    # import pDONR sequence
    pDONRSequenceMap_path = os.path.join(
        "supplementary_files", pDONRSequenceMap_filename)
    pDONR_SequenceMap = SeqIO.parse(pDONRSequenceMap_path, "genbank")
    pDONR_SeqRecord = next(pDONR_SequenceMap)
    pDONR_sequence = str(pDONR_SeqRecord.seq).lower()
    pDONR_features = pDONR_SeqRecord.features

    # print the pDONR sequence
    # print(pDONR_sequence)

    # check with an existing plasmid map
    check(pDONR_sequence, pDONR_features)

    key = "test"
    # modify sequence
    sgRNA = "TCAATCGATGATCTCCCCAT"
    pDONRsgRNAPDGFB_Sequence, pDONRsgRNAPDGFB_modifiedfeatures = modify_sequence_with_features(
        pDONR_sequence, BbsIoverhang, sgRNA, pDONR_features, key)
    print("Modified sequence.")

    # create a SeqRecord object
    # pDONRsgRNAPDGFB_SequenceObject = copy.copy(pDONR_SequenceMap)
    # pDONRsgRNAPDGFB_SequenceObject.seq = Seq(pDONRsgRNAPDGFB_Sequence)
    pDONRsgRNAPDGFB_SequenceObject = Seq(pDONRsgRNAPDGFB_Sequence)
    pDONRsgRNAPDGFB_SequenceMap = SeqRecord(
        pDONRsgRNAPDGFB_SequenceObject, id=key, description=key, annotations={"molecule_type": "DNA"})
    pDONRsgRNAPDGFB_SequenceMap.features = pDONRsgRNAPDGFB_modifiedfeatures

    # export to genbank
    output_folder = "FASTA_output"
    pDONRsgRNAPDGFB_filename = "pDONR-U6-{}gRNA-PGKpuro2APDGFB.gb".format(
        key)
    output_filepath = os.path.join(output_folder, pDONRsgRNAPDGFB_filename)
    # with open(output_filepath, "w") as output_handle:
    #     SeqIO.write(pDONRsgRNAPDGFB_SequenceObject, output_handle, "genbank")
    with open(output_filepath, "w") as output_handle:
        SeqIO.write(pDONRsgRNAPDGFB_SequenceMap, output_handle, "genbank")
    print("Exported genbank file.")


if __name__ == "__main__":
    gRNASequences_filename = "gRNASequences.xlsx"
    pDONRSequenceMap_filename = "0036 pDONR-U6gRNA-PGKpuro2APDGFB.gb"
    main(gRNASequences_filename, pDONRSequenceMap_filename)
