import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import datetime
import copy

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
    modified_sequence = sequence[:index] + sgRNA + \
        sequence[(index + len(old_sgRNA)):]
    return modified_sequence.lower(), index


def creat_new_sequence_map(sequence: str, sgRNA: str, old_sgRNA: str, RCAS_features: list[SeqFeature], sgRNA_identifier: str):
    """
    Creates a new sequence map object with the sequence containing the appropriate sgRNA sequence and labeled feature.
    """
    # Modify the sequence (e.g. Change the sgRNA)
    modified_sequence, index = modify_sequence(sequence, sgRNA, old_sgRNA)
    print("Modified the sequence for {} plasmid map.".format(sgRNA_identifier))

    # Add in the gRNA feature
    sgRNA_featureLoc = FeatureLocation(index, index + 20, strand=1)
    sgRNA_feature = SeqFeature(
        location=sgRNA_featureLoc, type="misc_feature", qualifiers={"label": "{} sgRNA".format(sgRNA_identifier),
                                                                    "color": "blue",
                                                                    "Description": "Cloned {} sgRNA, mouse genome".format(sgRNA_identifier), })
    RCAS_features.append(sgRNA_feature)
    print("Added the {} sgRNA feature to the plasmid map.".format(sgRNA_identifier))

    return modified_sequence, RCAS_features


def check(RCAS_sequence: str, old_sgRNASequence: str):
    """
    Uses a previously made plasmid map (CDKN2A sgRNA) to check if the
    algorithm is working correctly.
    """
    # Check the replacement sgRNA sequence is 20 bp
    if len(old_sgRNASequence) != 20:
        raise ValueError(
            "The inputted sgRNA sequence for the example map is the wrong length.")

    verif_sgRNA = "tggtgaagttcgtgcgatcc".lower()
    verif_SequenceMap_filename = "0053 RCASBP-Y DV_CDKN2AgRNA_PGKpuro2APDGFB.fa"
    verif_SequenceMap_path = os.path.join(
        "supplementary_files", verif_SequenceMap_filename)
    verif_SequenceMap = SeqIO.parse(verif_SequenceMap_path, "fasta")
    SeqRecord = next(verif_SequenceMap)
    verif_Sequence = str(SeqRecord.seq).lower()

    experimental_sequence, features = modify_sequence(
        RCAS_sequence, verif_sgRNA, old_sgRNASequence)

    if verif_Sequence == experimental_sequence:
        print("Algorithm passes check.")
    else:
        raise ValueError("There is a problem with the algorithm.")


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
    check(RCAS_sequence, old_sgRNASequence)

    # iterates through the dictionary to make plasmid maps for all the sgRNA sequences provided.
    for sgRNA_identifier in sgRNA_dic.keys():
        # Modify sequence
        sgRNA = sgRNA_dic[sgRNA_identifier]

        # Create the new sequence map
        RCASfeatures_copy = copy.deepcopy(RCAS_features)
        RCASsgRNAPDGFB_Sequence, RCASsgRNAPDGFB_modifiedfeatures = creat_new_sequence_map(
            RCAS_sequence, sgRNA, old_sgRNASequence, RCASfeatures_copy, sgRNA_identifier)
        print("Modified the sequence for {} plasmid map.".format(sgRNA_identifier))

        # create a SeqRecord object
        RCASsgRNAPDGFB_SequenceObject = Seq(RCASsgRNAPDGFB_Sequence)
        RCASsgRNAPDGFB_SequenceMap = SeqRecord(
            RCASsgRNAPDGFB_SequenceObject, id=sgRNA_identifier, description=sgRNA_identifier, annotations={"molecule_type": "DNA",
                                                                                                           "topology": "circular",
                                                                                                           "date": datetime.datetime.now().strftime("%Y-%m-%d")})
        RCASsgRNAPDGFB_SequenceMap.features = RCASsgRNAPDGFB_modifiedfeatures

        # export to genbank
        output_folder = "FASTA_output"
        RCASsgRNAPDGFB_filename = "RCASBP-Y DV_{}gRNA_PGKpuro2APDGFB.gb".format(
            sgRNA_identifier)
        output_filepath = os.path.join(output_folder, RCASsgRNAPDGFB_filename)
        with open(output_filepath, "w") as output_handle:
            SeqIO.write(RCASsgRNAPDGFB_SequenceMap, output_handle, "genbank")
        print("Exported genbank file for {}.\n".format(sgRNA_identifier))


if __name__ == "__main__":
    gRNASequences_filename = "gRNASequences.xlsx"
    RCAS_SequenceMap_filename = "0052 RCASBP-Y DV_ATMgRNA_PGKpuro2APDGFB.gb"
    old_sgRNASequence = "ataattcatgctgtcaccag".lower()
    main(gRNASequences_filename, RCAS_SequenceMap_filename,
         old_sgRNASequence)
