"""
Calls all south africa sequences 'C' subtype
"""

import pandas as pd
import pickle
import os

def read_fasta(fasta_file):
    seqs = {}
    max_len = 0
    with open(fasta_file, "r") as f:
        last_id = ""
        for i, line in enumerate(f):
            if i % 2 == 0:
                last_id = line.strip()[1:]  # removing ">"
                seqs[last_id] = ""
            else:
                seqs[last_id] = get_seq_name(line.strip())
                if len(line.strip()) > max_len:
                    max_len = len(line.strip())
    return seqs

def write_csv(seqs_dict, output_file, subtype_to_label):
    i = 0
    with open(output_file, "w") as f:
        # f.write("sequence,label\n")
        for seqid, seq in seqs_dict.items():
            subtype = get_subtype(seqid)
            if subtype == "unknown": continue
            try:
                label = subtype_to_label[subtype]
            except KeyError:
                #China
                if "01_AE" in subtype:
                    print(f"Subtype {subtype} not found in label map, using 01_AE instead")
                    subtype = "01_AE"
                    label = subtype_to_label[subtype]
                #Bulgaria
                elif "F1B" in subtype:
                    print(f"Subtype {subtype} not found in label map, using 12_BF1 instead")
                    subtype = "12_BF1"
                    label = subtype_to_label[subtype]
                elif "BF1" in subtype:
                    print(f"Subtype {subtype} not found in label map, using 17_BF1 instead")
                    subtype = "17_BF1"
                    label = subtype_to_label[subtype]
                elif "A1D" in subtype:
                    print(f"Subtype {subtype} not found in label map, using 35_A1D instead")
                    subtype = "35_A1D"
                    label = subtype_to_label[subtype]
                elif "A1C" in subtype:
                    print(f"Subtype {subtype} not found in label map, using A1C URF instead")
                    subtype = "A1C URF"
                    label = subtype_to_label[subtype]
                elif "BG" in subtype:
                    print(f"Subtype {subtype} not found in label map, using 14_BG instead")
                    subtype = "14_BG"
                    label = subtype_to_label[subtype]
                elif "BC" in subtype:
                    print(f"Subtype {subtype} not found in label map, using 31_BC instead")
                    subtype = "31_BC"
                    label = subtype_to_label[subtype]
                else:
                    print(f"Subtype {subtype} not found in label map, using CPZ instead")
                    subtype = "CPZ"
                    label = subtype_to_label[subtype]
                # else:
                #     print(f"Subtype {subtype} not found in label map")
                #     continue
            f.write(f"{seq.upper()},{label},{seqid}\n")
            i += 1
    print("Wrote", i, "sequences to", output_file)

def get_subtype(seq):
    # meta has columns "Accn ID only" and "Final subtype call"
    if seq in meta["Accn ID only"].values:
        subtype = meta[meta["Accn ID only"] == seq]["Final subtype call"].values[0]
    elif "africa" in seq.lower():
        subtype = "C"
    else:
        # raise ValueError(f"Sequence {seq} not found in metadata")
        print(f"Sequence {seq} not found in metadata")
        return "unknown"
    return subtype.strip()


"""

CUSTOM DEFINED FOR EACH DATASET

"""
def get_seq_name(seq_name):
    return seq_name.split(".")[-1]

def get_subtype(seq_name):
    return seq_name.split(".")[0]

with open(r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\unified_map_subtype_to_label.pkl", "rb") as f:
    subtype_to_label = pickle.load(f)

print("LABEL MAP:", subtype_to_label)


# meta_path= '../subtypes.csv'
# meta_path = "../../raw/HIV1_NXT_2021_pol_DNA_2809taxa_BS_notaln_1-2809_REGA_final_subtypes_DRMs.csv"
# seqs_path = '../../raw/HIV1_NXT_2021_pol_DNA_2809taxa_BS_notaln_subtype_not_in_name+SA_1148taxa.fasta'
# seqs_path = '../../raw/HIV-1_Subtype_REF_2020_genome_DNA_LANL_NOTaln_459taxa.fas'
# seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\relabeled\2253_filtered_NXTSA_with_subtypes_RELABELED.fasta"
# seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\all_lanl_full_genomes\test_generated\final_output.fasta"
# seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\train_4293\test_generated\final_output.fasta"
# seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\validation_datasets\LANL_full_genomes_minus_4293.fasta"
# seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol\test.fasta"
# seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol\test_generated\final_output.fasta"

seqs_path = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\new_datasets\Bulgaria\BulgariaValidation.fasta"

# out_dir = "hiv1_full_genome_not_align"
# out_dir = "filteredNXTSA"
# out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\all_lanl_full_genomes\test_generated"
# out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\train_4293"
# out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\validation_datasets"
# out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol"
out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\validation_datasets\Bulgaria"

# out_path= f'{out_dir}/LANL_minus_4293_curated.csv'
# out_path = f"{out_dir}/train.csv"
# out_path = f"{out_dir}/test.csv"
out_path = f"{out_dir}/BulgariaValidation.csv"

# if not outpathexists create it
if not os.path.exists(f"{out_dir}"):
    os.makedirs(f"{out_dir}")

# meta = pd.read_csv(meta_path, header=1)
seqs = read_fasta(seqs_path)

print("SEQS:", len(seqs))


# subtypes_not_in_label_map = set([i.strip() for i in meta["Final subtype call"].unique()]) - set(subtype_to_label.keys())
# # add these subtypes to the label map
# for subtype in subtypes_not_in_label_map:
#     subtype_to_label[subtype] = len(subtype_to_label)

# print("new label map:", subtype_to_label, len(subtype_to_label))

write_csv(seqs, out_path, subtype_to_label)

# with open(f"map_subtype_to_label_{len(subtype_to_label)}.pkl", "wb") as f:
#     pickle.dump(subtype_to_label, f)

