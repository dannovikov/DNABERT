from sklearn.model_selection import train_test_split as sklearn_train_test_split
import random
import pandas as pd
import pickle

RAND_SEED = 1
TEST_SIZE = 0.15

# NUM_CLASSES = 200
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\6105_pol_noalign_plus_last_additional_seqs_subtyped.fasta"
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\test_generated\7584_seqs_noalign_aug_shifted.fasta"
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\combined_full_and_pol\raw\combined.fasta"
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\relabeled\6105_pol_noalign_plus_last_additional_seqs_subtyped_NORELABLENEEDED.fasta"
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation\6105_pol_noalign\generated\augmin3_6131_pol.fasta"
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation\459_full_genome_noalign\generated\augmin3_497_full_genome.fasta"


# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\new_datasets\HIV-1_FULL_GENOME_ALL_LANL_2line_clean_known.fasta"
# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation\all_lanl_full_genome\test_generated\final_output.fasta"

# raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation\4293_from_5225_seqs_no_unknown_no_gaps\base_generated\final_output.fasta"
raw_sequences = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation\6105_pol_noalign\generated\augmin5_6131_pol.fasta"

# out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\all_lanl_full_genomes"
# out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\train_4293"

out_dir = r"C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol"

train_output = f"{out_dir}\\train.csv"
test_output = f"{out_dir}\\test.csv"

def get_subtype(seqid):
    return seqid.split(".")[0]

def strat_train_test_split(seqs_dict):
    # stratified train test split based on subtypes
    seqs = list(seqs_dict.keys())

    subtypes =  [get_subtype(seq) for seq in seqs]
    single_types = set([subtype for subtype in subtypes if subtypes.count(subtype) == 1])
    less_than_five_types = set([subtype for subtype in subtypes if 1 < subtypes.count(subtype) < 5])

    filtered_seqs = [seq for seq in seqs if get_subtype(seq) not in single_types and get_subtype(seq) not in less_than_five_types]
    removed_seqs =  [seq for seq in seqs if seq not in filtered_seqs]

    filtered_subtypes = [get_subtype(seq) for seq in filtered_seqs]

    train_seqs, test_seqs = sklearn_train_test_split(filtered_seqs, test_size=TEST_SIZE, stratify=filtered_subtypes, random_state=RAND_SEED)

    train_seqs_dict = {seq: seqs_dict[seq] for seq in train_seqs}
    test_seqs_dict = {seq: seqs_dict[seq] for seq in test_seqs}

    train_types = [get_subtype(seq) for seq in train_seqs]
    test_types = [get_subtype(seq) for seq in test_seqs]

    for seq in removed_seqs:
        # if the seq is from a single subtype, add it to the train set
        if get_subtype(seq) in single_types:
            train_seqs_dict[seq] = seqs_dict[seq]
            print(f"Adding {seq} to train set of type {get_subtype(seq)} because it is from a single type")
        # if the seq is from a subtype with less than 5 sequences, add 1 to the test set and the rest to the train set
        elif get_subtype(seq) in less_than_five_types:
            if get_subtype(seq) not in test_types:
                test_seqs_dict[seq] = seqs_dict[seq]
                test_types.append(get_subtype(seq))
                print(f"Adding {seq} to test set of type {get_subtype(seq)} because it is from a type with less than 5 sequences")
            else:
                train_seqs_dict[seq] = seqs_dict[seq]
                print(f"Adding {seq} to train set of type {get_subtype(seq)} because it is from a type with less than 5 sequences")


    print(f"Number of subtypes in train set: {len(set([get_subtype(seq) for seq in train_seqs_dict.keys()]))}")
    print(f"Number of subtypes in test set: {len(set([get_subtype(seq) for seq in test_seqs_dict.keys()]))}")

    # create label maps for the train and test sets
    sorted_subtypes = sorted(set(subtypes))
    map_subtype_to_label = {subtype: i for i, subtype in enumerate(sorted_subtypes)}
    test_map_subtype_to_label = {subtype: i for subtype, i in map_subtype_to_label.items() if subtype in test_types}

    return train_seqs_dict, test_seqs_dict, map_subtype_to_label, test_map_subtype_to_label


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
                seqs[last_id] = line.strip()
                if len(line.strip()) > max_len:
                    max_len = len(line.strip())
    return seqs

def write_csv(seqs_dict, output_file, subtype_to_label):
    with open(output_file, "w") as f:
        for seqid, seq in seqs_dict.items():
            subtype = get_subtype(seqid)
            label = subtype_to_label[subtype]
            f.write(f"{seq.upper()},{label}\n")

def write_fasta(seqs_dict, output_file, subtype_to_label):
    with open(output_file, "w") as f:
        for seqid, seq in seqs_dict.items():
            # subtype = get_subtype(seqid)
            # name = f">{subtype}.{seqid}"
            name = f">{seqid}"
            f.write(name + "\n")
            f.write(seq + "\n")
            

# with open(fr"{out_dir}\..\map_subtype_to_label.pkl", "rb") as f:
#     subtype_to_label = pickle.load(f)
seqs = read_fasta(raw_sequences)
train_seqs, test_seqs, label_map, test_label_map = strat_train_test_split(seqs)

#assert that each value in the testmap aligns with the label map
for subtype, label in test_label_map.items():
    assert subtype in label_map

# # extend the label map to NUM_CLASSES
# for i in range(len(label_map), NUM_CLASSES):
#     label_map[str(i)] = i

print(label_map)

# with open(fr"{out_dir}\map_subtype_to_label.pkl", "wb") as f:
#     pickle.dump(label_map, f)

# with open(fr"{out_dir}\map_subtype_to_label_test.pkl", "wb") as f:
#     pickle.dump(test_label_map, f)




# write_csv(train_seqs, train_output, label_map)
# write_csv(test_seqs, test_output, label_map)
write_fasta(train_seqs, train_output.replace(".csv", ".fasta"), label_map)
write_fasta(test_seqs, test_output.replace(".csv", ".fasta"), label_map)
