# nolint: line_length_linter
.libPaths("C:\\Users\\Dan\\Documents\\R_Packages")
library(devtools)
devtools::install_local("E:\\projects\\seqspawnR_with_loadingbars\\SeqSpawnR")
library(SeqSpawnR)
library(seqinr)

# setname <- "retry_071124\\combined_full_and_pol"
set <- "test"

# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\", setname, "\\original_", set, "_seqs.fasta")
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\6105_pol_noalign_plus_last_additional_seqs_subtyped.fasta")
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\subtype_format_HIV-1_Subtype_REF_2020_genome_DNA_LANL_NOTaln_459taxa.fas")

# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\combined_full_and_pol\\train.fasta")
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\combined_full_and_pol\\test.fasta")
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\base_augmentation\\6105_pol_noalign\\original_6105_pol_noalign_plus_last_additional_seqs_subtyped_NORELABLENEEDED.fasta") # nolint: line_length_linter.
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\base_augmentation\\459_full_genome_noalign\\original_subtype_format_HIV-1_Subtype_REF_2020_genome_DNA_LANL_NOTaln_459taxa.fas") # nolint: line_length_linter.

# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\aug_new_combined\\base_combined\\", set, ".fasta") # nolint: line_length_linter.

# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\new_datasets\\HIV-1_FULL_GENOME_ALL_LANL_2line_clean_known.fasta")

#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\all_lanl_full_genomes\train.fasta
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\all_lanl_full_genomes\\", set, ".fasta") # nolint: line_length_linter.

#C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\new_datasets\4293_from_5225_seqs_no_unknown_no_gaps.fasta
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\new_datasets\\4293_from_5225_seqs_no_unknown_no_gaps.fasta")
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\train_4293\train.fasta
# input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\train_4293\\base\\test.fasta") 

#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol\base_aug\train.fasta
input_file <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\_compare_tg_vs_sep\\training\\pol\\base_aug\\", set, ".fasta") # nolint: line_length_linter.

# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\", setname, "\\", set, "_generated\\")
# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\base_augmentation\\6105_pol_noalign\\generated\\") # nolint: line_length_linter.
# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\base_augmentation\\459_full_genome_noalign\\generated\\") # nolint: line_length_linter.
# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\aug_new_combined\\aug_test\\", set, "_generated\\") # nolint: line_length_linter.
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation\all_lanl_full_genome

# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\base_augmentation\\all_lanl_full_genome\\", set, "_generated\\") # nolint: line_length_linter.
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\all_lanl_full_genomes
# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\all_lanl_full_genomes\\", set, "_generated\\") # nolint: line_length_linter.
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\raw\base_augmentation
# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\raw\\base_augmentation\\4293_from_5225_seqs_no_unknown_no_gaps\\", set, "_generated\\") # nolint: line_length_linter.
# gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\train_4293\\", set, "_generated\\") # nolint: line_length_linter.
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol
gen_dir <- paste0("C:\\Users\\Dan\\Desktop\\CDC\\Projects\\dives\\data\\preproc\\retry_071124\\_compare_tg_vs_sep\\training\\pol\\", set, "_generated\\") # nolint: line_length_linter.


if (!dir.exists(gen_dir)) {
  dir.create(gen_dir, recursive = TRUE)
}

if (set == "train") {
  N <- 300
} else if (set == "base") {
  N <- 5
} else {
  N <- 5
}
# S <- 150 #Number of SNPs that generated sequences can have at most
S <- 45


# Read sequences from a FASTA file
fasta_sequences <- read.fasta(input_file, whole.header=TRUE, forceDNAtolower = FALSE)

# 1. Count unique subtypes and sequences
subtype_counter <- table(sapply(strsplit(names(fasta_sequences), "\\."), `[`, 1))
print("Subtypes and their counts:")
print(subtype_counter)

# 2. Calculate spawn quotas
#spawn_quotas <- (N - subtype_counter) %/% subtype_counter
spawn_quotas <- ceiling((N - subtype_counter) / subtype_counter)
print("Spawn quotas:")
print(spawn_quotas)

# 3. & 4. Modified loop for spawn quotas and file handling
for (name in names(fasta_sequences)) {
  subtype <- unlist(strsplit(name, "\\."))[1]
  
  output_file <- paste0(gen_dir, subtype, "_generated.fasta")
  
  ref_sequence <- toupper(paste(fasta_sequences[[name]], collapse = ""))
  
  write(paste(">", name, sep = ""), output_file, append = TRUE)
  write(ref_sequence, output_file, append = TRUE)
  
  quota <- spawn_quotas[subtype]
  if (quota < 1) {
    next
  }
  spawned <- SeqSpawnR::spawn_sequences(quota, seed = ref_sequence, snps = S)
  
  for (i in 1:length(spawned)) {
    new_name <- paste(name, "_spawned_", i, sep = "")
    spawned_sequence <- toupper(paste(spawned[[i]], collapse = ""))

    #shift <- sample(1:20, 1)
    #shifted_sequence <- substr(spawned_sequence, shift, nchar(spawned_sequence))
    shifted_sequence <- spawned_sequence
    
    write(paste(">", new_name, sep = ""), output_file, append = TRUE)
    write(shifted_sequence, output_file, append = TRUE)
  }
}

# 5. Concatenate all output files
output_files <- list.files(path = gen_dir, pattern = "*_generated.fasta", full.names = TRUE)
final_output <- paste0(gen_dir, "final_output.fasta")

for (file in output_files) {
  lines <- readLines(file)
  write(lines, final_output, append = TRUE, ncolumns = 1)
}

# 6. delete intermediate files
file.remove(output_files)
