cd finetune

export TOKENIZERS_PARALLELISM=false

# pol training
# export DATA_PATH=/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/unshifted_data
# export MAX_LENGTH=1000 # Please set the number as 0.25 * your sequence length. 
											# e.g., set it as 250 if your DNA sequences have 1000 nucleotide bases
											# This is because the tokenized will reduce the sequence length by about 5 times

# combined full and pol training
#data/preproc/retry_071124/combined_full_and_pol/final
# export DATA_PATH=/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/combined_full_and_pol/final

# export DATA_PATH=/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/aug_new_combined/data

# export DATA_PATH=/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/all_lanl_full_genomes/training
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\train_4293
# export DATA_PATH=/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/train_4293
#C:\Users\Dan\Desktop\CDC\Projects\dives\data\preproc\retry_071124\_compare_tg_vs_sep\training\pol
export DATA_PATH=/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/_compare_tg_vs_sep/training/combined

export MAX_LENGTH=2500
export LR=0.000005

# Training use DataParallel (added eval_accumulation_steps)
# --model_name_or_path zhihan1996/DNABERT-2-117M \
# --model_name_or_path /mnt/c/users/dan/desktop/cdc/projects/dives/retry_07_11_24/DNABERT_2/finetune/output/dnabert2-pol_only_18900/checkpoint-18900 \
python train.py \
    --model_name_or_path zhihan1996/DNABERT-2-117M \
    --data_path  ${DATA_PATH} \
    --kmer -1 \
    --run_name DNABERT2_${DATA_PATH} \
    --model_max_length ${MAX_LENGTH} \
    --per_device_train_batch_size 1 \
    --per_device_eval_batch_size 1 \
    --gradient_accumulation_steps 1 \
    --learning_rate ${LR} \
    --num_train_epochs 2 \
    --fp16 \
    --save_steps 5000 \
    --output_dir output/dnabert2 \
    --evaluation_strategy steps \
    --eval_steps 5000 \
    --warmup_steps 50 \
    --logging_steps 100 \
    --overwrite_output_dir True \
    --log_level info \
    --find_unused_parameters False \
    --dataloader_num_workers 0 \
    --eval_accumulation_steps 30 \
    
# Training use DistributedDataParallel (more efficient)
# export num_gpu=1 # please change the value based on your setup

# torchrun --nproc-per-node=${num_gpu} train.py \
#     --model_name_or_path zhihan1996/DNABERT-2-117M \
#     --data_path  ${DATA_PATH} \
#     --kmer -1 \
#     --run_name DNABERT2_${DATA_PATH} \
#     --model_max_length ${MAX_LENGTH} \
#     --per_device_train_batch_size 8 \
#     --per_device_eval_batch_size 16 \
#     --gradient_accumulation_steps 1 \
#     --learning_rate ${LR} \
#     --num_train_epochs 5 \
#     --fp16 \
#     --save_steps 200 \
#     --output_dir output/dnabert2 \
#     --evaluation_strategy steps \
#     --eval_steps 200 \
#     --warmup_steps 50 \
#     --logging_steps 100 \
#     --overwrite_output_dir True \
#     --log_level info \
#     --find_unused_parameters False

