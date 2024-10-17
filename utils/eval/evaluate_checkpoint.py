import pickle
import csv
import pandas as pd
import numpy as np
from sklearn.metrics import classification_report
from transformers import AutoTokenizer, AutoModelForSequenceClassification, Trainer
import transformers
from torch.utils.data import Dataset
from dataclasses import dataclass, field
from typing import Optional
import torch
import os

os.environ["TOKENIZERS_PARALLELISM"] = "false"



# -- name of the folder containing test.csv

# TEST_DIR = "NXT_and_SA"
# TEST_DIR = "combined_full_and_pol/final"
# TEST_DIR = "validation_datasets"
TEST_DIR = "validation_datasets/China"
# TEST_DIR = "validation_datasets/Bulgaria"

# TEST_FILE_NAME = "test_file_from_training.csv"
# TEST_FILE_NAME = "new_curated_5225_seqs.csv"
# TEST_FILE_NAME = "5225_seqs_no_gaps.csv"
# TEST_FILE_NAME = "4293_from_5225_seqs_no_unknown.csv"

# TEST_FILE_NAME = "LANL_minus_4293_curated.csv"
# TEST_FILE_NAME = "NXTSA_properly_labeled.csv"
TEST_FILE_NAME = "ChinaValidation.csv"
# TEST_FILE_NAME = "BulgariaValidation.csv"



# -- name of the folder in which to store the results
# MODEL_ID = "500aug-full-and-pol"
# MODEL_ID = "18900_pol_only"
# MODEL_ID = "4293_model_curated_full"

MODEL_ID = "new_full_model"
# MODEL_ID = "new_pol_only_model"
# MODEL_ID = "new_combined_model"


# TEST_NAME = TEST_DIR+"_"+MODEL_ID
# TEST_NAME = "NXTSA_pol_" + MODEL_ID
# TEST_NAME = "LANL_minus_4293_curated_" + MODEL_ID
# TEST_NAME = "NXTSA_" 
# TEST_NAME = "LANL_minus_4293_curated_" + MODEL_ID
# TEST_NAME = "new_NXTSA_validation_"
TEST_NAME = "ChinaValidation_" + MODEL_ID
# TEST_NAME = "BulgariaValidation_"
TEST_NAME += MODEL_ID
if not os.path.exists(TEST_NAME):
    os.makedirs(TEST_NAME)


# CHECKPOINT_PATH = f"/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/retry_07_11_24/DNABERT_2/finetune/output/dnabert2-500aug-full-and-pol-86/checkpoint-30000"
#C:\Users\Dan\Desktop\CDC\Projects\dives\retry_07_11_24\DNABERT_2\finetune\output\dnabert2-all_lanl_full\checkpoint-150000
# CHECKPOINT_PATH = f"/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/retry_07_11_24/DNABERT_2/finetune/output/dnabert2-all_lanl_full/checkpoint-150000"
# C:\Users\Dan\Desktop\CDC\Projects\dives\retry_07_11_24\DNABERT_2\finetune\output\dnabert2-pol_only_18900\checkpoint-18900
# CHECKPOINT_PATH = f"/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/retry_07_11_24/DNABERT_2/finetune/output/dnabert2-pol_only_18900/checkpoint-18900"

# NEW FULL MODEL
# CHECKPOINT_PATH = f"/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/retry_07_11_24/DNABERT_2/finetune/output/dnabert2_curated_full_4293/checkpoint-45000"


# NEW POL ONLY MODEL
#C:\Users\Dan\Desktop\CDC\Projects\dives\retry_07_11_24\DNABERT_2\finetune\output\__dnabert2pol\checkpoint-30000
# CHECKPOINT_PATH = f"/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/retry_07_11_24/DNABERT_2/finetune/output/__dnabert2pol/checkpoint-30000"


# NEW COMBINED MODEL
#C:\Users\Dan\Desktop\CDC\Projects\dives\retry_07_11_24\DNABERT_2\finetune\output\_dnabert2_combined\checkpoint-95000
CHECKPOINT_PATH = f"/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/retry_07_11_24/DNABERT_2/finetune/output/_dnabert2_combined/checkpoint-95000"

TEST_LABELS_PATH = "/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/unified_map_subtype_to_label.pkl"



@dataclass
class ModelArguments:
    model_name_or_path: Optional[str] = field(default="zhihan1996/DNABERT-2-117M")

@dataclass
class DataArguments:
    test_data_path: str = field(default=f'/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/{TEST_DIR}/{TEST_FILE_NAME}', metadata={"help": "Path to the test data."})
    kmer: int = field(default=-1, metadata={"help": "k-mer for input sequence. -1 means not using k-mer."})

@dataclass
class TrainingArguments(transformers.TrainingArguments):
    cache_dir: Optional[str] = field(default=None)
    run_name: str = field(default="DNABERT2")
    model_max_length: int = field(default=2500, metadata={"help": "Maximum sequence length."})
    output_dir: str = field(default="output")
    per_device_eval_batch_size: int = field(default=8)
    eval_accumulation_steps: int = field(default=1)

class SupervisedDataset(Dataset):
    def __init__(self, data_path: str, tokenizer: transformers.PreTrainedTokenizer, kmer: int = -1):
        super(SupervisedDataset, self).__init__()
        with open(data_path, "r") as f:
            data = list(csv.reader(f))

        self.texts = [d[0] for d in data]
        labels = [int(d[1]) for d in data]
        names = [d[2] if len(d) > 2 else "" for d in data]
        
        output = tokenizer(
            self.texts,
            return_tensors="pt",
            padding="longest",
            max_length=tokenizer.model_max_length,
            truncation=True,
        )
        self.input_ids = output["input_ids"]
        self.attention_mask = output["attention_mask"]
        self.labels = labels
        self.num_labels = len(set(labels))
        self.names = names

        self.current_batch_names = []
        self.current_batch_lengths = []

    def __len__(self):
        return len(self.input_ids)

    def __getitem__(self, i) -> dict:
        self.current_batch_names.append(self.names[i])
        self.current_batch_lengths.append(len(self.texts[i]))
        return dict(input_ids=self.input_ids[i], labels=self.labels[i], attention_mask=self.attention_mask[i])

# Closure to accumulate metrics batch-wise and process logits
def preprocess_with_data(label_dict, seq_names, label_map, dataset):
    accumulated_predictions = []
    accumulated_labels = []
    misclassifications = {}
    misclass_sequences = []
    confidences = {'id': [], 'true_class': [], 'predicted_class': [], 'true_class_confidence': [], 
                   'predicted_class_confidence': [], 'top_3_confidences': [], 'position_of_true_class': []}
    misclasses = set()

    batch_names_list = dataset.current_batch_names
    batch_lengths_list = dataset.current_batch_lengths

    def preprocess_logits_for_metrics(logits, labels):
        logits_np = logits[0].detach().cpu().numpy()
        labels_np = labels.detach().cpu().numpy()

        # assert len(batch_names_list) == 8
        

        # Compute predictions (argmax) for this batch
        predictions = np.argmax(logits_np, axis=1)

        # Accumulate predictions and labels for final metrics computation
        accumulated_predictions.extend(predictions)
        accumulated_labels.extend(labels_np)

        # Process each sample in the batch
        for i, logit in enumerate(logits_np):
            softmax_logit = torch.nn.functional.softmax(torch.tensor(logit), dim=0)

            # Track confidences
            confidences['id'].append(batch_names_list[i])
            confidences['true_class'].append(label_dict[labels_np[i]])
            confidences['predicted_class'].append(label_dict[predictions[i]])
            confidences['true_class_confidence'].append(softmax_logit[labels_np[i]].item())
            confidences['predicted_class_confidence'].append(softmax_logit[predictions[i]].item())
            confidences['top_3_confidences'].append([(label_dict[j], round(softmax_logit[j].item(), 4)) 
                                                      for j in np.argsort(logit)[-3:][::-1]])
            confidences['position_of_true_class'].append(np.where(np.argsort(logit)[::-1] == labels_np[i])[0][0] + 1)
            confidences['length'] = batch_lengths_list[i]

            # If the confidence of the true class is less than 0.6, add to misclasses
            if softmax_logit[labels_np[i]].item() < 0.6:
                misclasses.add(label_dict[labels_np[i]])

            # Collect misclassification information
            if labels_np[i] != predictions[i]:
                if labels_np[i] in misclassifications:
                    if predictions[i] in misclassifications[labels_np[i]]:
                        misclassifications[labels_np[i]][predictions[i]] += 1
                    else:
                        misclassifications[labels_np[i]][predictions[i]] = 1
                else:
                    misclassifications[labels_np[i]] = {predictions[i]: 1}
                misclass_sequences.append(seq_names[i])

        batch_names_list.clear()
        batch_lengths_list.clear()
        return torch.tensor(predictions)

    # Function to compute final metrics after accumulation
    def compute_final_metrics():
        # Generate dynamic label list based on predictions and true labels
        unique_true_classes = set(accumulated_labels)
        unique_pred_classes = set(accumulated_predictions)
        label_list = []
        for p in sorted(unique_true_classes.union(unique_pred_classes)):
            for key, value in label_map.items():
                if value == p:
                    label_list.append(key)
                    break

        # Compute classification report using the dynamically generated label_list
        report = classification_report(accumulated_labels, accumulated_predictions, target_names=label_list, output_dict=True, zero_division=0)

        return report, confidences, misclassifications, misclasses, misclass_sequences

    return preprocess_logits_for_metrics, compute_final_metrics

def evaluate_model():
    checkpoint_path = CHECKPOINT_PATH
    model_args = ModelArguments()
    data_args = DataArguments()
    training_args = TrainingArguments()

    # Load model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=training_args.cache_dir,
        model_max_length=training_args.model_max_length,
        padding_side="right",
        use_fast=True,
        trust_remote_code=True,
    )
    print("loading model")
    model = AutoModelForSequenceClassification.from_pretrained(
        checkpoint_path,
        trust_remote_code=True)
    print("model loaded")

    # Load label map
    with open(TEST_LABELS_PATH, 'rb') as f:
        label_map = pickle.load(f)

    print("label map loaded", label_map)

    # Load test dataset
    test_dataset = SupervisedDataset(data_path=data_args.test_data_path, tokenizer=tokenizer, kmer=data_args.kmer)
    print("test dataset loaded")

    # Prepare for dynamic label list generation and batch-wise processing
    label_dict = {value: key for key, value in label_map.items()}
    preprocess_fn, compute_final_metrics_fn = preprocess_with_data(label_dict, test_dataset.names, label_map, test_dataset)

    # Initialize trainer
    trainer = Trainer(
        model=model, 
        tokenizer=tokenizer,
        args=training_args,
        preprocess_logits_for_metrics=preprocess_fn
    )
    print("trainer initialized")

    # Run predictions
    print("starting predictions")
    outputs = trainer.predict(test_dataset)
    print("finished predictions")

    # Compute final metrics based on accumulated predictions and labels
    final_report, confidences, misclassifications, misclasses, misclass_sequences = compute_final_metrics_fn()

    # 1. Evaluation metrics per class
    metrics_df = pd.DataFrame(final_report).transpose()
    metrics_df.to_csv(f'/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/{TEST_NAME}/{TEST_NAME}_evaluation_metrics_per_class.csv', index=True)

    # 2. Confidences
    confidences_df = pd.DataFrame(confidences)
    confidences_df["match?"] = confidences_df["true_class"] == confidences_df["predicted_class"]
    confidences_df.to_csv(f'/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/{TEST_NAME}/{TEST_NAME}_confidences.csv', index=False)

    # 3. Misclassifications
    with open(f'/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/{TEST_NAME}/{TEST_NAME}_misclassifications.csv', 'w') as f:
        f.write('true_class,misclassified_as, misclassified_seqs\n')
        for true_class, misclassified_as in misclassifications.items():
            t_c = [key for key, value in label_map.items() if value == true_class][0]
            f.write(f"{t_c},[")
            for misclass, count in misclassified_as.items():
                m_c = [key for key, value in label_map.items() if value == misclass][0]
                f.write(f"{m_c}:{count} ")
            f.write('],')
            f.write(';'.join(misclass_sequences))
            f.write('\n')

    # 4. Misclassified classes
    with open(f'/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/{TEST_NAME}/{TEST_NAME}_misclassified_classes.txt', 'w') as f:
        for misclass in misclasses:
            f.write(f"{misclass}\n")

    print("Results saved to", f'/mnt/c/Users/Dan/Desktop/CDC/Projects/dives/data/preproc/retry_071124/{TEST_NAME}/{TEST_NAME}_evaluation_metrics_per_class.csv')

if __name__ == "__main__":
    evaluate_model()