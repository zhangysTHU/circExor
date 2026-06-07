    

               

import os

import copy

import random

import itertools

import numpy as np

import pandas as pd

import warnings

import joblib

from sklearn.model_selection import train_test_split,RandomizedSearchCV

import sklearn.metrics as metrics

from sklearn.linear_model import LogisticRegression

import matplotlib.pyplot as plt

from sklearn.impute import SimpleImputer

from sklearn.utils import resample

from scipy.stats import uniform, randint

import torch

from torch.utils.data import TensorDataset, DataLoader

from sklearn.metrics import accuracy_score, roc_auc_score

import torch.optim as optim

    

def ligate_sequence(seq, add_len):

    seq = str(seq)

    seqlen = len(seq)

    if seqlen >= add_len:

        extra = seq[:add_len]

    else:

        repeat_times = add_len // seqlen

        remainder = add_len % seqlen

        extra = seq * repeat_times + seq[:remainder]

    return seq + extra

    

def balance_dataset_by_tag(df, tag_column='tag', random_state=42):

    """
    df : pd.DataFrame
        Input dataset to be processed.
    tag_column : str, default 'tag'
        Column used to split the dataset; default is 'tag'.
    random_state : int, default 42
    Returns:
    pd.DataFrame
        Balanced dataset with equal numbers of rows for tag 0 and tag 1.
    """

                               

    df_tag_0 = df[df[tag_column] == 0]

    df_tag_1 = df[df[tag_column] == 1]

                              

    if len(df_tag_0) > len(df_tag_1):

        df_tag_0_downsampled = resample(df_tag_0, 

                                        replace=False,           

                                        n_samples=len(df_tag_1),               

                                        random_state=random_state)          

        df_balanced = pd.concat([df_tag_0_downsampled, df_tag_1])

    else:

        df_tag_1_downsampled = resample(df_tag_1, 

                                        replace=False,           

                                        n_samples=len(df_tag_0),               

                                        random_state=random_state)          

        df_balanced = pd.concat([df_tag_0, df_tag_1_downsampled])

              

    df_balanced = df_balanced.sample(frac=1, random_state=random_state).reset_index(drop=True)

    return df_balanced

    

def one_hot_encode(seq, max_len=1000):

    """
    Trim or pad the input sequence to a fixed length and convert it to one-hot encoding.
    """

                

    if len(seq) > max_len:

        seq = seq[:max_len]

    else:

        seq = seq + 'N' * (max_len - len(seq))

        

                 

    NT_dict = {

        'A': [1, 0, 0, 0],

        'T': [0, 1, 0, 0],

        'C': [0, 0, 0, 1],

        'G': [0, 0, 1, 0],

        'N': [0, 0, 0, 0]

    }

    

                

    encoded_seq = [NT_dict.get(base.upper(), [0, 0, 0, 0]) for base in seq]

    

    return encoded_seq

def convert_dataframe_to_tensor(data, max_len=1500):

    import torch

    

                              

    seq_list = data["Sequence"].apply(lambda x: one_hot_encode(x, max_len=max_len)).tolist()

    

                           

    label_list = data["tag"].astype(float).tolist()

    

    seq_tensor = torch.tensor(seq_list, dtype=torch.float32)

                                                          

    seq_tensor = seq_tensor.transpose(1, 2)

    

    label_tensor = torch.tensor(label_list, dtype=torch.long)

    

    return {"seq": seq_tensor, "label": label_tensor}

    

import torch

                                

import torch.nn as nn

class circRNA_GRU(nn.Module):

    def __init__(self, seq_len=1500, num_classes=2):

        super(circRNA_GRU, self).__init__()

        

                                                 

        self.input_size = 4                                         

        self.hidden_size = 128            

        self.num_layers = 3               

        self.bidirectional = True      

        

                             

                                                               

        self.gru = nn.GRU(

            input_size=self.input_size,

            hidden_size=self.hidden_size,

            num_layers=self.num_layers,

            batch_first=True,

            bidirectional=self.bidirectional

        )

        

                        

                                                                        

        num_directions = 2 if self.bidirectional else 1

        mlp_input_dim = self.num_layers * num_directions * self.hidden_size                    

        head_hidden_dim = 256

        dropout_rate = 0.1

        

                                

        self.classifier = nn.Sequential(

                      

            nn.Linear(mlp_input_dim, head_hidden_dim),

            nn.BatchNorm1d(head_hidden_dim),

            nn.ReLU(),

            nn.Dropout(dropout_rate),

            

                                                

            nn.Linear(head_hidden_dim, head_hidden_dim),

            nn.BatchNorm1d(head_hidden_dim),

            nn.ReLU(),

            nn.Dropout(dropout_rate),

            

                                               

            nn.Linear(head_hidden_dim, num_classes)

        )

    def forward(self, x):

                                                     

                                                   

        x = x.permute(0, 2, 1) 

        

                

                                                                

                                                                 

        _, h = self.gru(x)

        

                                  

                                                                            

        batch_size = h.size(1)

                               

        h_flatten = h.transpose(0, 1).contiguous().view(batch_size, -1)

        

               

        logits = self.classifier(h_flatten)

        return logits

    

                              

def evaluate_performance(y_test, y_pred, y_prob):

           

    auroc = metrics.roc_auc_score(y_test,y_prob)

    auroc_curve = metrics.roc_curve(y_test, y_prob)

           

    auprc=metrics.average_precision_score(y_test, y_prob) 

    auprc_curve=metrics.precision_recall_curve(y_test, y_prob)

             

    accuracy=metrics.accuracy_score(y_test,y_pred) 

        

    mcc=metrics.matthews_corrcoef(y_test,y_pred)

    

    recall=metrics.recall_score(y_test, y_pred)

    precision=metrics.precision_score(y_test, y_pred)

    f1=metrics.f1_score(y_test, y_pred)

    class_report=metrics.classification_report(y_test, y_pred,target_names = ["control","case"])

    model_perf = {"auroc":auroc,"auroc_curve":auroc_curve,

                  "auprc":auprc,"auprc_curve":auprc_curve,

                  "accuracy":accuracy, "mcc": mcc,

                  "recall":recall,"precision":precision,"f1":f1,

                  "class_report":class_report}

        

    return model_perf

    

                             

def eval_output(model_perf,path):

    with open(os.path.join(path,"Evaluate_Result_TestSet.txt"),'w') as f:

        f.write("AUROC=%s\tAUPRC=%s\tAccuracy=%s\tMCC=%s\tRecall=%s\tPrecision=%s\tf1_score=%s\n" %

               (model_perf["auroc"],model_perf["auprc"],model_perf["accuracy"],model_perf["mcc"],model_perf["recall"],model_perf["precision"],model_perf["f1"]))

        f.write("\n######NOTE#######\n")

        f.write("#According to help_documentation of sklearn.metrics.classification_report:in binary classification, recall of the positive class is also known as sensitivity; recall of the negative class is specificity#\n\n")

        f.write(model_perf["class_report"])

    

def train_and_evaluate_model(model, train_loader, val_loader, epochs, device, output_dir):

    criterion = nn.CrossEntropyLoss()

    optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)            

    

                

    for epoch in range(epochs):

        model.train()

        train_loss = 0.0

        

                      

        for batch_x, batch_y in train_loader:

            batch_x, batch_y = batch_x.to(device), batch_y.to(device)

            

            optimizer.zero_grad()

            outputs = model(batch_x)

            loss = criterion(outputs, batch_y)

            loss.backward()

            optimizer.step()

            

            train_loss += loss.item() * batch_x.size(0)

        

        train_loss /= len(train_loader.dataset)

        

                      

        model.eval()

        val_loss = 0.0

        all_preds = []

        all_labels = []

        all_probs = []

        

        with torch.no_grad():

            for batch_x, batch_y in val_loader:

                batch_x, batch_y = batch_x.to(device), batch_y.to(device)

                outputs = model(batch_x)

                loss = criterion(outputs, batch_y)

                val_loss += loss.item() * batch_x.size(0)

                

                probs = torch.softmax(outputs, dim=1)[:, 1]                

                _, preds = torch.max(outputs, 1)          

                

                all_probs.extend(probs.cpu().numpy())

                all_preds.extend(preds.cpu().numpy())

                all_labels.extend(batch_y.cpu().numpy())

                

        val_loss /= len(val_loader.dataset)

        val_acc = accuracy_score(all_labels, all_preds)

        

        try:

            val_auc = roc_auc_score(all_labels, all_probs)

        except ValueError:

            val_auc = 0.0                       

            

        print(f"Epoch [{epoch+1}/{epochs}] - Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f} | Val Acc: {val_acc:.4f} | Val AUC: {val_auc:.4f}")

                                   

    print("\nTraining completed. Evaluating on validation set...")

    model_perf = evaluate_performance(all_labels, all_preds, all_probs)

    

               

    if not os.path.exists(output_dir):

        os.makedirs(output_dir)

    eval_output(model_perf, output_dir)

    print(f"Evaluation results saved to {output_dir}/Evaluate_Result.txt")

    

                 

    torch.save(model.state_dict(), os.path.join(output_dir, "circRNA_GRU_1500nt.pth"))

    

    return model_perf

    

             

SEED = 100

random.seed(SEED)

np.random.seed(SEED)

warnings.filterwarnings(action='ignore')

            

output_dir = "circExor/DL_models/circRNA_DL_Model_tridivided_intra5fold_Output/GRU"

if not (os.path.exists(output_dir)):

    os.mkdir(output_dir)

    

           

dataset = pd.read_csv(

    'circExor/reference_preprocessing/circRNA/output_with_sequences.csv',

    sep='\t',

    index_col=False

)

dataset_filtered = dataset

                           

add_len = 4 

dataset_filtered["Sequence"] = dataset_filtered["Sequence"].apply(

    lambda x: ligate_sequence(x, add_len)

)

           

                                                                              

                   

                   

                                

    

           

dataset_filtered = balance_dataset_by_tag(dataset_filtered, tag_column='tag', random_state=42)

                         

train_df, temp_df = train_test_split(dataset_filtered, test_size=0.4, random_state=SEED, stratify=dataset_filtered['tag'])

val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=SEED, stratify=temp_df['tag'])

                                

train_tensors = convert_dataframe_to_tensor(train_df, max_len=1500)

val_tensors = convert_dataframe_to_tensor(val_df, max_len=1500)

test_tensors = convert_dataframe_to_tensor(test_df, max_len=1500)

      

torch.save(train_tensors, os.path.join(output_dir, "train_tensors_1500nt.pt"))

torch.save(val_tensors, os.path.join(output_dir, "val_tensors_1500nt.pt"))

torch.save(test_tensors, os.path.join(output_dir, "test_tensors_1500nt.pt"))

    

                    

x_train, y_train = train_tensors["seq"], train_tensors["label"]

x_val, y_val = val_tensors["seq"], val_tensors["label"]

                              

print(f"Train shapes - X: {x_train.shape}, y: {y_train.shape}")

print(f"Val shapes   - X: {x_val.shape}, y: {y_val.shape}")

                                      

val_dataset = TensorDataset(x_val, y_val)

val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False)

                                          

from sklearn.model_selection import KFold

import pandas as pd

import numpy as np

      

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print(f"Using device: {device}")

kf = KFold(n_splits=5, shuffle=True, random_state=42)

                            

cv_metrics_list = [] 

fold_aucs = []

print("\n--- Starting 5-Fold Cross Validation on Training Set ---")

for fold, (train_idx, val_idx) in enumerate(kf.split(x_train)):

    print(f"\n--- Fold {fold + 1} ---")

    

                         

    x_train_fold = x_train[train_idx]

    y_train_fold = y_train[train_idx]

    x_val_fold = x_train[val_idx]

    y_val_fold = y_train[val_idx]

    

                            

    train_dataset_fold = TensorDataset(x_train_fold, y_train_fold)

    val_dataset_fold = TensorDataset(x_val_fold, y_val_fold)

    

    train_loader_fold = DataLoader(train_dataset_fold, batch_size=64, shuffle=True)

    val_loader_fold = DataLoader(val_dataset_fold, batch_size=64, shuffle=False)

    

                    

    model_fold = circRNA_GRU(seq_len=1500, num_classes=2).to(device)

    

                   

    fold_output_dir = os.path.join(output_dir, f"GRU_Fold_{fold+1}")

    if not os.path.exists(fold_output_dir):

        os.makedirs(fold_output_dir, exist_ok=True)

        

                                            

    fold_metrics = train_and_evaluate_model(

        model=model_fold, 

        train_loader=train_loader_fold, 

        val_loader=val_loader_fold, 

        epochs=15, 

        device=device, 

        output_dir=fold_output_dir 

    )

    

    fold_aucs.append(fold_metrics['auroc'])

    

                                      

    cv_metrics_list.append({

        "Fold": f"Fold_{fold+1}",

        "AUROC": fold_metrics.get('auroc', 0),

        "AUPRC": fold_metrics.get('auprc', 0),

        "Accuracy": fold_metrics.get('accuracy', 0),

        "Precision": fold_metrics.get('precision', 0),

        "Recall": fold_metrics.get('recall', 0),

        "F1_Score": fold_metrics.get('f1', 0),

        "MCC": fold_metrics.get('mcc', 0)

    })

                                       

cv_metrics_df = pd.DataFrame(cv_metrics_list)

cv_metrics_df.to_csv(os.path.join(output_dir, "circRNA_GRU_5fold_cv_metrics.csv"), index=False)

print(f"\n5-Fold CV metrics saved to: {os.path.join(output_dir, 'circRNA_GRU_5fold_cv_metrics.csv')}")

print("\n--- 5-Fold Cross Validation Results ---")

print(f"Average AUROC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}")

print("If you are satisfied with the hyper-parameters, we will now train on FULL x_train data...")

                

print("\n--- Starting Final Training on FULL Training Set ---")

train_dataset_full = TensorDataset(x_train, y_train)

train_loader_full = DataLoader(train_dataset_full, batch_size=64, shuffle=True)

final_model = circRNA_GRU(seq_len=1500, num_classes=2).to(device)

final_output_dir = os.path.join(output_dir, "GRU_Final_Model")

if not os.path.exists(final_output_dir):

    os.makedirs(final_output_dir, exist_ok=True)

perf_metrics = train_and_evaluate_model(

    model=final_model, 

    train_loader=train_loader_full, 

    val_loader=val_loader,                   

    epochs=20, 

    device=device, 

    output_dir=final_output_dir

)

        

print("\nValidation Set Final Performance (Trained on Full Train Set):")

print(f"AUROC: {perf_metrics['auroc']:.4f}")

print(f"Accuracy: {perf_metrics['accuracy']:.4f}")

print(f"F1 Score: {perf_metrics['f1']:.4f}")

