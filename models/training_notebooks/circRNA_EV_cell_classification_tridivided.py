# %%
# Python import
import os
import copy
import random
import itertools
import numpy as np
import pandas as pd
import lightgbm as lgb
import warnings
import joblib
from sklearn.model_selection import train_test_split,RandomizedSearchCV
import sklearn.metrics as metrics
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.utils import resample
import lightgbm as lgb
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
from scipy.stats import uniform, randint
from ngboost import NGBClassifier
from catboost import CatBoostClassifier
from ngboost.distns import k_categorical
from xgboost import XGBClassifier

# %%
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

# %%
def balance_dataset_by_tag(df, tag_column='tag', random_state=42):
    """
    df : pd.DataFrame
        需要处理的原始数据集。
    tag_column : str, 默认 'tag'
        用于分割数据的列名，默认是 'tag'。
    random_state : int, 默认 42
    返回:
    pd.DataFrame
        平衡后的数据集，tag列为0和1的行数一致。
    """
    # 分割数据集为 tag=0 和 tag=1 的两部分
    df_tag_0 = df[df[tag_column] == 0]
    df_tag_1 = df[df[tag_column] == 1]

    # 判断哪个数据集较大，然后对较大的数据集进行下采样
    if len(df_tag_0) > len(df_tag_1):
        df_tag_0_downsampled = resample(df_tag_0, 
                                        replace=False,  # 不允许重复采样
                                        n_samples=len(df_tag_1),  # 调整为与tag=1一致
                                        random_state=random_state)  # 设置随机种子
        df_balanced = pd.concat([df_tag_0_downsampled, df_tag_1])
    else:
        df_tag_1_downsampled = resample(df_tag_1, 
                                        replace=False,  # 不允许重复采样
                                        n_samples=len(df_tag_0),  # 调整为与tag=0一致
                                        random_state=random_state)  # 设置随机种子
        df_balanced = pd.concat([df_tag_0, df_tag_1_downsampled])

    # 打乱数据集的顺序
    df_balanced = df_balanced.sample(frac=1, random_state=random_state).reset_index(drop=True)

    return df_balanced


# %%
# Count the frequency of k-mer in each RNA sequence
# k-mer was normalized by total k-mer count of each RNA sequence
def _count_kmer(Dataset, k):  # k = 3, 4, 5
    
    # copy dataset
    dataset = copy.deepcopy(Dataset)
    # alphabet of nucleotide
    nucleotide = ['A', 'C', 'G', 'T']
    
    # generate k-mers
    #  k == 5:
    five = list(itertools.product(nucleotide, repeat=5))
    pentamer = [''.join(n) for n in five]
    
    #  k == 4:
    four = list(itertools.product(nucleotide, repeat=4))
    tetramer = [''.join(n) for n in four]

    # k == 3:
    three = list(itertools.product(nucleotide, repeat=3))
    threemer = [''.join(n) for n in three]
    
    # input features can be combinations of different k values
    if k == 34:
        table_kmer = dict.fromkeys(threemer, 0)
        table_kmer.update(dict.fromkeys(tetramer, 0))
    elif k == 45:
        table_kmer = dict.fromkeys(tetramer, 0)
        table_kmer.update(dict.fromkeys(pentamer, 0))
    elif k == 345:
        table_kmer = dict.fromkeys(threemer, 0)
        table_kmer.update(dict.fromkeys(tetramer, 0))
        table_kmer.update(dict.fromkeys(pentamer, 0))

    # count k-mer for each sequence
    for mer in table_kmer.keys():
        table_kmer[mer] = dataset["Sequence"].apply(lambda x: x.count(mer))
    
    # for k-mer raw count without normalization, index: nuc:1 or cyto:0
    rawcount_kmer_df = pd.DataFrame(table_kmer)
    df1_rawcount = pd.concat([rawcount_kmer_df, dataset["RNA_Symbol"]], axis=1)
    df1_rawcount.index = dataset["tag"]

    # for k-mer frequency with normalization, index: nuc:1 or cyto:0
    freq_kmer_df = rawcount_kmer_df.apply(lambda x: x / x.sum(), axis=1)
    df1 = pd.concat([freq_kmer_df, dataset["RNA_Symbol"]], axis=1)
    df1.index = dataset["tag"]

    return df1, df1_rawcount


# %%
#Evaluate performance of model
def evaluate_performance(y_test, y_pred, y_prob):
    # AUROC
    auroc = metrics.roc_auc_score(y_test,y_prob)
    auroc_curve = metrics.roc_curve(y_test, y_prob)
    # AUPRC
    auprc=metrics.average_precision_score(y_test, y_prob) 
    auprc_curve=metrics.precision_recall_curve(y_test, y_prob)
    #Accuracy
    accuracy=metrics.accuracy_score(y_test,y_pred) 
    #MCC
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

# %%
# Output result of evaluation
def eval_output(model_perf,path):
    with open(os.path.join(path,"Evaluate_Result_TestSet.txt"),'w') as f:
        f.write("AUROC=%s\tAUPRC=%s\tAccuracy=%s\tMCC=%s\tRecall=%s\tPrecision=%s\tf1_score=%s\n" %
               (model_perf["auroc"],model_perf["auprc"],model_perf["accuracy"],model_perf["mcc"],model_perf["recall"],model_perf["precision"],model_perf["f1"]))
        f.write("\n######NOTE#######\n")
        f.write("#According to help_documentation of sklearn.metrics.classification_report:in binary classification, recall of the positive class is also known as sensitivity; recall of the negative class is specificity#\n\n")
        f.write(model_perf["class_report"])

# %%
# Plot AUROC of model
def plot_AUROC(model_perf,path):
    #get AUROC,FPR,TPR and threshold
    roc_auc = model_perf["auroc"]
    fpr,tpr,threshold = model_perf["auroc_curve"]
    #return AUROC info
    temp_df = pd.DataFrame({"FPR":fpr,"TPR":tpr})
    temp_df.to_csv(os.path.join(path,"AUROC_info.txt"),header = True,index = False, sep = '\t')
    #plot
    plt.figure()
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='AUROC (area = %0.2f)' % roc_auc) 
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("AUROC of Models")
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(path,"AUROC_TestSet.pdf"),format = "pdf")

# %%
def logistic_regression_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42):
    """
    执行Logistic Regression二分类任务，包含超参数调优、训练、评估和结果保存。

    参数：
    x_train : 训练集特征
    y_train : 训练集标签
    x_test : 测试集特征
    y_test : 测试集标签
    output_dir : 输出目录路径
    SEED : 随机种子，默认为42
    """
    print("\n*** Logistic Regression  ***")

    # Logistic Regression params
    lr_param_dict = {
        "penalty": ["l2"],
        "C": [1e-3, 5e-3, 1e-2, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000],
        "solver": ["liblinear"],
        "random_state": [SEED]
    }

    # Initiate model
    lr_model = LogisticRegression()

    # Adjust hyper-parameters with 5-fold cross validation
    lr_rscv = RandomizedSearchCV(lr_model, lr_param_dict, n_iter=100, cv=5, verbose=0,
                                 scoring="roc_auc", random_state=SEED, n_jobs=30)
    lr_rscv.fit(x_train, y_train)

    # Output path
    path = os.path.join(output_dir, "LogisticRegression")
    if not os.path.exists(path):
        os.mkdir(path)

    # Model performance (AUROC) on cross-validation dataset
    lr_cv_perf = np.array([lr_rscv.cv_results_["split%s_test_score" % str(i)] for i in range(5)])[:, lr_rscv.best_index_]

    # Get best model with score [max(mean(auc(5 cross-validation)))]
    lr_best_model = lr_rscv.best_estimator_

    # Get predict_class(y_pred) and predict_probality_for_case(y_prob) of TestSet
    y_pred = lr_best_model.predict(x_test)
    y_prob = lr_best_model.predict_proba(x_test)[:, 1]

    # Get model performance
    model_perf = evaluate_performance(y_test, y_pred, y_prob)

    # Output result of evaluation
    eval_output(model_perf, path)

    # You can make bar plot consisted of accuracy, sensitivity, specificity, auroc, f1 score, MCC, precision, recall, auprc according to the "Evaluate_Result_TestSet.txt"
    # Plot AUROC
    plot_AUROC(model_perf, path)

    # Save model
    joblib.dump(lr_best_model, os.path.join(path, "best_LogisticRegression_model.pkl"))

    return model_perf

# %%
def SVM_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42):
    #Construct SVM model
    print("\n*** SVM ***")

    # SVM params
    SVM_param_dict = {
        'kernel':('linear', 'rbf'), 
        'C':[0.01,0.1,1,10, 100], 
        'gamma':[0.001, 0.005, 0.1 ,0.5,1, 2],
        "probability":[True],
        "random_state":[SEED]
    }

    #Initiate model
    SVM_model = svm.SVC()
    #Adjust hyper-parameters with 5-fold cross validation
    SVM_rscv = RandomizedSearchCV(SVM_model, SVM_param_dict, n_iter=100,cv = 5,verbose = 0,
                            scoring = "roc_auc",random_state=SEED,n_jobs =30)
    SVM_rscv.fit(x_train, y_train) 


    #Evaluate best SVM model
    #Output path
    path = os.path.join(output_dir,"SVM")
    if not (os.path.exists(path)):
        os.mkdir(path)
        
    # Model performance(AUROC) on cross-validation dataset
    SVM_cv_perf = np.array([ SVM_rscv.cv_results_["split%s_test_score"%str(i)] for i in range(5)])[:,SVM_rscv.best_index_]

    #Get best model with score [max(mean(auc(5 cross validation)))]
    svm_best_model = SVM_rscv.best_estimator_
    #Get predict_class(y_pred) and predict_probality_for_case(y_prob) of TestSet
    y_pred = svm_best_model.predict(x_test)
    y_prob = svm_best_model.predict_proba(x_test)[:,1]

    #Get model performance
    model_perf = evaluate_performance(y_test,y_pred,y_prob)
    #Output result of evaluation
    eval_output(model_perf,path)
    #You can make bar plot consisted of accuracy,sensitivity,specificity,auroc,f1 score,MCC,precision,recall,auprc according to the "Evaluate_Result_TestSet.txt"
    # Plot AUROC
    plot_AUROC(model_perf,path)

    #save model
    joblib.dump(svm_best_model,os.path.join(path,"best_SVM_model.pkl"))
    #load model
    #svm_best_model = joblib.load(os.path.join(path,"best_SVM_model.pkl"))
    ######################################
    return model_perf



# %%
def lightGBM_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_iter = 100):
    print("\n*** LightGBM ***")

    # LightGBM params
    lgb_param_dict = {
        "learning_rate":[0.1, 0.05, 0.02, 0.01],
        "num_leaves": range(10,36,5),
        "max_depth" : [2,3,4,5,10,20,40,50],
        "min_child_samples": range(1, 45, 2),
        "colsample_bytree" : [i / 10 for i in range(2,11)],
        "metric" : ["binary_logloss"],
        "n_jobs":[1],
        "n_estimators" : range(100,2500,100),
        "subsample" :  [i / 10 for i in range(2, 11)],
        "subsample_freq" : [0, 1, 2],
        "reg_alpha" : [0, 0.001, 0.005, 0.01, 0.1],
        "reg_lambda" : [0, 0.001, 0.005, 0.01, 0.1],
        "objective":["binary"],
        "random_state":[SEED]
    }

    #Initiate model
    lgb_model = lgb.LGBMClassifier()
    #Adjust hyper-parameters with 5-fold cross validation
    lgb_rscv = RandomizedSearchCV(lgb_model, lgb_param_dict, n_iter=n_iter,cv = 5,verbose = 0,
                            scoring = "roc_auc",random_state=SEED,n_jobs = 30)
    lgb_rscv.fit(x_train, y_train)   


    #Evaluate best LightGBM model
    #Output path
    path = os.path.join(output_dir,"LightGBM")
    if not (os.path.exists(path)):
        os.mkdir(path)
        
    # Model performance(AUROC) on cross-validation dataset
    lgb_cv_perf = np.array([ lgb_rscv.cv_results_["split%s_test_score"%str(i)] for i in range(5)])[:,lgb_rscv.best_index_]

    #Get best model with score [max(mean(auc(5 cross validation)))]
    lgb_best_model = lgb_rscv.best_estimator_
    #Get predict_class(y_pred) and predict_probality_for_case(y_prob) of TestSet
    y_pred = lgb_best_model.predict(x_test)
    y_prob = lgb_best_model.predict_proba(x_test)[:,1]

    #Get model performance
    model_perf = evaluate_performance(y_test,y_pred,y_prob)
    #Output result of evaluation
    eval_output(model_perf,path)
    #You can make bar plot consisted of accuracy,sensitivity,specificity,auroc,f1 score,MCC,precision,recall,auprc according to the "Evaluate_Result_TestSet.txt"
    # Plot AUROC
    plot_AUROC(model_perf,path)

    #save model
    joblib.dump(lgb_best_model,os.path.join(path,"best_LightGBM_model.pkl"))
    return model_perf

# %%
def catboost_classification_random_search(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_iter=100):
    print("\n*** CatBoost (RandomizedSearchCV) ***")
    
    # CatBoost hyperparameter space for random search
    cb_param_dict = {
        'iterations': randint(100, 2000),
        'learning_rate': uniform(0.01, 0.19),  # 0.01 to 0.2
        'depth': randint(4, 10),  # Tree depth
        'l2_leaf_reg': uniform(0, 10),  # L2 regularization
        'border_count': [32, 64, 128, 256],  # Number of splits for numerical features
        'bagging_temperature': uniform(0, 1),  # Bayesian bagging
        'random_strength': uniform(0, 10),  # Score randomization
        'scale_pos_weight': [1, 2, 5, 10],  # For imbalanced datasets
        'loss_function': ['Logloss'],
        'eval_metric': ['AUC'],
        'random_state': [SEED],
        'verbose': [False]
    }
    
    # Initiate model
    cb_model = CatBoostClassifier()
    
    # RandomizedSearchCV with 5-fold cross-validation
    cb_rscv = RandomizedSearchCV(
        cb_model,
        cb_param_dict,
        n_iter=n_iter,
        cv=5,
        verbose=1,
        scoring='roc_auc',
        random_state=SEED,
        n_jobs=30 # Use all available CPU cores
    )
    
    cb_rscv.fit(x_train, y_train, eval_set=(x_test, y_test), early_stopping_rounds=50, verbose=False)
    
    # Output path
    path = os.path.join(output_dir, "CatBoost_RandomSearch")
    if not os.path.exists(path):
        os.mkdir(path)
    
    # Model performance (AUROC) on cross-validation dataset
    cb_cv_perf = np.array([cb_rscv.cv_results_[f"split{i}_test_score"] for i in range(5)])[:, cb_rscv.best_index_]
    
    # Get best model
    cb_best_model = cb_rscv.best_estimator_
    
    # Get predictions
    y_pred = cb_best_model.predict(x_test)
    y_prob = cb_best_model.predict_proba(x_test)[:, 1]
    
    # Get model performance (assuming evaluate_performance is defined)
    model_perf = evaluate_performance(y_test, y_pred, y_prob)
    
    # Output result of evaluation (assuming eval_output is defined)
    eval_output(model_perf, path)
    
    # Plot AUROC (assuming plot_AUROC is defined)
    plot_AUROC(model_perf, path)
    
    # Save model
    joblib.dump(cb_best_model, os.path.join(path, "best_CatBoost_model_random_search.pkl"))
    
    return model_perf

# %%
def RF_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_ter = 100):
    # Construct Random Forest model
    print("\n*** Random Forest  ***")

    # Random Forest params
    rf_param_dict = {
        "n_estimators": [10, 50, 100, 200, 500, 1000],
        "max_depth": [None, 10, 20, 50, 100],
        "min_samples_split": [2, 5, 10],
        "min_samples_leaf": [1, 2, 4],
        "max_features": ['auto', 'sqrt', 'log2'],
        "random_state": [SEED]
    }

    # Initiate model
    rf_model = RandomForestClassifier()
    # Adjust hyper-parameters with 5-fold cross-validation
    rf_rscv = RandomizedSearchCV(rf_model, rf_param_dict, n_iter=n_ter, cv=5, verbose=0,
                                scoring="roc_auc", random_state=SEED, n_jobs=30)
    rf_rscv.fit(x_train, y_train)

    # Evaluate best Random Forest model
    # Output path
    path = os.path.join(output_dir, "RandomForest")
    if not (os.path.exists(path)):
        os.mkdir(path)

    # Model performance(AUROC) on cross-validation dataset
    rf_cv_perf = np.array([rf_rscv.cv_results_["split%s_test_score" % str(i)] for i in range(5)])[:, rf_rscv.best_index_]

    # Get best model with score [max(mean(auc(5 cross validation)))]
    rf_best_model = rf_rscv.best_estimator_
    # Get predict_class(y_pred) and predict_probability_for_case(y_prob) of TestSet
    y_pred = rf_best_model.predict(x_test)
    y_prob = rf_best_model.predict_proba(x_test)[:, 1]

    # Get model performance
    model_perf = evaluate_performance(y_test, y_pred, y_prob)
    # Output result of evaluation
    eval_output(model_perf, path)
    # You can make bar plot consisted of accuracy, sensitivity, specificity, auroc, f1 score, MCC, precision, recall, auprc according to the "Evaluate_Result_TestSet.txt"
    # Plot AUROC
    plot_AUROC(model_perf, path)

    # Save model
    joblib.dump(rf_best_model, os.path.join(path, "best_RandomForest_model.pkl"))
    return model_perf


# %%
def ngboost_classification_random_search(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_iter=50):
    print("\n*** NGBoost (RandomizedSearchCV) ***")
    
    # NGBoost hyperparameter space for random search
    ngb_param_dict = {
        "n_estimators": randint(50, 500),       # Number of boosting iterations
        "learning_rate": uniform(0.01, 0.19),    # Step size shrinkage (0.01 to 0.2)
        "minibatch_frac": uniform(0.1, 0.9),     # Fraction of data used for each iteration
        "col_sample": uniform(0.5, 0.5),        # Feature subsampling ratio (0.5 to 1.0)
        "natural_gradient": [True, False],      # Whether to use natural gradient
        "verbose": [False],
        "random_state": [np.random.RandomState(SEED)] 
    }
    
    # Initiate model (NGBoost classifier)
    ngb_model = NGBClassifier(Dist=k_categorical(2))  # Binary classification
    
    # RandomizedSearchCV with 3-fold CV (NGBoost is slower than XGBoost/CatBoost)
    ngb_rscv = RandomizedSearchCV(
        ngb_model,
        ngb_param_dict,
        n_iter=n_iter,
        cv=5,                   # Fewer folds due to computational cost
        verbose=1,
        scoring="roc_auc",
        random_state=SEED,
        n_jobs=30         
    )
    
    # Fit model (no early stopping in NGBoost)
    ngb_rscv.fit(x_train, y_train)
    
    # Output directory
    path = os.path.join(output_dir, "NGBoost_RandomSearch")
    if not os.path.exists(path):
        os.mkdir(path)
    
    # Get best model
    ngb_best_model = ngb_rscv.best_estimator_
    
    # Predictions
    y_pred = ngb_best_model.predict(x_test)            # Class predictions
    y_prob = ngb_best_model.predict_proba(x_test)[:, 1] # Probability estimates
    
    # Evaluate performance (assuming evaluate_performance() is defined)
    model_perf = evaluate_performance(y_test, y_pred, y_prob)
    
    # Save outputs
    eval_output(model_perf, path)                      # Save metrics
    plot_AUROC(model_perf, path)                      # Plot AUROC
    joblib.dump(ngb_best_model, os.path.join(path, "best_NGBoost_model.pkl"))
    
    return model_perf

# %%
# Random seed
SEED = 100
random.seed(SEED)
np.random.seed(SEED)

warnings.filterwarnings(action='ignore')

# Output dir
output_dir = "/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/ML_models/circRNA_ML_Model_tridivided_Output"
if not (os.path.exists(output_dir)):
    os.mkdir(output_dir)

# %%
# 先划分数据集后计算
dataset = pd.read_csv(
    '/lulabdata3/huangkeyun/zhangys/RNA_locator/ML_python_scripts/reference_preprocessing/circRNA/output_with_sequences.csv',
    sep='\t',
    index_col=False
)
dataset_filtered = dataset
# 将circRNA序列首尾拼接保证kmer的环形分析
add_len = 4 
dataset_filtered["Sequence"] = dataset_filtered["Sequence"].apply(
    lambda x: ligate_sequence(x, add_len)
)


dataset_filtered = balance_dataset_by_tag(dataset_filtered, tag_column='tag', random_state=42)

# 3. 划分训练集、验证集和测试集（6:2:2）
train_df, temp_df = train_test_split(dataset_filtered, test_size=0.4, random_state=SEED, stratify=dataset_filtered['tag'])
val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=SEED, stratify=temp_df['tag'])

# 4. 分别提取训练集、验证集和测试集的 k-mer 特征
df_kmer_train, df_kmer_train_raw = _count_kmer(train_df, 345)
df_kmer_val, df_kmer_val_raw = _count_kmer(val_df, 345)
df_kmer_test, df_kmer_test_raw = _count_kmer(test_df, 345)

# 保存特征
df_kmer_train.to_csv(os.path.join(output_dir, "train_kmer345_freq.tsv"), sep='\t')
df_kmer_train_raw.to_csv(os.path.join(output_dir, "train_kmer345_rawcount.tsv"), sep='\t')
df_kmer_val.to_csv(os.path.join(output_dir, "val_kmer345_freq.tsv"), sep='\t')
df_kmer_val_raw.to_csv(os.path.join(output_dir, "val_kmer345_rawcount.tsv"), sep='\t')
df_kmer_test.to_csv(os.path.join(output_dir, "test_kmer345_freq.tsv"), sep='\t')
df_kmer_test_raw.to_csv(os.path.join(output_dir, "test_kmer345_rawcount.tsv"), sep='\t')


# %%
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# 构建 SeqRecord 列表
fasta_records = [
    SeqRecord(Seq(row["Sequence"]), id=row["RNA_Symbol"], description=row["Subcellular_Localization"])
    for _, row in test_df.iterrows()
]

# 保存为 fasta 文件
fasta_path = os.path.join(output_dir, "test_set_sequences.fasta")
SeqIO.write(fasta_records, fasta_path, "fasta")

# %%
# 获取训练、验证和测试数据的特征矩阵 + 标签
x_train = df_kmer_train.drop(columns=["RNA_Symbol"]).values
x_val = df_kmer_val.drop(columns=["RNA_Symbol"]).values
x_test = df_kmer_test.drop(columns=["RNA_Symbol"]).values

y_train = train_df["tag"].values
y_val = val_df["tag"].values
y_test = test_df["tag"].values

# 6. 缺失值处理
imputer = SimpleImputer(strategy='mean')
x_train = imputer.fit_transform(x_train)
x_val = imputer.transform(x_val)
x_test = imputer.transform(x_test)

# %%
# kmer = 345, model = LR
#Construct Logistic Regression model
LR_perf = logistic_regression_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

SVM_perf = SVM_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

RF_perf = RF_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

lightGBM_perf = lightGBM_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

CatBoost_perf = catboost_classification_random_search(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

NGBoost_perf = ngboost_classification_random_search(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

