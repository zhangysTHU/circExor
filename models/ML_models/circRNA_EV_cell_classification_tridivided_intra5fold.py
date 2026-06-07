    

               

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

    

                                                   

                                                                

def _count_kmer(Dataset, k):               

    

                  

    dataset = copy.deepcopy(Dataset)

                            

    nucleotide = ['A', 'C', 'G', 'T']

    

                     

              

    five = list(itertools.product(nucleotide, repeat=5))

    pentamer = [''.join(n) for n in five]

    

              

    four = list(itertools.product(nucleotide, repeat=4))

    tetramer = [''.join(n) for n in four]

             

    three = list(itertools.product(nucleotide, repeat=3))

    threemer = [''.join(n) for n in three]

    

                                                              

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

                                   

    for mer in table_kmer.keys():

        table_kmer[mer] = dataset["Sequence"].apply(lambda x: x.count(mer))

    

                                                                       

    rawcount_kmer_df = pd.DataFrame(table_kmer)

    df1_rawcount = pd.concat([rawcount_kmer_df, dataset["RNA_Symbol"]], axis=1)

    df1_rawcount.index = dataset["tag"]

                                                                    

    freq_kmer_df = rawcount_kmer_df.apply(lambda x: x / x.sum(), axis=1)

    df1 = pd.concat([freq_kmer_df, dataset["RNA_Symbol"]], axis=1)

    df1.index = dataset["tag"]

    return df1, df1_rawcount

    

                              

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

    

                     

def plot_AUROC(model_perf,path):

                                    

    roc_auc = model_perf["auroc"]

    fpr,tpr,threshold = model_perf["auroc_curve"]

                      

    temp_df = pd.DataFrame({"FPR":fpr,"TPR":tpr})

    temp_df.to_csv(os.path.join(path,"AUROC_info.txt"),header = True,index = False, sep = '\t')

         

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

    

def logistic_regression_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42):

    """
    Run binary Logistic Regression with hyperparameter tuning, training, evaluation, and result saving.

    Parameters：
    x_train : Training features
    y_train : Training labels
    x_test : Test features
    y_test : Test labels
    output_dir : Output directory path
    SEED : Random seed; default is 42
    """

    print("\n*** Logistic Regression  ***")

                                

    lr_param_dict = {

        "penalty": ["l2"],

        "C": [1e-3, 5e-3, 1e-2, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000],

        "solver": ["liblinear"],

        "random_state": [SEED]

    }

                    

    lr_model = LogisticRegression()

                                                          

    lr_rscv = RandomizedSearchCV(lr_model, lr_param_dict, n_iter=100, cv=5, verbose=0,

                                 scoring="roc_auc", random_state=SEED, n_jobs=30)

    lr_rscv.fit(x_train, y_train)

                 

    path = os.path.join(output_dir, "LogisticRegression")

    if not os.path.exists(path):

        os.makedirs(path, exist_ok=True)

                                              

    from sklearn.model_selection import cross_validate

    import joblib

    import pandas as pd

    

                                     

    best_params = lr_rscv.best_params_

    cv_template_model = LogisticRegression(**best_params)

    cv_results = cross_validate(cv_template_model, x_train, y_train, cv=5, return_estimator=True)

    

    cv_metrics_list = []            

    

    print("\n--- 5-Fold Cross Validation Models Performance ---")

    for i, model in enumerate(cv_results['estimator']):

              

        model_path = os.path.join(path, f"cv_fold_{i+1}_LR_model.pkl")

        joblib.dump(model, model_path)

        

                     

        y_pred_cv = model.predict(x_test)

        y_prob_cv = model.predict_proba(x_test)[:, 1]

        model_perf_cv = evaluate_performance(y_test, y_pred_cv, y_prob_cv)

        

               

        cv_metrics_list.append({

            "Fold": f"Fold_{i+1}",

            "AUROC": model_perf_cv['auroc'],

            "AUPRC": model_perf_cv['auprc'],

            "Accuracy": model_perf_cv['accuracy'],

            "Precision": model_perf_cv['precision'],

            "Recall": model_perf_cv['recall'],

            "F1_Score": model_perf_cv['f1'],

            "MCC": model_perf_cv['mcc']

        })

        

                                

        print(f"Fold {i+1} Model Performance:")

        print(f"  AUROC     : {model_perf_cv['auroc']:.4f}")

        print(f"  AUPRC     : {model_perf_cv['auprc']:.4f}")

        print(f"  Accuracy  : {model_perf_cv['accuracy']:.4f}")

        print(f"  Precision : {model_perf_cv['precision']:.4f}")

        print(f"  Recall    : {model_perf_cv['recall']:.4f}")

        print(f"  F1 Score  : {model_perf_cv['f1']:.4f}")

        print(f"  MCC       : {model_perf_cv['mcc']:.4f}")

        print("-" * 50)

        

                   

    cv_metrics_df = pd.DataFrame(cv_metrics_list)

    cv_metrics_df.to_csv(os.path.join(path, "5fold_cv_metrics.csv"), index=False)

                                                                   

                                                           

    lr_cv_perf = np.array([lr_rscv.cv_results_["split%s_test_score" % str(i)] for i in range(5)])[:, lr_rscv.best_index_]

                                                                    

    lr_best_model = lr_rscv.best_estimator_

                                                                                 

    y_pred = lr_best_model.predict(x_test)

    y_prob = lr_best_model.predict_proba(x_test)[:, 1]

                           

    model_perf = evaluate_performance(y_test, y_pred, y_prob)

                                 

    eval_output(model_perf, path)

                

    joblib.dump(lr_best_model, os.path.join(path, "best_LogisticRegression_model.pkl"))

    return model_perf

    

def SVM_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42):

                        

    print("\n*** SVM ***")

                

    SVM_param_dict = {

        'kernel':('linear', 'rbf'), 

        'C':[0.01,0.1,1,10, 100], 

        'gamma':[0.001, 0.005, 0.1 ,0.5,1, 2],

        "probability":[True],

        "random_state":[SEED]

    }

                   

    SVM_model = svm.SVC()

                                                         

    SVM_rscv = RandomizedSearchCV(SVM_model, SVM_param_dict, n_iter=100,cv = 5,verbose = 0,

                            scoring = "roc_auc",random_state=SEED,n_jobs =30)

    SVM_rscv.fit(x_train, y_train) 

                            

                

    path = os.path.join(output_dir,"SVM")

    if not (os.path.exists(path)):

        os.makedirs(path, exist_ok=True)

        

                                                 

    from sklearn.model_selection import cross_validate

    import joblib

    import pandas as pd

    

                                     

    best_params = SVM_rscv.best_params_

    cv_template_model = svm.SVC(**best_params)

    cv_results = cross_validate(cv_template_model, x_train, y_train, cv=5, return_estimator=True)

    

    cv_metrics_list = []            

    

    print("\n--- 5-Fold Cross Validation Models Performance ---")

    for i, model in enumerate(cv_results['estimator']):

              

        model_path = os.path.join(path, f"cv_fold_{i+1}_SVM_model.pkl")

        joblib.dump(model, model_path)

        

                     

        y_pred_cv = model.predict(x_test)

        y_prob_cv = model.predict_proba(x_test)[:, 1]

        model_perf_cv = evaluate_performance(y_test, y_pred_cv, y_prob_cv)

        

               

        cv_metrics_list.append({

            "Fold": f"Fold_{i+1}",

            "AUROC": model_perf_cv['auroc'],

            "AUPRC": model_perf_cv['auprc'],

            "Accuracy": model_perf_cv['accuracy'],

            "Precision": model_perf_cv['precision'],

            "Recall": model_perf_cv['recall'],

            "F1_Score": model_perf_cv['f1'],

            "MCC": model_perf_cv['mcc']

        })

        

                    

        print(f"Fold {i+1} Model Performance:")

        print(f"  AUROC     : {model_perf_cv['auroc']:.4f}")

        print(f"  AUPRC     : {model_perf_cv['auprc']:.4f}")

        print(f"  Accuracy  : {model_perf_cv['accuracy']:.4f}")

        print(f"  Precision : {model_perf_cv['precision']:.4f}")

        print(f"  Recall    : {model_perf_cv['recall']:.4f}")

        print(f"  F1 Score  : {model_perf_cv['f1']:.4f}")

        print(f"  MCC       : {model_perf_cv['mcc']:.4f}")

        print("-" * 50)

        

                   

    cv_metrics_df = pd.DataFrame(cv_metrics_list)

    cv_metrics_df.to_csv(os.path.join(path, "5fold_cv_metrics.csv"), index=False)

                                                                   

                                                          

    SVM_cv_perf = np.array([ SVM_rscv.cv_results_["split%s_test_score"%str(i)] for i in range(5)])[:,SVM_rscv.best_index_]

                                                                   

    svm_best_model = SVM_rscv.best_estimator_

                                                                                

    y_pred = svm_best_model.predict(x_test)

    y_prob = svm_best_model.predict_proba(x_test)[:,1]

                          

    model_perf = evaluate_performance(y_test,y_pred,y_prob)

                                

    eval_output(model_perf,path)

                                                                                                                                                                 

                

                                 

               

    joblib.dump(svm_best_model,os.path.join(path,"best_SVM_model.pkl"))

               

                                                                          

                                          

    return model_perf

    

def lightGBM_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_iter = 100):

    print("\n*** LightGBM ***")

                     

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

                   

    lgb_model = lgb.LGBMClassifier()

                                                         

    lgb_rscv = RandomizedSearchCV(lgb_model, lgb_param_dict, n_iter=n_iter,cv = 5,verbose = 0,

                            scoring = "roc_auc",random_state=SEED,n_jobs = 30)

    lgb_rscv.fit(x_train, y_train)   

                                 

                

    path = os.path.join(output_dir,"LightGBM")

    if not (os.path.exists(path)):

        os.makedirs(path, exist_ok=True)

        

                                                 

    from sklearn.model_selection import cross_validate

    import joblib

    import pandas as pd

    

                                     

    best_params = lgb_rscv.best_params_

    cv_template_model = lgb.LGBMClassifier(**best_params)

    cv_results = cross_validate(cv_template_model, x_train, y_train, cv=5, return_estimator=True)

    

    cv_metrics_list = []            

    

    print("\n--- 5-Fold Cross Validation Models Performance ---")

    for i, model in enumerate(cv_results['estimator']):

              

        model_path = os.path.join(path, f"cv_fold_{i+1}_LightGBM_model.pkl")

        joblib.dump(model, model_path)

        

                     

        y_pred_cv = model.predict(x_test)

        y_prob_cv = model.predict_proba(x_test)[:, 1]

        model_perf_cv = evaluate_performance(y_test, y_pred_cv, y_prob_cv)

        

               

        cv_metrics_list.append({

            "Fold": f"Fold_{i+1}",

            "AUROC": model_perf_cv['auroc'],

            "AUPRC": model_perf_cv['auprc'],

            "Accuracy": model_perf_cv['accuracy'],

            "Precision": model_perf_cv['precision'],

            "Recall": model_perf_cv['recall'],

            "F1_Score": model_perf_cv['f1'],

            "MCC": model_perf_cv['mcc']

        })

        

                    

        print(f"Fold {i+1} Model Performance:")

        print(f"  AUROC     : {model_perf_cv['auroc']:.4f}")

        print(f"  AUPRC     : {model_perf_cv['auprc']:.4f}")

        print(f"  Accuracy  : {model_perf_cv['accuracy']:.4f}")

        print(f"  Precision : {model_perf_cv['precision']:.4f}")

        print(f"  Recall    : {model_perf_cv['recall']:.4f}")

        print(f"  F1 Score  : {model_perf_cv['f1']:.4f}")

        print(f"  MCC       : {model_perf_cv['mcc']:.4f}")

        print("-" * 50)

        

                   

    cv_metrics_df = pd.DataFrame(cv_metrics_list)

    cv_metrics_df.to_csv(os.path.join(path, "5fold_cv_metrics.csv"), index=False)

                                                                   

                                                          

    lgb_cv_perf = np.array([ lgb_rscv.cv_results_["split%s_test_score"%str(i)] for i in range(5)])[:,lgb_rscv.best_index_]

                                                                   

    lgb_best_model = lgb_rscv.best_estimator_

                                                                                

    y_pred = lgb_best_model.predict(x_test)

    y_prob = lgb_best_model.predict_proba(x_test)[:,1]

                          

    model_perf = evaluate_performance(y_test,y_pred,y_prob)

                                

    eval_output(model_perf,path)

                                                                                                                                                                 

                

                                 

               

    joblib.dump(lgb_best_model,os.path.join(path,"best_LightGBM_model.pkl"))

    return model_perf

    

def catboost_classification_random_search(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_iter=100):

    print("\n*** CatBoost (RandomizedSearchCV) ***")

    

                                                     

    cb_param_dict = {

        'iterations': randint(100, 2000),

        'learning_rate': uniform(0.01, 0.19),               

        'depth': randint(4, 10),              

        'l2_leaf_reg': uniform(0, 10),                     

        'border_count': [32, 64, 128, 256],                                           

        'bagging_temperature': uniform(0, 1),                    

        'random_strength': uniform(0, 10),                       

        'scale_pos_weight': [1, 2, 5, 10],                           

        'loss_function': ['Logloss'],

        'eval_metric': ['AUC'],

        'random_state': [SEED],

        'verbose': [False]

    }

    

                    

    cb_model = CatBoostClassifier()

    

                                                     

    cb_rscv = RandomizedSearchCV(

        cb_model,

        cb_param_dict,

        n_iter=n_iter,

        cv=5,

        verbose=1,

        scoring='roc_auc',

        random_state=SEED,

        n_jobs=30                              

    )

    

    cb_rscv.fit(x_train, y_train, eval_set=(x_test, y_test), early_stopping_rounds=50, verbose=False)

    

                 

    path = os.path.join(output_dir, "CatBoost_RandomSearch")

    if not os.path.exists(path):

        os.makedirs(path, exist_ok=True)

        

                                                 

    from sklearn.model_selection import cross_validate

    import joblib

    import pandas as pd

    

                                     

    best_params = cb_rscv.best_params_

    cv_template_model = CatBoostClassifier(**best_params)

    cv_template_model.set_params(verbose=False)                   

    

    cv_results = cross_validate(cv_template_model, x_train, y_train, cv=5, return_estimator=True)

    

    cv_metrics_list = []            

    

    print("\n--- 5-Fold Cross Validation Models Performance ---")

    for i, model in enumerate(cv_results['estimator']):

              

        model_path = os.path.join(path, f"cv_fold_{i+1}_CatBoost_model.pkl")

        joblib.dump(model, model_path)

        

                     

        y_pred_cv = model.predict(x_test)

        y_prob_cv = model.predict_proba(x_test)[:, 1]

        model_perf_cv = evaluate_performance(y_test, y_pred_cv, y_prob_cv)

        

               

        cv_metrics_list.append({

            "Fold": f"Fold_{i+1}",

            "AUROC": model_perf_cv['auroc'],

            "AUPRC": model_perf_cv['auprc'],

            "Accuracy": model_perf_cv['accuracy'],

            "Precision": model_perf_cv['precision'],

            "Recall": model_perf_cv['recall'],

            "F1_Score": model_perf_cv['f1'],

            "MCC": model_perf_cv['mcc']

        })

        

                    

        print(f"Fold {i+1} Model Performance:")

        print(f"  AUROC     : {model_perf_cv['auroc']:.4f}")

        print(f"  AUPRC     : {model_perf_cv['auprc']:.4f}")

        print(f"  Accuracy  : {model_perf_cv['accuracy']:.4f}")

        print(f"  Precision : {model_perf_cv['precision']:.4f}")

        print(f"  Recall    : {model_perf_cv['recall']:.4f}")

        print(f"  F1 Score  : {model_perf_cv['f1']:.4f}")

        print(f"  MCC       : {model_perf_cv['mcc']:.4f}")

        print("-" * 50)

        

                   

    cv_metrics_df = pd.DataFrame(cv_metrics_list)

    cv_metrics_df.to_csv(os.path.join(path, "5fold_cv_metrics.csv"), index=False)

                                                                   

    

                                                           

    cb_cv_perf = np.array([cb_rscv.cv_results_[f"split{i}_test_score"] for i in range(5)])[:, cb_rscv.best_index_]

    

                    

    cb_best_model = cb_rscv.best_estimator_

    

                     

    y_pred = cb_best_model.predict(x_test)

    y_prob = cb_best_model.predict_proba(x_test)[:, 1]

    

                                                                      

    model_perf = evaluate_performance(y_test, y_pred, y_prob)

    

                                                                   

    eval_output(model_perf, path)

    

                                                 

    plot_AUROC(model_perf, path)

    

                

    joblib.dump(cb_best_model, os.path.join(path, "best_CatBoost_model_random_search.pkl"))

    

    return model_perf

    

def RF_classification(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_ter = 100):

                                   

    print("\n*** Random Forest  ***")

                          

    rf_param_dict = {

        "n_estimators": [10, 50, 100, 200, 500, 1000],

        "max_depth": [None, 10, 20, 50, 100],

        "min_samples_split": [2, 5, 10],

        "min_samples_leaf": [1, 2, 4],

        "max_features": ['auto', 'sqrt', 'log2'],

        "random_state": [SEED]

    }

                    

    rf_model = RandomForestClassifier()

                                                          

    rf_rscv = RandomizedSearchCV(rf_model, rf_param_dict, n_iter=n_ter, cv=5, verbose=0,

                                scoring="roc_auc", random_state=SEED, n_jobs=30)

    rf_rscv.fit(x_train, y_train)

                                       

                 

    path = os.path.join(output_dir, "RandomForest")

    if not (os.path.exists(path)):

        os.makedirs(path, exist_ok=True)

        

                                                 

    from sklearn.model_selection import cross_validate

    import joblib

    import pandas as pd

    

                                     

    best_params = rf_rscv.best_params_

    cv_template_model = RandomForestClassifier(**best_params)

    cv_results = cross_validate(cv_template_model, x_train, y_train, cv=5, return_estimator=True)

    

    cv_metrics_list = []            

    

    print("\n--- 5-Fold Cross Validation Models Performance ---")

    for i, model in enumerate(cv_results['estimator']):

              

        model_path = os.path.join(path, f"cv_fold_{i+1}_RF_model.pkl")

        joblib.dump(model, model_path)

        

                     

        y_pred_cv = model.predict(x_test)

        y_prob_cv = model.predict_proba(x_test)[:, 1]

        model_perf_cv = evaluate_performance(y_test, y_pred_cv, y_prob_cv)

        

               

        cv_metrics_list.append({

            "Fold": f"Fold_{i+1}",

            "AUROC": model_perf_cv['auroc'],

            "AUPRC": model_perf_cv['auprc'],

            "Accuracy": model_perf_cv['accuracy'],

            "Precision": model_perf_cv['precision'],

            "Recall": model_perf_cv['recall'],

            "F1_Score": model_perf_cv['f1'],

            "MCC": model_perf_cv['mcc']

        })

        

                    

        print(f"Fold {i+1} Model Performance:")

        print(f"  AUROC     : {model_perf_cv['auroc']:.4f}")

        print(f"  AUPRC     : {model_perf_cv['auprc']:.4f}")

        print(f"  Accuracy  : {model_perf_cv['accuracy']:.4f}")

        print(f"  Precision : {model_perf_cv['precision']:.4f}")

        print(f"  Recall    : {model_perf_cv['recall']:.4f}")

        print(f"  F1 Score  : {model_perf_cv['f1']:.4f}")

        print(f"  MCC       : {model_perf_cv['mcc']:.4f}")

        print("-" * 50)

        

                   

    cv_metrics_df = pd.DataFrame(cv_metrics_list)

    cv_metrics_df.to_csv(os.path.join(path, "5fold_cv_metrics.csv"), index=False)

                                                                   

                                                          

    rf_cv_perf = np.array([rf_rscv.cv_results_["split%s_test_score" % str(i)] for i in range(5)])[:, rf_rscv.best_index_]

                                                                    

    rf_best_model = rf_rscv.best_estimator_

                                                                                   

    y_pred = rf_best_model.predict(x_test)

    y_prob = rf_best_model.predict_proba(x_test)[:, 1]

                           

    model_perf = evaluate_performance(y_test, y_pred, y_prob)

                                 

    eval_output(model_perf, path)

                                                                                                                                                                          

                

    plot_AUROC(model_perf, path)

                

    joblib.dump(rf_best_model, os.path.join(path, "best_RandomForest_model.pkl"))

    return model_perf

    

def ngboost_classification_random_search(x_train, y_train, x_test, y_test, output_dir, SEED=42, n_iter=50):

    print("\n*** NGBoost (RandomizedSearchCV) ***")

    

                                                    

    ngb_param_dict = {

        "n_estimators": randint(50, 500),                                      

        "learning_rate": uniform(0.01, 0.19),                                       

        "minibatch_frac": uniform(0.1, 0.9),                                               

        "col_sample": uniform(0.5, 0.5),                                                

        "natural_gradient": [True, False],                                       

        "verbose": [False],

        "random_state": [np.random.RandomState(SEED)] 

    }

    

                                         

    ngb_model = NGBClassifier(Dist=k_categorical(2))                         

    

                                       

    ngb_rscv = RandomizedSearchCV(

        ngb_model,

        ngb_param_dict,

        n_iter=n_iter,

        cv=5,                   

        verbose=1,

        scoring="roc_auc",

        random_state=SEED,

        n_jobs=30         

    )

    

                                              

    ngb_rscv.fit(x_train, y_train)

    

                      

    path = os.path.join(output_dir, "NGBoost_RandomSearch")

    if not os.path.exists(path):

        os.makedirs(path, exist_ok=True)

        

                                                 

    from sklearn.model_selection import cross_validate

    import joblib

    import pandas as pd

    

                                     

    best_params = ngb_rscv.best_params_

                                                        

    cv_template_model = NGBClassifier(Dist=k_categorical(2), **best_params)

    cv_results = cross_validate(cv_template_model, x_train, y_train, cv=5, return_estimator=True)

    

    cv_metrics_list = []            

    

    print("\n--- 5-Fold Cross Validation Models Performance ---")

    for i, model in enumerate(cv_results['estimator']):

              

        model_path = os.path.join(path, f"cv_fold_{i+1}_NGBoost_model.pkl")

        joblib.dump(model, model_path)

        

                     

        y_pred_cv = model.predict(x_test)

        y_prob_cv = model.predict_proba(x_test)[:, 1]

        model_perf_cv = evaluate_performance(y_test, y_pred_cv, y_prob_cv)

        

               

        cv_metrics_list.append({

            "Fold": f"Fold_{i+1}",

            "AUROC": model_perf_cv['auroc'],

            "AUPRC": model_perf_cv['auprc'],

            "Accuracy": model_perf_cv['accuracy'],

            "Precision": model_perf_cv['precision'],

            "Recall": model_perf_cv['recall'],

            "F1_Score": model_perf_cv['f1'],

            "MCC": model_perf_cv['mcc']

        })

        

                    

        print(f"Fold {i+1} Model Performance:")

        print(f"  AUROC     : {model_perf_cv['auroc']:.4f}")

        print(f"  AUPRC     : {model_perf_cv['auprc']:.4f}")

        print(f"  Accuracy  : {model_perf_cv['accuracy']:.4f}")

        print(f"  Precision : {model_perf_cv['precision']:.4f}")

        print(f"  Recall    : {model_perf_cv['recall']:.4f}")

        print(f"  F1 Score  : {model_perf_cv['f1']:.4f}")

        print(f"  MCC       : {model_perf_cv['mcc']:.4f}")

        print("-" * 50)

        

                   

    cv_metrics_df = pd.DataFrame(cv_metrics_list)

    cv_metrics_df.to_csv(os.path.join(path, "5fold_cv_metrics.csv"), index=False)

                                                                   

    

                                                           

    ngb_cv_perf = np.array([ngb_rscv.cv_results_["split%s_test_score" % str(i)] for i in range(5)])[:, ngb_rscv.best_index_]

                    

    ngb_best_model = ngb_rscv.best_estimator_

    

                 

    y_pred = ngb_best_model.predict(x_test)                                

    y_prob = ngb_best_model.predict_proba(x_test)[:, 1]                        

    

                          

    model_perf = evaluate_performance(y_test, y_pred, y_prob)

    

                  

    eval_output(model_perf, path)                                     

                                                                    

    joblib.dump(ngb_best_model, os.path.join(path, "best_NGBoost_model.pkl"))

    

    return model_perf

    

             

SEED = 100

random.seed(SEED)

np.random.seed(SEED)

warnings.filterwarnings(action='ignore')

            

output_dir = "circExor/ML_models/circRNA_ML_Model_tridivided_intra5fold_Output"

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

                              

df_kmer_train, df_kmer_train_raw = _count_kmer(train_df, 345)

df_kmer_val, df_kmer_val_raw = _count_kmer(val_df, 345)

df_kmer_test, df_kmer_test_raw = _count_kmer(test_df, 345)

      

df_kmer_train.to_csv(os.path.join(output_dir, "train_kmer345_freq.tsv"), sep='\t')

df_kmer_train_raw.to_csv(os.path.join(output_dir, "train_kmer345_rawcount.tsv"), sep='\t')

df_kmer_val.to_csv(os.path.join(output_dir, "val_kmer345_freq.tsv"), sep='\t')

df_kmer_val_raw.to_csv(os.path.join(output_dir, "val_kmer345_rawcount.tsv"), sep='\t')

df_kmer_test.to_csv(os.path.join(output_dir, "test_kmer345_freq.tsv"), sep='\t')

df_kmer_test_raw.to_csv(os.path.join(output_dir, "test_kmer345_rawcount.tsv"), sep='\t')

    

dataset_filtered['tag'].value_counts()

    

test_df['tag'].value_counts()

    

from Bio.SeqRecord import SeqRecord

from Bio.Seq import Seq

from Bio import SeqIO

import os

fasta_records = [

    SeqRecord(

        Seq(row["Sequence"]), 

        id=f"{row['RNA_Symbol']}_{row['tag']}", 

        description=""

    )

    for _, row in test_df.iterrows()

]

              

fasta_path = os.path.join(output_dir, "test_set_sequences.fasta")

SeqIO.write(fasta_records, fasta_path, "fasta")

    

                        

x_train = df_kmer_train.drop(columns=["RNA_Symbol"]).values

x_val = df_kmer_val.drop(columns=["RNA_Symbol"]).values

x_test = df_kmer_test.drop(columns=["RNA_Symbol"]).values

y_train = train_df["tag"].values

y_val = val_df["tag"].values

y_test = test_df["tag"].values

          

imputer = SimpleImputer(strategy='mean')

x_train = imputer.fit_transform(x_train)

x_val = imputer.transform(x_val)

x_test = imputer.transform(x_test)

    

                        

                                    

LR_perf = logistic_regression_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

    

SVM_perf = SVM_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

    

RF_perf = RF_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

    

lightGBM_perf = lightGBM_classification(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

    

CatBoost_perf = catboost_classification_random_search(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

    

NGBoost_perf = ngboost_classification_random_search(x_test=x_val, x_train=x_train,y_test=y_val,y_train=y_train, output_dir=output_dir)

