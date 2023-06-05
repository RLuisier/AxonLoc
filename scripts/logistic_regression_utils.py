import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys
from sklearn.base import clone
from sklearn.model_selection import RandomizedSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.feature_selection import chi2
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
sys.path.insert(1, "/".join(os.getcwd().split('/')[:-1]))
from scripts.ml_utils import *


def get_best_logit_model(X_train, y_train, 
                         parameters_lr={'solver': ['lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga'], 
                                        'penalty': ['l2', 'l1', 'elasticnet'],
                                        'C': [100, 10, 1.0, 0.1, 0.01]}, class_weight='balanced', 
                         scoring='f1', random_state=1, max_iter=1000,  n_jobs=-1,  cv=20,  verbose=1,  n_iter=75):
    
    """Returns the best model logistic model found by the RandomSearchCV.
    
    Args:
        X_train (array-like, sparse matrix) of shape (n_samples, n_features): Training vector, 
        where n_samples is the number of samples and n_features is the number of features.
        y_train (array-like) of shape (n_samples,): Target vector relative to X.
        parameters_lr (dict or list of dicts, optional): Dicitonnary with parameters names (str) as keys and distributions or lists of parameters to try. 
        Defaults to {'solver': ['lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga'], 'penalty': ['l2', 'l1', 'elasticnet'], 'C': [100, 10, 1.0, 0.1, 0.01]}.
        class_weight (dict or 'balanced'): Weights associated with classes in the form {class_label: weight}. The "balanced" mode
        uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data. Defaults to 'balanced'.
        scoring (str, optional): Scoring method. Defaults to 'f1'.
        random_state (int, optional): Used when solver=='sag', 'saga', or 'liblinear' to shuffle the data. Defaults to 1.
        max_iter (int, optional): Maximun number of iterations taken for the solvers to converge. Defaults to 1000.
        n_jobs (int, optional): Number of jobs to run in parallel. None means 1. -1 means using all processors. Defaults to -1.
        cv (int, optional): _description_. Defaults to 20.
        verbose (int, optional): Controls the verbosity: the higher, the more messages. Defaults to 1.
        n_iter (int, optional): Number of parameter settings that are sampled. n_iter trades off runtime vs quality of the solution. Defaults to 75.
    """

    # Instantiate the randomized search model
    lr_clf = RandomizedSearchCV(LogisticRegression(random_state=random_state, max_iter=max_iter, class_weight=class_weight),
                                parameters_lr, 
                                scoring= scoring,
                                n_jobs=n_jobs,
                                cv=cv,
                                verbose=verbose, 
                                n_iter=n_iter)

    # Train the Logistic Regression Classifier
    lr_clf = lr_clf.fit(X_train, y_train.values.ravel())

    # Find the best model and print it
    log_reg_model = lr_clf.best_estimator_
    return log_reg_model


def plot_coefficients_from_logistic_model(model, features_names, color='gray', 
                                          negative_regulators = [], 
                                          feature_name_col='feature', condition_name='', top_features=40, 
                                          save_path=None, exp_name='', format='png'):
    
    """Plot the coefficients of the logistic regression model in a barplot, sorted by ascending order.

    Args:
        model (sklearn.linear_model.LogisticRegression): The logistic regression model.
        features_names (list): Names of the features associated with the coefficients. Should be in the same order as in the initial training vector.
        color (str, optional): Color for the barplot. Defaults to 'gray'.
        negative_regulators (list, optional): Feature names to write in bold. Defaults to [].
        feature_name_col (str, optional): Name describing the features type. Defaults to 'feature'.
        condition_name (str, optional): Name of the condition for the model, used for the title. Defaults to ''.
        top_features (int, optional): Number of positive and negative features to plot. Defaults to 40.
        save_path (_type_, optional): Path to save the plot. Defaults to None.
        exp_name (str, optional): Name of the experiment. Defaults to ''.
        format (str, optional): Format type to save the figure. Defaults to 'png'.
    """
    if feature_name_col is None:
        feature_name_col = 'feature'
    coeff = pd.DataFrame(model.coef_[0], columns=['coef'], index=features_names)['coef'].reset_index().rename(columns={'index': feature_name_col}).sort_values('coef')
    
    fig, ax = plt.subplots(1, 1, figsize=(20, 5))

    # Plot the coefficients
    if (len(coeff) > 2*top_features):
        coeff = pd.concat([coeff[:top_features], coeff[-top_features:]])
    sns.barplot(coeff, x=feature_name_col, y='coef', ax=ax, color=color)
    plt.xticks(rotation=90);
    ax.set_title('Coefficients of the logistic regression - '+condition_name, 
                 weight='bold', 
                 fontsize=17);
    ax.set_xlabel(feature_name_col)
    ax.axhline(y=0, color='black')
    ylim = max(-ax.get_ylim()[0], ax.get_ylim()[1])
    plt.ylim(-ylim, ylim)
    
    # Write the negative regulators in bold
    tick_labels = []
    for feature_name_col in list(coeff[feature_name_col].values):
        feature_tick_name = []
        for feature in feature_name_col.split('*'):
            if feature in negative_regulators:
                feature_tick_name.append(r"$\bf{"+feature+"}$")

            else:
                feature_tick_name.append(feature)
        tick_labels.append('*'.join(feature_tick_name))
    plt.tick_params(axis='x', labelcolor="black", labelrotation=0)    
    ax.set_xticks(ticks = np.arange(len(coeff)), 
                labels = tick_labels, rotation = 90)
    
    if save_path is not None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(os.path.join(save_path, exp_name+"_"+condition_name+"_logit_coefficients."+format), bbox_inches='tight', format=format)
    
   


def train_logistic_all_conditions(datasets_per_condition, dataset_key = 'dataset', suffix='', include_plots=False, step_learning=10, feature_name_col='feature',
                                  conditions_names=None, palette=None, negative_regulators=[], class_weight='balanced', threshold=0.5, 
                                  parameters_lr={'solver': ['lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga'], 
                                        'penalty': ['l2', 'l1', 'elasticnet'],
                                        'C': [100, 10, 1.0, 0.1, 0.01]}, min_train_size=10, top_features=40, save_path=None, exp_name='',  format='png', 
                                  scoring='f1', random_state=1, max_iter=1000,  n_jobs=-1,  cv=20,  verbose=1,  n_iter=75):
    
    """Based on a dictionnary containing datasets for multiple treatments/conditions, find the best models through Random Search, train them and
    plot the resulting learning curves and coefficients for each treatment/condition.

    Args:
        datasets_per_condition (dict): dictionnary with conditions/treatments in keys and dictionnaries in values.
        Here is its structure: {treatment1: 
                                        {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test},
                                treatment2:
                                        {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test}}
        dataset_key (str, optional): Key of the treatment dictionnaries containing the dataset informations. Defaults to 'dataset'.
        suffix (str, optional): Suffix that will be added to the created keys in the treatment dictionnaries. Defaults to ''.
        include_plots (bool, optional): Whether to plot the learning curves and coefficients. Defaults to False.
        step_learning (int, optional): Step to plot the learning curves. Defaults to 10.
        feature_name_col (str, optional): Name describing the features type. Defaults to 'feature'.
        conditions_names (list, optional): list of names if you need to use other names than the keys of datasets_per_condition to describe
        the treatments. Defaults to None.
        palette (dict or list, optional): Color palette to use. Can be a list or a dict with the same keys than datasets_per_condition and the colors in values.
        Defaults to None.
        negative_regulators (list, optional): Feature names to write in bold in the coefficients barplot. Defaults to [].
        class_weight (dict or 'balanced'): Weights associated with classes in the form {class_label: weight}. The "balanced" mode
        uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data. Defaults to 'balanced'.
        threshold (float, optional): Decision threshold in the logistic regression. Defaults to 0.5.
        parameters_lr (dict or list of dicts, optional): Dicitonnary with parameters names (str) as keys and distributions or lists of parameters to try. 
        Defaults to {'solver': ['lbfgs', 'newton-cg', 'liblinear', 'sag', 'saga'], 'penalty': ['l2', 'l1', 'elasticnet'], 'C': [100, 10, 1.0, 0.1, 0.01]}.
        min_train_size (int, optional): Minimum training size for the first point to plot in the learning curves. Defaults to 10.
        top_features (int, optional): Number of positive and negative features to plot. Defaults to 40.
        save_path (_type_, optional): Path to save the plot. Defaults to None.
        exp_name (str, optional): Name of the experiment. Defaults to ''.
        format (str, optional): Format type to save the figure. Defaults to 'png'.
        scoring (str, optional): Scoring method. Defaults to 'f1'.
        random_state (int, optional): Used when solver=='sag', 'saga', or 'liblinear' to shuffle the data. Defaults to 1.
        max_iter (int, optional): Maximun number of iterations taken for the solvers to converge. Defaults to 1000.
        n_jobs (int, optional): Number of jobs to run in parallel. None means 1. -1 means using all processors. Defaults to -1.
        cv (int, optional): _description_. Defaults to 20.
        verbose (int, optional): Controls the verbosity: the higher, the more messages. Defaults to 1.
        n_iter (int, optional): Number of parameter settings that are sampled. n_iter trades off runtime vs quality of the solution. Defaults to 75.
    """
    
    if conditions_names is None:
        conditions_names = list(datasets_per_condition.keys())
        
    for i, condition in enumerate(datasets_per_condition.keys()):
        X_train, X_test, y_train, y_test = datasets_per_condition[condition][dataset_key].values()
        if 'class_weight' in datasets_per_condition[condition].keys():
            class_weight = datasets_per_condition[condition]['class_weight']
        else:
            class_weight = class_weight
        logit_model = get_best_logit_model(X_train, y_train, class_weight=class_weight, parameters_lr=parameters_lr, scoring=scoring, random_state=random_state,
                                           max_iter=max_iter, n_jobs=n_jobs, cv=cv, verbose=verbose, n_iter=n_iter)
        datasets_per_condition[condition]['logit_model'+suffix] = logit_model
        y_pred = (logit_model.predict_proba(X_test)[:,1] > threshold).astype(int)
        datasets_per_condition[condition]['logit_scores'+suffix] = classification_report(y_test, y_pred)
        datasets_per_condition[condition]['logit_coef'+suffix] = pd.DataFrame(logit_model.coef_[0], columns=['coef'], index=X_train.columns)
        datasets_per_condition[condition]['chi2'+suffix] =  pd.DataFrame(chi2(pd.concat([X_train, X_test]), 
                                                                            np.hstack((y_train.values.astype(int), 
                                                                                        y_test.values.astype(int)))), 
                                                                        index=['scores', 'pvalues'], 
                                                                        columns=X_train.columns).T
        print("-----------------", conditions_names[i], "-----------------")
        print(datasets_per_condition[condition]['logit_scores'+suffix])
        
    if include_plots:
        plot_curves_all_conditions(datasets_per_condition,
                    dataset_key=dataset_key, 
                    model_key='logit_model'+suffix,
                    step_learning=step_learning,
                    conditions_names=conditions_names,
                    threshold=threshold,
                    min_train_size=min_train_size, 
                    save_path=save_path,
                    exp_name=exp_name,
                    format=format)

        for i, condition in enumerate(datasets_per_condition.keys()):
            if palette is None:
                color = 'gray'
            else:
                color = palette[condition]
            plot_coefficients_from_logistic_model(datasets_per_condition[condition]['logit_model'+suffix], 
                                    color=color, 
                                    negative_regulators=negative_regulators,
                                    feature_name_col=feature_name_col,
                                    condition_name = conditions_names[i],
                                    features_names=datasets_per_condition[condition][dataset_key]['X_train'].columns,
                                    top_features=top_features, save_path=save_path, exp_name=exp_name, format=format)
    


def metrics_versus_thresholds(precision, recall, thresholds_pr, fpr, thresholds_roc, actual_threshold=0.5):
    """Plot the sensitivity, specificity and precision according to the decision threshold.

    Args:
        precision (list or array): precision values as obtained from sklearn.metrics.precision_recall_curve
        recall (list or array): recall values as obtained from sklearn.metrics.precision_recall_curve
        thresholds_pr (list or array): thresholds as obtained from sklearn.metrics.precision_recall_curve
        fpr (list or array): false positive rates as obtained from sklearn.metrics.roc_curve
        thresholds_roc (list or array): thresholds as obtained from sklearn.metrics.roc_curve
        actual_threshold (float, optional): Actual threshold to visualize. Defaults to 0.5.
    """
    
    plt.figure(figsize=(15, 5))
    plt.subplot(121)
    plt.plot(thresholds_pr, recall, label='Sensitivity=TP/(TP+FN)')
    plt.plot(thresholds_pr, precision, label='Precision=TP/(TP+FP)')
    plt.plot(thresholds_roc, 1-fpr, label='Specificity=TN/(TN+FP)')
    plt.xlabel("Threshold")
    plt.axvline(x=actual_threshold, linestyle='dashed', color='black', linewidth=2)
    plt.legend(loc='best')
    plt.title("Classification metrics as a function of the threshold", weight='bold', fontsize=16)
    plt.xlim(min(min(thresholds_roc), min(thresholds_pr)), min(max(thresholds_roc), max(thresholds_pr)))
    
    f1_score = (2*precision*recall)/(precision+recall)
    # f1_score_nt3 = (2*datasets_per_condition['transport_nt3']['PR_curve'][0][:-1]*datasets_per_condition['transport_nt3']['PR_curve'][1][:-1])/(datasets_per_condition['transport_nt3']['PR_curve'][0][:-1]+datasets_per_condition['transport_nt3']['PR_curve'][1][:-1])
    plt.subplot(122)
    plt.plot(thresholds_pr, f1_score, label='F1-score')
    plt.xlabel("Threshold")
    plt.ylabel("F1-score")
    plt.axvline(x=actual_threshold, color='black', linestyle='dashed')
    plt.ylim(0, 1)
    plt.title("F1-score as a function of the threshold", weight='bold', fontsize=16);
    
def metrics_versus_thresholds_all_conditions(datasets_per_condition, model_key, conditions_names=None):
    
    """Plot metrics_versus_thresholds for all conditions/treatments of the dictionnary datasets_per_condition.

    Args:
        datasets_per_condition (dict): dictionnary with conditions/treatments in keys and dictionnaries containing dataset and models in values.
        Here is its structure: {treatment1: 
                                {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test},
                                {model_key: sklearn.linear_model.LogisticRegression model}
                        treatment2:
                                {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test}, 
                                {model_key: sklearn.linear_model.LogisticRegression model}}
        conditions_names (list, optional): short treatment names for titles. Defaults to None, will use the datasets_per_condition keys.
        model_key (str): Key containing the models. 
        
    """
    if conditions_names is None:
        conditions_names = list(datasets_per_condition.keys())
        
    for i, condition in enumerate(list(datasets_per_condition.keys())):
        metrics_versus_thresholds(precision=datasets_per_condition[condition]['PR_curve_'+model_key][0][:-1],
                          recall=datasets_per_condition[condition]['PR_curve_'+model_key][1][:-1],
                          thresholds_pr=datasets_per_condition[condition]['PR_curve_'+model_key][2],
                          fpr=datasets_per_condition[condition]['ROC_curve_'+model_key][0],
                          thresholds_roc=datasets_per_condition[condition]['ROC_curve_'+model_key][2])
        plt.suptitle(conditions_names[i], weight='bold', fontsize=18)
        plt.tight_layout()
        
def get_significant_features_between_models(model1, X_train1, y_train1, model2, X_train2, y_train2, feature_names, n_permutations=2000,
                                            model1_name='model 1', model2_name='model 2'):
    
    """Get statistics on the difference between two models with same features.
    For each feature, a null distribution of the differences is created by retraining both models with shuffled labels, and computing the difference between the
    coefficients of both model. This process is repeated many times (n_permutations). Then, the observed coefficient difference between the two models
    is compared with the null distribution of differences for the given feature. An empirical p-value is computed.

    Args:
        model1 (sklearn.linear_model.LogisticRegression): First model to compare
        X_train1 (array-like): Training vector for the first model
        y_train1 (array-like): Target vector relative to X_train1
        model2 (sklearn.linear_model.LogisticRegression): Second model to compare
        X_train2 (array-like): Training vector for the second model
        y_train2 (array-like): Target vetor relative to X_train2
        feature_names (list): Names of the features.
        n_permutations (int, optional): Number of permutations to generate the null distributions. Defaults to 2000.
        model1_name (str, optional): Name of the first model. Defaults to 'model 1'.
        model2_name (str, optional): Name of the second model. Defaults to 'model 2'.

    Returns:
        pd.DataFrame: a table with statistics for each feature. In rows are the features and here is a description of each column:
            nb_inf (int): Number of samples in the null distribution inferior to the observed value
            nb_sup (int): Numner of samples in the null distribution superior to the observed value
            p_val_empirical_inf (float): Empirical pvalue for one-sided test "less" computed as nb_inf/n_permutations
            p_val_empirical_sup (float): Empirical pvalue for one-sided test "greater" computed as nb_sup/n_permutations
            p_val_empirical (float): Empirical pvalue for two-sided test computed as the minimum between p_val_inf and p_val_sup.
            -log10(p-val): -log10 of the empirical pvalue p_val
            null_mean: mean of the null distribution
            null_std: standard deviation of the null distribution
            z_score: z_score of the observed value calculated from the null distribution statistics
            p_value: p_value derived from z_score
            coef_model1: coefficient value in model1
            coef_model2: coefficient value in model2
            coef_diff: observed coefficient difference

    """
    observed_differences = (model1.coef_ - model2.coef_).ravel()
    observed_differences = pd.DataFrame(observed_differences, columns=['diff_coef'], index=feature_names)
    
    permutations_differences = []
    

    for _ in range(n_permutations):
        model1_perm = clone(model1)
        y_train1_perm = np.random.permutation(y_train1)
        model1_perm.fit(X_train1, y_train1_perm)
        
        model2_perm = clone(model2)
        y_train2_perm = np.random.permutation(y_train2)
        model2_perm.fit(X_train2, y_train2_perm)
        
        diff = (model1_perm.coef_ - model2_perm.coef_).ravel()
        permutations_differences.append(diff)
    perm_distributions = pd.DataFrame(permutations_differences, columns=feature_names)
    
    total_inf = []
    total_sup = []
    for feature in list(feature_names):
        total_inf.append((perm_distributions[feature] < observed_differences.loc[feature, 'diff_coef']).sum())
        total_sup.append((perm_distributions[feature] > observed_differences.loc[feature, 'diff_coef']).sum())
        
    statistics = pd.DataFrame({'nb_inf': total_inf, 'nb_sup': total_sup}, index=feature_names)
    statistics['p_val_empirical_inf'] = statistics['nb_inf']/n_permutations
    statistics['p_val_empirical_sup'] =  statistics['nb_sup']/n_permutations
    statistics['p_val_empirical'] = statistics.apply(lambda x: min(x['p_val_empirical_inf'], x['p_val_empirical_sup']), axis=1)
    statistics['p_val_empirical'].replace(0, 1e-300, inplace=True)
    statistics['-log10(p-val)'] = -np.log10(statistics['p_val_empirical'])
    statistics["null_mean"] = perm_distributions.mean()
    statistics['null_std'] = perm_distributions.std()
    statistics[f"coef_{model1_name}"] = model1.coef_.ravel()
    statistics[f"coef_{model2_name}"] = model2.coef_.ravel()
    statistics["coef_diff"] = observed_differences
    statistics['z_score'] = (statistics['coef_diff'] - statistics['null_mean'])/statistics['null_std']
    statistics['p_value'] = statistics['z_score'].apply(lambda x: norm.sf(abs(x)))


    return statistics




def plot_pvalues(stats, feature_name='feature', p_val_col='p_value', sig_threshold=0.05, keep_all=True):
    """Plot the -log10(p-values).

    Args:
        stats (pd.DataFrame): Table containing the statistics for eacf feature.
        feature_name (str, optional): Name describing the features type. Defaults to 'feature'.
        p_val_col (str, optional): Name of the column in stat containing the p-values. Defaults to 'p_value'.
        sig_threshold (float, optional): Significance threshold to apply for the p-values. Defaults to 0.05.
        keep_all (bool, optional): Whether to keep all the p-values or to plot only the signifivant ones. Defaults to True.
    """
    stats = stats.reset_index().rename(columns={'index': feature_name})
    stats['-log10(p_value)'] = -np.log10(stats[p_val_col])
    if keep_all:
        plt.figure(figsize=(7, len(stats)*0.25))
        sns.barplot(stats.sort_values('-log10(p_value)', ascending=False), x='-log10(p_value)', y=feature_name, color='gray')
        plt.axvline(x=-np.log10(sig_threshold), color='black', linestyle='dashed')
    else:
        plt.figure(figsize=(7, len(stats[stats['p_value'] < sig_threshold])*0.25))
        sns.barplot(stats[stats['p_value'] < sig_threshold].sort_values('-log10(p_value)', ascending=False), x='-log10(p_value)', y=feature_name, color='gray')
    plt.title("Significance of the coefficient differences \n between the two models", weight='bold', fontsize=15)
    
    
    

    
def plot_features_coefficients_per_condition(datasets_per_condition, 
                                            condition1_name, 
                                            condition2_name, 
                                            coef_key='logit_coef', 
                                            model1_name='model 1', 
                                            model2_name='model 2',
                                            palette=None, feature_name='feature', 
                                            scale=False):
    """Plot potitive and negative sorted coefficients for both models, with the corresponding coefficient of the other model.

    Args:
        datasets_per_condition (dict): dictionnary with conditions/treatments in keys and dictionnaries containing dataset and models in values.
        Here is its structure: {treatment1: 
                                {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test},
                                {model_key: sklearn.linear_model.LogisticRegression model}
                        treatment2:
                                {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test}, 
                                {model_key: sklearn.linear_model.LogisticRegression model}}
        condition1_name (_type_): First condition : key in datasets_per_condition 
        condition2_name (_type_): Second condition: key in datasets_per_condition
        coef_key (str, optional): key of the treatment dictionnaries containing the coefficients from the logistic regression. Defaults to 'logit_coef'.
        model1_name (str, optional): Short name for the first condition. Defaults to 'model 1'.
        model2_name (str, optional): Short name for the second condition. Defaults to 'model 2'.
        palette (dict or list, optional): Color palette. Defaults to None.
        feature_name (str, optional): Name describing the features type. Defaults to 'feature'.
        scale (bool, optional): Whether to scale or not the coefficients. Defaults to False.
    """
    coeffs_all = datasets_per_condition[condition1_name][coef_key].rename(columns={'coef': model1_name}).merge(datasets_per_condition[condition2_name][coef_key].rename(columns={'coef': model2_name}),
                                                                                                    left_index=True, 
                                                                                                    right_index=True)
    if scale:
        coeffs_all[model1_name] = coeffs_all[model1_name].apply(lambda x: 
            (abs(x)-abs(coeffs_all[model1_name]).min())/(abs(coeffs_all[model1_name]).max()-abs(coeffs_all[model1_name]).min()))*np.sign(coeffs_all[model1_name])

        coeffs_all[model2_name] = coeffs_all[model2_name].apply(lambda x: 
            (abs(x)-abs(coeffs_all[model2_name]).min())/(abs(coeffs_all[model2_name]).max()-abs(coeffs_all[model2_name]).min()))*np.sign(coeffs_all[model2_name])

    coeffs_melted = pd.melt(coeffs_all, var_name='condition', value_name='coef')
    coeffs_melted[feature_name] = list(coeffs_all.index) * 2

    fig, ax = plt.subplots(1, 2, figsize=(10, 7), sharex=True)
    # Positive regulators model 1
    # RBPs to keep:
    positive_model1 = coeffs_melted[(coeffs_melted['condition'] == model1_name) & (coeffs_melted['coef'] > 0)].sort_values('coef', ascending=False)[feature_name].head(15).values
    coeffs_plot = coeffs_melted[coeffs_melted[feature_name].isin(positive_model1)].sort_values(['condition', 'coef'], ascending=[True, False])
    sns.barplot(coeffs_plot, y=feature_name, x='coef', hue='condition', palette=palette, ax=ax[0])
    ax[0].set_title("Positive regulators in "+model1_name, weight='bold', fontsize=15)
    ax[0].set_ylabel("Relative importance of the coefficient")
    # positive regulators in model 2
    positive_model2 = coeffs_melted[(coeffs_melted['condition'] == model2_name) & (coeffs_melted['coef'] > 0)].sort_values('coef', ascending=False)[feature_name].head(15).values
    coeffs_plot = coeffs_melted[coeffs_melted[feature_name].isin(positive_model2)].sort_values(['condition', 'coef'], ascending=[False, False])
    sns.barplot(coeffs_plot, y=feature_name, x='coef', hue='condition', palette=palette, ax=ax[1])
    ax[1].set_title("Positive regulators in "+model2_name, weight='bold', fontsize=15)
    plt.tight_layout()

    fig, ax = plt.subplots(1, 2, figsize=(10, 7), sharex=True)

    # negative regulators in model 1
    negative_model1 = coeffs_melted[(coeffs_melted['condition'] == model1_name) & (coeffs_melted['coef'] < 0)].sort_values('coef', ascending=True)[feature_name].head(15).values
    coeffs_plot = coeffs_melted[coeffs_melted[feature_name].isin(negative_model1)].sort_values(['condition', 'coef'], ascending=[True, True])
    sns.barplot(coeffs_plot, y=feature_name, x='coef', hue='condition', palette=palette, ax=ax[0])
    ax[0].set_title("Negative regulators in "+model1_name, weight='bold', fontsize=15)

    # negative regulators in model 2
    negative_model2 = coeffs_melted[(coeffs_melted['condition'] == model2_name) & (coeffs_melted['coef'] < 0)].sort_values('coef', ascending=True)[feature_name].head(15).values
    coeffs_plot = coeffs_melted[coeffs_melted[feature_name].isin(negative_model2)].sort_values(['condition', 'coef'], ascending=[False, True])
    sns.barplot(coeffs_plot, y=feature_name, x='coef', hue='condition', palette=palette, ax=ax[1])
    ax[1].set_title("Negative regulators in "+model2_name, weight='bold', fontsize=15)

    plt.tight_layout()

def get_significant_features_one_model(model, 
                                        X_train, 
                                        y_train, 
                                        feature_names, 
                                        n_permutations=2000):
    """To identify the statistically significant predictors of a model, we compared the observed predictor coefficients to their respective null distributions. 
    Therefore, to generate a null distribution for each predictor, we shuffled the response variable label randomly, 
    re-trained the model and re-calculated the predictor coefficients. We repeat this process n_permutations times. 
    To compare the observed predictor coefficients with the null distributions generated, we computed their z-scores as: z = (x-mu)/sigma, 
    where x is the observed predictor coefficient, mu and sigma are the mean and standard deviation of the null distribution for this predictor. 

    Args:
        model (sklearn.linear_model.LogisticRegression): The model to evaluate
        X_train (array-like): Training vector
        y_train (array-like): Target vector relative to X_train.
        feature_names (list): Feature names
        n_permutations (int, optional): Number of permutations to compute the null distributions. Defaults to 2000.

    Returns:
        pd.DataFrame: a table with statistics for each feature. In rows are the features and here is a description of each column:
            nb_inf (int): Number of samples in the null distribution inferior to the observed value
            nb_sup (int): Numner of samples in the null distribution superior to the observed value
            p_val_empirical_inf (float): Empirical pvalue for one-sided test "less" computed as nb_inf/n_permutations
            p_val_empirical_sup (float): Empirical pvalue for one-sided test "greater" computed as nb_sup/n_permutations
            p_val_empirical (float): Empirical pvalue for two-sided test computed as the minimum between p_val_inf and p_val_sup.
            -log10(p-val): -log10 of the empirical pvalue p_val
            coef: coefficient value in model
            null_mean: mean of the null distribution
            null_std: standard deviation of the null distribution
            z_score: z_score of the observed value calculated from the null distribution statistics
            p_value: p_value derived from z_score


    """
    observed = (model.coef_).ravel()
    observed = pd.DataFrame(observed, columns=['coef'], index=feature_names)

    permutations = []
    

    for _ in range(n_permutations):
        model_perm = clone(model)
        y_train_perm = np.random.permutation(y_train)
        model_perm.fit(X_train, y_train_perm)
        
        permutations.append(model_perm.coef_.ravel())
    perm_distributions = pd.DataFrame(permutations, columns=feature_names)
    
    total_inf = []
    total_sup = []
    for feature in list(feature_names):
        total_inf.append((perm_distributions[feature] < observed.loc[feature, 'coef']).sum())
        total_sup.append((perm_distributions[feature] > observed.loc[feature, 'coef']).sum())
        
    statistics = pd.DataFrame({'nb_inf': total_inf, 'nb_sup': total_sup}, index=feature_names)
    statistics['p_val_empirical_inf'] = statistics['nb_inf']/n_permutations
    statistics['p_val_empirical_sup'] =  statistics['nb_sup']/n_permutations
    statistics['p_val_empirical'] = statistics.apply(lambda x: min(x['p_val_empirical_inf'], x['p_val_empirical_sup']), axis=1)
    statistics['p_val_empirical'].replace(0, 1e-300, inplace=True)
    statistics['-log10(p-val)'] = -np.log10(statistics['p_val_empirical'])
    statistics["coef"] = model.coef_.ravel()
    statistics["null_mean"] = perm_distributions.mean()
    statistics['null_std'] = perm_distributions.std()
    statistics['z_score'] = (statistics['coef'] - statistics['null_mean'])/statistics['null_std']
    statistics['p_value'] = statistics['z_score'].apply(lambda x: norm.sf(abs(x)))

    return statistics
        

def extract_common_regulators(stats_dict, sig_threshold=0.05):
    """From a dictionnary containing multiple statistics table for features from different models, extract the common significant predictors.

    Args:
        stats_dict (_type_): _description_
        sig_threshold (float, optional): _description_. Defaults to 0.05.

    Returns:
        pd.DataFrame, pd.DataFrame: Statistics table for the common positive regulators, Statistics table for the common negative regulators.
    """
    for i, condition in enumerate(stats_dict.keys()):
        significant_regulators = stats_dict[condition][stats_dict[condition]['p_value'] < sig_threshold]
        significant_regulators.columns = [col+'_'+condition for col in significant_regulators.columns]
        if i == 0:
            common_regulators = significant_regulators.copy()
        else:
            common_regulators = common_regulators.merge(significant_regulators,
                                                             left_index=True, 
                                                             right_index=True)
            
    z_scores_columns = [col for col in common_regulators.columns if col.startswith('z_score')]
    
    common_positive_regulators = common_regulators.copy()
    common_negative_regulators = common_regulators.copy()
    for col in z_scores_columns:
        common_positive_regulators = common_positive_regulators[common_positive_regulators[col] > 0]
        common_negative_regulators = common_negative_regulators[common_negative_regulators[col] < 0]
        
    return common_positive_regulators, common_negative_regulators

def plot_zscores(stats, feature_name='feature', ascending=True, palette=None):
    """From a statistics table, plot the z_scores colored by condition.

    Args:
        stats (pd.DataFrame): a statistics table
        feature_name (str, optional): Name describing the features type. Defaults to 'feature'.
        ascending (bool, optional): Whether to sort in ascending order. Defaults to True.
        palette (dict or list, optional): color palette. Defaults to None.
    """
    z_scores_columns = [col for col in stats.columns if col.startswith('z_score')]
    features = list(stats.index)
    
    if len(z_scores_columns) > 1:
        df_plot = pd.melt(stats[z_scores_columns], var_name='condition', value_name='z_score')
        df_plot[feature_name] = features*len(z_scores_columns)
        df_plot['condition'] = df_plot['condition'].apply(lambda x: x.split('z_score_')[1])

    else:
        df_plot = stats.reset_index().rename(columns={'index': feature_name})
    
    sns.barplot(data=df_plot.sort_values('z_score', ascending=ascending), 
                x='z_score', 
                y=feature_name, 
                hue='condition', 
                palette=palette)    

# def plot_zscores(stats, feature_name='feature', ascending=True, palette=None, sig_z_threshold=None):
#     z_scores_columns = [col for col in stats.columns if col.startswith('z_score')]
#     features = list(stats.index)
#     plt.figure(figsize=(5, 0.25*len(stats)))
    
#     if len(z_scores_columns) > 1:
#         df_plot = pd.melt(stats[z_scores_columns], var_name='condition', value_name='z_score')
#         df_plot[feature_name] = features*len(z_scores_columns)
#         df_plot['condition'] = df_plot['condition'].apply(lambda x: x.split('z_score_')[1])

#     else:
#         df_plot = stats.reset_index().rename(columns={'index': feature_name})
#         df_plot['condition'] = ''
    
#     sns.barplot(data=df_plot.sort_values('z_score', ascending=ascending), 
#                 x='z_score', 
#                 y=feature_name, 
#                 hue='condition', 
#                 palette=palette)
    
#     if sig_z_threshold is not None:
#         plt.axvline(x=sig_z_threshold, linestyle='dashed', color='black')
#         plt.axvline(x=-sig_z_threshold, linestyle='dashed', color='black')

    
def extract_not_common_regulators(stats1, stats2, sig_threshold=0.05):
    """From two statistics tables corresponding to two different models with the same features, extract the specific regulators to each model.
    Args:
        stats1 (pd.DataFrame): Statistics table for features in model 1 
        stats2 (pd.DataFrame): Statistics table for features in model 2
        sig_threshold (float, optional): Significance threshold (in terms of p-values). Defaults to 0.05.

    Returns:
        pd.DataFrame, pd.DataFrame: Statistics table with the features being predictors specific to model 1, 
                                    Statistics table with the features being predictors specific to model 2
    """
    significant_regulators1 = list(stats1[stats1['p_value'] < sig_threshold].index)
    significant_regulators2 = list(stats2[stats2['p_value'] < sig_threshold].index)
    
    regulators1_only = [item for item in significant_regulators1 if item not in significant_regulators2]
    regulators2_only = [item for item in significant_regulators2 if item not in significant_regulators1]
    
    return stats1.loc[regulators1_only], stats2.loc[regulators2_only]


def merge_models_statistics(stats_dict):
    """Merge statistics summary tables from a dictionnary of multiple tables.

    Args:
        stats_dict (dict): Dictionnary with conditions/treatments in keys and statistics tables for features in values.

    Returns:
        pd.DataFrame: The merged table.
    """
    for i, condition in enumerate(stats_dict.keys()):
        regulators = stats_dict[condition].copy()
        regulators.columns = [col+'_'+condition for col in regulators.columns]
        if i == 0:
            all_regulators = regulators.copy()
        else:
            all_regulators = all_regulators.merge(regulators,
                                                             left_index=True, 
                                                             right_index=True)
    return all_regulators


def plot_common_and_specific_regulators_zscores(datasets_per_condition,
                                                feature_name='feature',
                                         dataset_key='dataset',
                                         model_key='logit_model',
                                         conditions_names=None, 
                                         model1_name=None, 
                                         model2_name=None,
                                         n_permutations=2000,
                                         sig_threshold=0.05,
                                         save_path=None,
                                         format='png', 
                                         palette_conditions=None,
                                         transparent=True, 
                                         return_groups=False):
    """

    Args:
        datasets_per_condition (dict): dictionnary with conditions/treatments in keys and dictionnaries containing dataset and models in values.
        Here is its structure: {treatment1: 
                                {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test},
                                {model_key: sklearn.linear_model.LogisticRegression model}
                        treatment2:
                                {dataset_key: {'X_train': X_train, 'X_test': X_test, 'y_train': y_train, 'y_test': y_test}, 
                                {model_key: sklearn.linear_model.LogisticRegression model}}
        conditions_names (list, optional): short treatment names for titles. Defaults to None, will use the datasets_per_condition keys.
        model1_name (str, optional): Name of the first model to compare with model2 to extract specific regulators.
        Has to be an item contained in conditions_names or in keys of datasets_per_condition if conditions_names is No. Defaults to None.
        model2_name (str, optional): Name of the first model to compare with model1 to extract specific regulators.
        Has to be an item contained in consitions_names or in keys of datasets_per_condition if conditions_names is None. Defaults to None.
        n_permutations (int, optional): Number of permutations to perform to create the null distributions. Default to 2000.
        sig_threshold (float, optional): Significance threshold.
    """
    if conditions_names is None:
        conditions_names = list(datasets_per_condition.keys())

    statistics = {}

    for condition, condition_name in zip(datasets_per_condition.keys(), conditions_names):
        X_train, X_test, y_train, y_test = datasets_per_condition[condition][dataset_key].values()
        model = datasets_per_condition[condition][model_key]
        stats = get_significant_features_one_model(model=model, 
                                                X_train=X_train,
                                                y_train=y_train,
                                                feature_names=X_train.columns,
                                                n_permutations=n_permutations)

        statistics[condition_name] = stats

    merged_stats = merge_models_statistics(statistics)
    model1_specific, model2_specific = extract_not_common_regulators(stats1=statistics[model1_name], stats2=statistics[model2_name], sig_threshold=sig_threshold)
    pos_common, neg_common = extract_common_regulators(statistics, sig_threshold=sig_threshold)

    plt.figure(figsize=(14, 21))
    plt.subplot(321)
    if len(pos_common) > 0:
        plot_zscores(pos_common, feature_name=feature_name, ascending=False, palette=palette_conditions)
    plt.title("Common positive regulators", weight='bold', fontsize=17)
    plt.ylabel(None)
    plt.subplot(322)
    if len(neg_common) > 0:
        plot_zscores(neg_common, feature_name=feature_name, ascending=True, palette=palette_conditions)
    plt.title("Common negative regulators", weight='bold', fontsize=17)
    plt.ylabel(None)
    plt.subplot(323)
    if len(merged_stats.loc[list(model1_specific[model1_specific['z_score'] > 0].index)]):
        plot_zscores(merged_stats.loc[list(model1_specific[model1_specific['z_score'] > 0].index)], feature_name=feature_name, ascending=False, palette=palette_conditions)
    plt.title(model1_name+"-specific positive regulators", weight='bold', fontsize=17)
    plt.ylabel(None)
    plt.subplot(324)
    if len(merged_stats.loc[list(model1_specific[model1_specific['z_score'] < 0].index)]):
        plot_zscores(merged_stats.loc[list(model1_specific[model1_specific['z_score'] < 0].index)], feature_name=feature_name, ascending=True, palette=palette_conditions)
    plt.title(model1_name+"-specific negative regulators", weight='bold', fontsize=17)
    plt.ylabel(None)
    plt.subplot(325)
    if len(merged_stats.loc[list(model2_specific[model2_specific['z_score'] > 0].index)]):
        plot_zscores(merged_stats.loc[list(model2_specific[model2_specific['z_score'] > 0].index)], feature_name=feature_name, ascending=False, palette=palette_conditions)
    plt.title(model2_name+"-specific positive regulators", weight='bold', fontsize=17)
    plt.ylabel(None)
    plt.subplot(326)
    if len(merged_stats.loc[list(model2_specific[model2_specific['z_score'] < 0].index)]):
        plot_zscores(merged_stats.loc[list(model2_specific[model2_specific['z_score'] < 0].index)], feature_name=feature_name, ascending=True, palette=palette_conditions)
    plt.title(model2_name+"-specific negative regulators", weight='bold', fontsize=17)
    plt.ylabel(None)
    plt.tight_layout()
    
    if return_groups:
        dic_summary = {}
        dic_summary['Positive common'] = pos_common
        dic_summary['Negative common'] = neg_common
        dic_summary[f"Positive {model1_name}-specific"] = merged_stats.loc[list(model1_specific[model1_specific['z_score'] > 0].index)]
        dic_summary[f"Negative {model1_name}-specific"] = merged_stats.loc[list(model1_specific[model1_specific['z_score'] < 0].index)]
        dic_summary[f"Positive {model2_name}-specific"] = merged_stats.loc[list(model2_specific[model2_specific['z_score'] > 0].index)]
        dic_summary[f"Negative {model2_name}-specific"] = merged_stats.loc[list(model2_specific[model2_specific['z_score'] < 0].index)]
        return dic_summary
    
    if save_path is not None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(os.path.join(save_path, "significant_regulators_zscores."+format), 
                    format=format, bbox_inches='tight', transparent=transparent)      
        
        
# def plot_stats_difference_between_models(stats, model1_name='model 1', model2_name='model 2', palette=None, sig_threshold=0.05, feature_name = 'feature'):
#     # negative_regulators_model2 = stats[(stats['coef_'+model2_name] < 0) & (stats['p_val_sup'] < sig_threshold) & (abs(stats['coef_'+model2_name]) > abs(stats['coef_'+model1_name]))]
#     # negative_regulators_model1 = stats[(stats['coef_'+model1_name] < 0) & (stats['p_val_inf'] < sig_threshold) & (abs(stats['coef_'+model2_name]) < abs(stats['coef_'+model1_name]))]
#     # positive_regulators_model2 = stats[(stats['coef_'+model2_name] > 0) & (stats['p_val_inf'] < sig_threshold) & (abs(stats['coef_'+model2_name]) > abs(stats['coef_'+model1_name]))]
#     # positive_regulators_model1 = stats[(stats['coef_'+model1_name] > 0) & (stats['p_val_sup'] < sig_threshold) & (abs(stats['coef_'+model2_name]) < abs(stats['coef_'+model1_name]))]
#     negative_regulators_model2 = stats[(stats['p_val_sup'] < sig_threshold) & (abs(stats['coef_'+model2_name]) > abs(stats['coef_'+model1_name]))]
#     negative_regulators_model1 = stats[(stats['p_val_inf'] < sig_threshold) & (abs(stats['coef_'+model2_name]) < abs(stats['coef_'+model1_name]))]
#     positive_regulators_model2 = stats[(stats['p_val_inf'] < sig_threshold) & (abs(stats['coef_'+model2_name]) > abs(stats['coef_'+model1_name]))]
#     positive_regulators_model1 = stats[(stats['p_val_sup'] < sig_threshold) & (abs(stats['coef_'+model2_name]) < abs(stats['coef_'+model1_name]))]
#     fig, ax = plt.subplots(2, 2, figsize=(10, 10))
#     ## Negative regulators model 1
#     if len(negative_regulators_model1) > 0:
#         negative_regulators_model1.sort_values('coef_'+model1_name, inplace=True)
#         features = list(negative_regulators_model1.index)
#         df_plot = pd.melt(negative_regulators_model1[['coef_'+model1_name, 'coef_'+model2_name]], value_name='coef', var_name='condition')
#         df_plot['condition'] = df_plot['condition'].apply(lambda x: x.split('coef_')[1])
#         df_plot[feature_name] = features * 2
#         sns.barplot(data=df_plot, x='coef', y=feature_name, ax=ax[0][0], hue='condition', palette=palette)
#         ax[0][0].set_title("Negative regulators "+model1_name, weight='bold')
#         ax[0][0].axvline(x=0, color='black', linestyle='dashed')
#     else:
#         ax[0][0].set_xlim(-0.1, 0.1)
#     # Negative regulators model 2 
#     if len(negative_regulators_model2) > 0:
#         negative_regulators_model2.sort_values('coef_'+model2_name, inplace=True)
#         features = list(negative_regulators_model2.index)
#         df_plot = pd.melt(negative_regulators_model2[['coef_'+model1_name, 'coef_'+model2_name]], value_name='coef', var_name='condition')
#         df_plot['condition'] = df_plot['condition'].apply(lambda x: x.split('coef_')[1])
#         df_plot[feature_name] = features * 2
#         sns.barplot(data=df_plot, x='coef', y=feature_name, ax=ax[0][1], hue='condition', palette=palette)
#         ax[0][1].set_title("Negative regulators "+model2_name, weight='bold')  
#         ax[0][1].axvline(x=0, color='black', linestyle='dashed') 
#     else:
#         ax[0][1].set_xlim(-0.1, 0.1)
#     # Positive regulators model 1
#     if len(positive_regulators_model1) > 0:
#         positive_regulators_model1.sort_values('coef_'+model1_name, ascending=False, inplace=True)
#         features = list(positive_regulators_model1.index)
#         df_plot = pd.melt(positive_regulators_model1[['coef_'+model1_name, 'coef_'+model2_name]], value_name='coef', var_name='condition')
#         df_plot['condition'] = df_plot['condition'].apply(lambda x: x.split('coef_')[1])
#         df_plot[feature_name] = features * 2
#         sns.barplot(data=df_plot, x='coef', y=feature_name, ax=ax[1][0], hue='condition', palette=palette)
#         ax[1][0].set_title("Positive regulators "+model1_name, weight='bold')  
#         ax[1][0].axvline(x=0, color='black', linestyle='dashed')
#     else:
#         ax[1][0].set_xlim(-0.1, 0.1)
#     # Positive regulators model 2 
#     if len(positive_regulators_model2) > 0:
#         positive_regulators_model2.sort_values('coef_'+model2_name, ascending=False, inplace=True)
#         features = list(positive_regulators_model2.index)
#         df_plot = pd.melt(positive_regulators_model2[['coef_'+model1_name, 'coef_'+model2_name]], value_name='coef', var_name='condition')
#         df_plot['condition'] = df_plot['condition'].apply(lambda x: x.split('coef_')[1])
#         df_plot[feature_name] = features * 2
#         sns.barplot(data=df_plot, x='coef', y=feature_name, ax=ax[1][1], hue='condition', palette=palette)
#         ax[1][1].set_title("Positive regulators "+model2_name, weight='bold')
#         ax[1][1].axvline(x=0, color='black', linestyle='dashed')
#     else:
#         ax[1][1].set_xlim(-0.1, 0.1)
#     ax[0][0].set_xlim(min(ax[0][0].get_xlim()[0], ax[0][1].get_xlim()[0]), max(ax[0][0].get_xlim()[1], ax[0][1].get_xlim()[1]))
#     ax[0][1].set_xlim(min(ax[0][0].get_xlim()[0], ax[0][1].get_xlim()[0]), max(ax[0][0].get_xlim()[1], ax[0][1].get_xlim()[1]))
#     ax[1][0].set_xlim(min(ax[1][0].get_xlim()[0], ax[1][1].get_xlim()[0]), max(ax[1][0].get_xlim()[1], ax[1][1].get_xlim()[1]))
#     ax[1][1].set_xlim(min(ax[1][0].get_xlim()[0], ax[1][1].get_xlim()[0]), max(ax[1][0].get_xlim()[1], ax[1][1].get_xlim()[1]))
#     plt.suptitle("Direction of the significant differences in coefficients between the two models", weight='bold', fontsize=15)
#     plt.tight_layout()


# def plot_summary_coefficients_differences(datasets_per_condition, condition1_name, condition2_name, model_key='logit_model', 
#                                           dataset_key='dataset',
#                                           model1_name = 'model 1', model2_name = 'model 2',
#                                           n_permutations=2000, sig_threshold=0.05, palette=None, feature_name='feature', keep_all=True):
    
#     model1 = datasets_per_condition[condition1_name][model_key]
#     X_train1=datasets_per_condition[condition1_name][dataset_key]['X_train']
#     y_train1=datasets_per_condition[condition1_name][dataset_key]['y_train']
    
#     model2 = datasets_per_condition[condition2_name][model_key]
#     X_train2=datasets_per_condition[condition2_name][dataset_key]['X_train']
#     y_train2=datasets_per_condition[condition2_name][dataset_key]['y_train']
    
#     feature_names = X_train1.columns
    
#     stats = get_significant_features_between_models(model1=model1,
#                                                     X_train1=X_train1, 
#                                                     y_train1=y_train1, 
#                                                     model2=model2,
#                                                     X_train2=X_train2,
#                                                     y_train2=y_train2,
#                                                     feature_names=feature_names,
#                                                     model1_name=model1_name,
#                                                     model2_name=model2_name,
#                                                     n_permutations=n_permutations)
    
    
    
#     plot_pvalues(stats, feature_name=feature_name, sig_threshold=sig_threshold, keep_all=keep_all)
    
#     plot_stats_difference_between_models(stats=stats, 
#                                          model1_name=model1_name, 
#                                          model2_name=model2_name, 
#                                          palette=palette, 
#                                          feature_name=feature_name, 
#                                          sig_threshold=sig_threshold)
    
#     return stats