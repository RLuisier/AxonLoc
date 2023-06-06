import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve, precision_recall_fscore_support
from sklearn.metrics import mean_squared_error, r2_score, accuracy_score, f1_score
import matplotlib.pyplot as plt
import os 
import shap

def plot_learning_curves(model, X_train, X_test, 
                         y_train, y_test, step=1, regressor=False, 
                         metric='accuracy', threshold=0.5, min_train_size=10):
    """Plot learning curves for training and testing set.

    Args:
        model (sklearn model): The model to evaluate
        X_train (array-like): Training vetor
        X_test (array-like): Testing vector
        y_train (array-like): Target vector relative to X_train.
        y_test (array-like): Target vector relative to X_test.
        step (int, optional): Step to compute accuracy along training set size. Defaults to 1.
        min_train_size (int, optional): Minimum training size for the first point of the learning curves.
    """
    train_accuracy, test_accuracy = [], []
    train_f1, test_f1 = [], []
    for m in np.arange(min_train_size, len(X_train), step=step):
        model.fit(X_train[:m], y_train[:m])
        y_train_predict = (model.predict_proba(X_train[:m])[:,1] > threshold).astype(int)
        y_test_predict = (model.predict_proba(X_test)[:,1] > threshold).astype(int)
        if regressor:
            train_accuracy.append(r2_score(y_train[:m], y_train_predict))
            test_accuracy.append(r2_score(y_test, y_test_predict))
        else:
            train_accuracy.append(accuracy_score(y_train[:m], y_train_predict))
            test_accuracy.append(accuracy_score(y_test, y_test_predict))
            train_f1.append(f1_score(y_train[:m], y_train_predict))
            test_f1.append(f1_score(y_test, y_test_predict))

        # train_accuracy.append(sum(y_train_predict == y_train[:m])/len(y_train[:m]))
        # test_accuracy.append(sum(y_test_predict == y_test)/len(y_test))
    
    if metric == 'f1-score':
        plt.plot(np.arange(10, len(X_train), step=step), train_f1, "r-+", linewidth=2, label="train")
        plt.plot(np.arange(10, len(X_train), step=step), test_f1, "b-", linewidth=3, label="test")
        plt.xlabel("Training set size", fontsize=18)
        plt.ylabel("f1-score (positive class)", fontsize=18)
        plt.legend(loc="upper right", fontsize=14)
        if not regressor:
            plt.ylim(0,1)
        plt.grid()
        plt.title('Learning curve')
    else:   
        plt.plot(np.arange(10, len(X_train), step=step),train_accuracy, "r-+", linewidth=2, label="train")
        plt.plot(np.arange(10, len(X_train), step=step),test_accuracy, "b-", linewidth=3, label="test")
        
        if type(y_train.sum()) == pd.Series:
            y_train_sum = y_train.sum().values[0]
        else:
            y_train_sum = y_train.sum()
        plt.axhline(y=max(y_train_sum/(len(y_train)), (len(y_train)-y_train_sum)/(len(y_train))), linestyle='dashed', color='black')
        plt.xlabel("Training set size", fontsize=18)
        plt.ylabel("accuracy", fontsize=18)
        plt.legend(loc="upper right", fontsize=14)
        if not regressor:
            plt.ylim(0,1)
        plt.grid()
        plt.title('Learning curve')



def plot_roc_curves(model, X_train, X_test, y_train, y_test, return_values=False, threshold=0.5):
    """Plot ROC curves for training and testing set.

    Args:
        model (sklearn model): The model to evaluate
        X_train (array-like): Training vetor
        X_test (array-like): Testing vector
        y_train (array-like): Target vector relative to X_train.
        y_test (array-like): Target vector relative to X_test.
    """

    y_scores_nn=model.predict_proba(X_train)[:,1]    
    y_scores_nn_test=model.predict_proba(X_test)[:,1] 
    auc_train=roc_auc_score(y_train, y_scores_nn)
    auc_test=roc_auc_score(y_test, y_scores_nn_test)

    fpr_nn, tpr_nn, thresholds = roc_curve(y_train,y_scores_nn)
    fpr_nn_test, tpr_nn_test, thresholds_test = roc_curve(y_test,y_scores_nn_test)

    t_train = min(thresholds, key=lambda x:abs(x-threshold))
    t_train_idx = list(thresholds).index(t_train)
    t_test = min(thresholds_test, key=lambda x:abs(x-threshold))
    t_test_idx = list(thresholds_test).index(t_test)
    plt.plot(fpr_nn, tpr_nn, "r", linewidth=2, label="train (AUC="+"%0.2f" % auc_train+")", zorder=1)
    plt.plot(fpr_nn_test, tpr_nn_test, "b", linewidth=2, label="test (AUC="+"%0.2f" % auc_test+")", zorder=1)
    plt.scatter([fpr_nn[t_train_idx]], [tpr_nn[t_train_idx]], s=65, marker='+', color='black', linewidths=4, label='threshold=0.5', zorder=2)
    plt.scatter([fpr_nn_test[t_test_idx]], [tpr_nn_test[t_test_idx]], s=65, marker='+', color='black', linewidths=4, zorder=2)
    
    plt.plot([0, 1], [0, 1], 'k--') # dashed diagonal
    plt.xlabel('FPR',fontsize=16)
    plt.ylabel('TPR',fontsize=16)
    plt.grid()
    plt.legend(loc="best")
    
    if return_values:
        return fpr_nn, tpr_nn, thresholds
    

def plot_precision_recall_curves(model, X_train, X_test, y_train, y_test, 
                                 threshold=0.5, return_values=False):
    """Plot ROC curves for training and testing set.

    Args:
        model (sklearn model): The model to evaluate
        X_train (array-like): Training vetor
        X_test (array-like): Testing vector
        y_train (array-like): Target vector relative to X_train.
        y_test (array-like): Target vector relative to X_test.
    """

    y_scores_nn=model.predict_proba(X_train)[:,1]    
    y_scores_nn_test=model.predict_proba(X_test)[:,1] 

    precision, recall, thresholds = precision_recall_curve(y_train,y_scores_nn, pos_label=1)
    precision_test, recall_test, thresholds_test = precision_recall_curve(y_test,y_scores_nn_test, pos_label=1)

    # Get the threshold closest to the one we use
    t_train = min(thresholds, key=lambda x:abs(x-threshold))
    t_train_idx = list(thresholds).index(t_train)
    t_test = min(thresholds_test, key=lambda x:abs(x-threshold))
    t_test_idx = list(thresholds_test).index(t_test)
    plt.plot(recall, precision, "r", linewidth=2, label="train", zorder=1)
    plt.plot(recall_test, precision_test, "b", linewidth=2, label="test", zorder=1)
    plt.scatter([recall[t_train_idx]], [precision[t_train_idx]], s=65, marker='+', color='black', linewidths=4, zorder=2)
    plt.scatter([recall_test[t_test_idx]], [precision_test[t_test_idx]], s=65, marker='+', color='black', linewidths=4, label='threshold=0.5', zorder=2)
    plt.plot([0, 1], [0.5, 0.5], 'k--') # dashed diagonal
    plt.xlabel('Recall (Positive class)',fontsize=16)
    plt.ylabel('Precision (Positive class)',fontsize=16)
    plt.grid()
    plt.legend(loc="best")
    plt.ylim(-0.1,1.1)
    plt.xlim(-0.1,1.1)
    
    if return_values:
        return precision, recall, thresholds



def plot_curves_all_conditions(datasets_per_condition, 
                               dataset_key = 'oversampled_dataset', 
                               model_key = 'MLP_model',
                               conditions_names=None, 
                               step_learning=100, 
                               save_path=None,
                               exp_name='',
                               return_values=False,
                               threshold=0.5,
                               min_train_size=10, 
                               format='png'):
    """For a given dictionary with multiple conditions in keys, plot the learning and ROC
    curves subsequentially.

    Args:
        datasets_per_condition (dict): dictionnary with conditions/treatments in keys and dictionnaries in values.
        These dictionnaries contains dataset_key and model_key in keys, containing the dataset and model to evaluate respectively.
        dataset_key (str, optional): Key of the subdictionnaries containing the dataset to evaluate. Defaults to 'oversampled_dataset'.
        model_key (str, optional): Key of the subdictionnaries containing the model to evaluate. Defaults to 'MLP_model'.
        conditions_names (list, optional): Names of the treatment as you want them to be written in plot titles. Defaults to None, will use the dictionnary keys.
        step_learning (int, optional): The step length you want to use to compute learning curves. Defaults to 100.
    """
    
    nb_conditions = len(datasets_per_condition.keys())
    if conditions_names is None:
        conditions_names = list(datasets_per_condition.keys())

    plt.figure(figsize=(28, 7*nb_conditions))

    for i, condition in enumerate(datasets_per_condition.keys()):

        X_train, X_test, y_train, y_test = datasets_per_condition[condition][dataset_key].values()
        # Plot learning curves
        plt.subplot(nb_conditions, 4, (i*4)+1)
        plot_learning_curves(model=datasets_per_condition[condition][model_key], 
                            X_train=X_train, 
                            X_test=X_test, 
                            y_train=y_train, 
                            y_test=y_test, 
                            step=step_learning,
                            threshold=threshold,
                            min_train_size=min_train_size)
        datasets_per_condition[condition][model_key].fit(X_train, y_train)
        final_accuracy = datasets_per_condition[condition][model_key].score(X_test, y_test)
        plt.text(s="Final accuracy: {:.2f}".format(final_accuracy), y=0.2, x=len(X_train)/2, fontsize=14, weight='bold')
        plt.title("Learning curves - "+conditions_names[i], weight='bold', fontsize=17)
        
        
        # Plot learning curves
        plt.subplot(nb_conditions, 4, (i*4)+2)
        plot_learning_curves(model=datasets_per_condition[condition][model_key], 
                            X_train=X_train, 
                            X_test=X_test, 
                            y_train=y_train, 
                            y_test=y_test, 
                            step=step_learning,
                            metric='f1-score',
                            threshold=threshold)
        datasets_per_condition[condition][model_key].fit(X_train, y_train)
        y_pred = (datasets_per_condition[condition][model_key].predict_proba(X_test)[:,1] > threshold).astype(int)
        precision, recall, final_f1, _ = precision_recall_fscore_support(y_test, y_pred, pos_label=1, average='binary')
        plt.text(s="Final f1_score: {:.2f}".format(final_f1), y=0.2, x=len(X_train)/2, fontsize=14, weight='bold')
        plt.title("Learning curves - "+conditions_names[i], weight='bold', fontsize=17)
        
        # Plot ROC curves
        plt.subplot(nb_conditions, 4, (i*4)+3)
        datasets_per_condition[condition]['ROC_curve_'+model_key] = plot_roc_curves(model=datasets_per_condition[condition][model_key], 
                            X_train=X_train, 
                            X_test=X_test, 
                            y_train=y_train, 
                            y_test=y_test, return_values=True)

        plt.title('ROC curves - '+conditions_names[i], weight='bold', fontsize=17)
        
        # Plot the precision-recall curves
        plt.subplot(nb_conditions, 4, (i*4)+4)
        datasets_per_condition[condition]['PR_curve_'+model_key] = plot_precision_recall_curves(model=datasets_per_condition[condition][model_key], 
                            X_train=X_train, 
                            X_test=X_test, 
                            y_train=y_train, 
                            y_test=y_test, return_values=True)
        
                
        plt.text(s="Recall: {:.2f}".format(recall), y=0.2, x=0.4, fontsize=14, weight='bold')
        plt.text(s="Precision: {:.2f}".format(precision), y=0.1, x=0.4, fontsize=14, weight='bold')

        plt.title('Precision recall curves - '+conditions_names[i], weight='bold', fontsize=17)
              
    if save_path is not None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(os.path.join(save_path, exp_name+"_learning_ROC_PR_curves."+format), bbox_inches='tight', format=format)


def get_outliers(serie):
    """Get outliers of a ps.Series.
    Args:
        serie (pd.Series): pd.Serie with numerical values

    Returns:
        ps.Series: Two series with the lower outliers and upper outliers.
    """
    Q1, Q3 = serie.quantile([0.25, 0.75]).values
    lower_bound = Q1 - 1.5*(Q3-Q1)
    upper_bound = Q3 + 1.5*(Q3-Q1)
    return serie[serie < lower_bound], serie[serie > upper_bound]

def set_binary_labels(df_data, col_scores, threshold_up = 1, threshold_down=-1, label_name='label'):
    """From a dataset with a continuous column, create a binary label column
    and return the corresponding dataset based on the thresholds. If the threshold_up and threshold_down are not
    equal, the sampled in between will be excluded.

    Args:
        df_data (pd.DataFrame): dataset containing the continuous column to build the label on
        col_scores (str): name of the column with continuous values
        threshold_up (int, optional): Threshold above which the samples will be labeled as belonging to the positive class. Defaults to 1.
        threshold_down (int, optional): Threshold under which the samples will be labeled as belonging to the negative class. Defaults to -1.
        label_name (str, optional): Name of the column with the binary label. Defaults to 'label'.

    Returns:
        pd.DataFrame: The dataset filtered and labeled
    """
    interesting_data = df_data[(df_data[col_scores] > threshold_up) | (df_data[col_scores] < threshold_down)]
    interesting_data[label_name] = interesting_data[col_scores].apply(lambda x: x > threshold_up).astype(int)
    return interesting_data


def get_interactions_dataset_from_combs(original_X, list_interactions, split_symbol='*'):
    """Create interactions terms from a list of combinations (list of tuples or list of strings) 
    and return the new dataset.

    Args:
        original_X (pd.DataFrame): The original dataset.
        list_interactions (list): List of tuples representing the interactions to create, 
        or list of strings such as ["RBP1*RBP2", "RBP2*RBP4"]. 
        split_symbol (str, optional): symbol to split and retrieve the two interacting terms in 
        the strings. Defaults to '*'.
    Returns:
        pd.DataFrame: The dataset with the created interaction terms.
    """
    interaction_X = original_X.copy()
    
    for interaction in list_interactions:
        if type(interaction) == str:
            features = interaction.split(split_symbol)
            interaction_X[interaction.replace(split_symbol, '*')] = interaction_X.loc[:, features].prod(axis=1)
        else:
            feature1 = interaction[0]
            feature2 = interaction[1]
            interaction_X[f"{feature1}*{feature2}"] = interaction_X[feature1]*interaction_X[feature2]
    return interaction_X



def create_interaction_dataset_from_list(X_train, 
                               X_test, 
                               y_train, 
                               y_test, 
                               list_interactions, 
                               split_symbol='*',
                               keep_old=True):
    """From a dataset, create a new one with interactions terms.

    Args:
        X_train (pd.DataFrame): Training set
        X_test (pd.DataFrame): Testing set
        y_train (pd.series): Training label
        y_test (pd.series): Testing label
        list_interactions (list): List of tuples representing the interactions to create, 
        or list of strings such as ["RBP1*RBP2", "RBP2*RBP4"]. 
        split_symbol (str, optional): symbol to split and retrieve the two interacting terms in 
        the strings. Defaults to '*'.
        keep_old (bool, optional): Whether to keep or not the original columns. Defaults to True.

    Returns:
        : The updated dataset
    """

    X_train_int = get_interactions_dataset_from_combs(original_X=X_train, 
                                                      list_interactions=list_interactions,
                                                      split_symbol=split_symbol)
    X_test_int = get_interactions_dataset_from_combs(original_X=X_test, 
                                                    list_interactions=list_interactions,
                                                    split_symbol=split_symbol)
    if not keep_old:
        X_train_int.drop(X_train.columns, axis=1, inplace=True)
        X_test_int.drop(X_test.columns, axis=1, inplace=True)
    
    return X_train_int, X_test_int, y_train, y_test
        
    
def plot_shapley(shap_values, X_test, feature_names, max_display=44, show=False, class_label=1, 
                 path_to_save=None, negative_regulators=[], condition_name='', split_symbol='*', xlim=None, sort=True, 
                 cmap='bwr', custom_legend=True, 
                 format='png', color_bar=False, beeswarm_only=False):
    """Plot Shapley values as beeswarm plot and barplot.
    
    Args:
        shap_values (array-like): vector of shapley values. This is a list of matrices of SHAP values: (# samples x # features)
        X_test (array-like): Testing vector
        feature_names (_type_): Feature names
        max_display (int, optional): How many top features to include in the plot. Defaults to 44.
        show (bool, optional): Whether to show the plot. Defaults to False.
        class_label (int, optional): Class label for which you want to plot the SHAP values. Defaults to 1.
        path_to_save (str, optional): Path to save the figures. Defaults to None.
        negative_regulators (list, optional): List of feature names to write in bold. Defaults to [].
        condition_name (str, optional): Name of the condition that will be used for the titles. Defaults to ''.
        split_symbol (str, optional): Symbol to split feature names (in case of interaction terms). Defaults to '*'.
        xlim (float, optional): The plot will be displayed between -xlim and xlin on the x-axis. Defaults to None.
        sort (bool, optional): Whether to sort the features per feature importance in the beeswarm plot. Defaults to True.
        cmap (str, optional): Color map to use for the beeswarm plot. Defaults to 'bwr'.
        custom_legend (bool, optional): Whether to use a custom legend (for categorical values). Defaults to True.
        format (str, optional): The format in which to save the figures. Defaults to 'png'.
        color_bar (boo, optional): Whether to display the color bar on the beeswarm plots.
        beeswarm_only: whether to plot the beeswarm only or to display also the barplot.
    """
    
    fig, ax = plt.subplots(figsize=(10, 5))
    shap.summary_plot(shap_values[class_label], X_test, feature_names=feature_names, max_display=max_display, 
                      show=show, sort=sort, color_bar=color_bar, cmap=cmap)
    tick_labels = []
    for feature_name in ax.get_yticklabels():
        feature_name_short = feature_name.get_text().split('=')[0]
        feature_tick_name = []
        for feature in feature_name_short.split(split_symbol):
            if feature in negative_regulators:
                feature_tick_name.append(r"$\bf{"+feature+"}$")
            else:
                feature_tick_name.append(feature)
        tick_labels.append(split_symbol.join(feature_tick_name))
    ax.set_yticks(ticks = np.arange(len(feature_names)), 
            labels = tick_labels)
    if xlim is not None:
        plt.xlim(xlim)
        
    if custom_legend:
        from matplotlib.patches import Patch
        cmap_col = plt.get_cmap(cmap)
        custom_patches = [Patch(facecolor=cmap_col(0.), label='0'),
                        Patch(facecolor=cmap_col(1.), label='1')]
        
        ax.legend(handles=custom_patches, loc='best')
    ax.set_title(condition_name,  weight='bold')

    if path_to_save is not None:
        fig.savefig(os.path.join(path_to_save, "beeswarm_"+condition_name.replace(' ', '_')+"."+format), bbox_inches='tight', format=format)


    if not beeswarm_only:
        fig, ax = plt.subplots()
        shap.summary_plot(shap_values[class_label], X_test, feature_names=feature_names, plot_type='bar', max_display=max_display, 
                        show=show, sort=sort)
        tick_labels = []
        for feature_name in ax.get_yticklabels():
            feature_name_short = feature_name.get_text().split('=')[0]
            feature_tick_name = []
            for feature in feature_name_short.split(split_symbol):
                if feature in negative_regulators:
                    feature_tick_name.append(r"$\bf{"+feature+"}$")
                else:
                    feature_tick_name.append(feature)
            tick_labels.append(split_symbol.join(feature_tick_name))
        ax.set_yticks(ticks = np.arange(len(feature_names)), 
                labels = tick_labels);
        ax.set_title(condition_name,  weight='bold')

        if path_to_save is not None:
            fig.savefig(os.path.join(path_to_save, "barplot_"+condition_name.replace(' ', '_')+"."+format), bbox_inches='tight', format=format)

def plot_one_score_according_to_classif(df_data, LS_score_col, y_true, y_pred, sample_names, 
                                        return_classif=False, y_label=None):
    """_Plot the scores distributions colored by percentage of false positives, false negatives, true positives and true 
    negatives.

    Args:
        df_data (pd.DataFrame): The DataFrame with samples in index and the wanted score in index.
        LS_score_col (str): Name of the column bearing the score you want to plot.
        y_true (array-like): True labels
        y_pred (array-like): Predicted lables
        sample_names (array-like): Sample names in the testing set.
    """

    test_results = pd.DataFrame({"y_pred": y_pred, "y_true": y_true}, index=sample_names)
    
        # False negatives: those that are classified as Negatives but are infact positives
    ## y_pred = 0; y_true = 1
    falses_negatves = test_results[(test_results['y_pred'] == 0) & (test_results['y_true'] == 1)]


    # False positives: those that are classified as Positives but are in fact negatives
    ## y_pred = 1; y_true = 0
    false_positives = test_results[(test_results['y_pred'] == 1) & (test_results['y_true'] == 0)]

    ## True positives
    true_positives = test_results[(test_results['y_pred'] == 1) & (test_results['y_true'] == 1)]

    ## True negatives
    true_negatives = test_results[(test_results['y_pred'] == 0) & (test_results['y_true'] == 0)]
    
    sns.histplot(df_data.loc[true_positives.index], x=LS_score_col, binwidth=0.1, label='true positives', alpha=0.2, stat='percent', element='step')
    sns.histplot(df_data.loc[true_negatives.index], x=LS_score_col, binwidth=0.1, label='true negatives', alpha=0.2, stat='percent', element='step')
    sns.histplot(df_data.loc[falses_negatves.index], x = LS_score_col, binwidth=0.1, label='false negatives', alpha=0.2, stat='percent', element='step')
    sns.histplot(df_data.loc[false_positives.index], x=LS_score_col, binwidth=0.1, label='false positives', alpha=0.2, stat='percent', element='step')
    if y_label is None:
        y_label = LS_score_col
    plt.ylabel(y_label)
    plt.legend()
    
    if return_classif:
        return [[true_positives, falses_negatves], 
                [false_positives, true_negatives]]

def plot_one_score_according_to_classif_per_condition(df_data, datasets_per_condition, dataset_key='dataset', model_key='logit_model',
                                                      conditions_names=None, save_path=None, format='pdf', exp_name='', transparent=True, return_classif=False, ylabels=None):
    """Plot the scores distrbutions colored by percentage of palse positives, false negatives, true positives and true negatives.

    Args:
        df_data (pd.DataFrame): The DataFrame with samples in index and the wanted score in index.
        datasets_per_condition (dict): dictionnary with conditions/treatments in keys and dictionnaries in values.
        dataset_key (str, optional): Key name of the dictionnaries in datasets_per_condition to retrieve the datasets. Defaults to 'dataset'.
        model_key (str, optional): Key name of the dictionnaries in datasets_per_condition to retrieve the models. Defaults to 'logit_model'.
        conditions_names (list, optional): short treatment names for titles. Defaults to None, will use the treatment keys.
    """
    
    if conditions_names is None:
        conditions_names = list(datasets_per_condition.keys())
    if ylabels is None:
        ylabels = list(datasets_per_condition.keys())
    
    plt.figure(figsize=(5*len(datasets_per_condition), 5))
    for i, condition in enumerate(datasets_per_condition.keys()):
        plt.subplot(1, len(datasets_per_condition), i+1)
        _, X_test, _, y_test = datasets_per_condition[condition][dataset_key].values()
        model = datasets_per_condition[condition][model_key]
        y_pred = model.predict(X_test)
        if return_classif:
            
            
            datasets_per_condition[condition][model_key+"_"+dataset_key+"_classif_matrix"] = plot_one_score_according_to_classif(df_data=df_data, 
                                                                                                                                 LS_score_col=condition, 
                                                                                                                                 y_true=y_test, 
                                                                                                                                 y_pred=y_pred, 
                                                                                                                                 sample_names=X_test.index, 
                                                                                                                                 return_classif=return_classif)
        else:
            plot_one_score_according_to_classif(df_data=df_data, 
                                                LS_score_col=condition, 
                                                y_true=y_test, 
                                                y_pred=y_pred, 
                                                sample_names=X_test.index, 
                                                return_classif=return_classif, 
                                                y_label=ylabels[i])
            
        plt.title(conditions_names[i], weight='bold')
    plt.tight_layout()
    
    if save_path is not None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        plt.savefig(os.path.join(save_path, exp_name+"_scores_distributions_classification."+format), 
                    bbox_inches='tight', format=format, transparent=transparent)