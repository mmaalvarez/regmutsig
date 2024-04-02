# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
# Data Processing
import pandas as pd
import numpy as np

# Modelling
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint

# Tree Visualisation
from sklearn.tree import export_graphviz
from IPython.display import Image
import graphviz

# arguments
import sys
import argparse
import os
from glob import glob

# +
## parse arguments

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir_path', type=str)
parser.add_argument('-p', '--NMF_parameters', type=str, help='nFact and k values')
parser.add_argument('-g', '--alteration', type=str, help='Alteration for which we want to predict whether a sample is altered or WT')
parser.add_argument('-s', '--seed', type=int, help='Seed for all functions that involve some randomness: train_test_split, RandomForestClassifier, RandomizedSearchCV')

# if interactive, pass some values manually
if 'ipykernel' in sys.modules:
    args = parser.parse_args(['-i', '/home/jovyan/fsupek_data/users/malvarez/projects/RepDefSig/models/model_2/2_signature_extraction_and_plotting/regional+SBS/NMF/random_forest/rf_inputs', 
                              '-p', 'nFact-11_k-11', 
                              '-g', 'Nitrosamines', 
                              '-s', '0'])
else:
    args = parser.parse_args()
    
input_dir_path, NMF_parameters, alteration, seed = args.input_dir_path, args.NMF_parameters, args.alteration, args.seed

test_size = 0.2
# -

# Load data

# +
# predict single prediction_scope
for f in glob(f'{input_dir_path}/res_{NMF_parameters}/*/{alteration}_exposures.tsv'):
    prediction_scope = f.split('/')[-2]
    exposures = pd.read_csv(f, sep="\t")

signature_names = exposures.drop(['id','altered_sample'], axis=1).keys().tolist()
# -

# Downsample the normal samples (i.e. with a '0') to match the number of case samples ('1's)
#
# Since I do `sample(frac=0.99-fraction_cases)` there should be slightly more 0's than 1's, still

number_samples = exposures.shape[0]
number_cases = exposures[exposures['altered_sample']==1].shape[0]
fraction_cases = number_cases/number_samples
exposures_normal_cases_matched = exposures.drop(exposures[exposures['altered_sample'] == 0].sample(frac=0.99-fraction_cases).index)

# Split data into features (X) and target (y)

X = exposures_normal_cases_matched.drop(['id','altered_sample'], axis=1)
y = exposures_normal_cases_matched['altered_sample']

# Split both X and y into training and test data
# - Training data is used to fit the model. The algorithm uses the training data to learn the relationship between the features and the target
# - Test data is used to evaluate the performance of the model
#   - test_size indicates the fraction of the original data to be used as test (its complement will be the fraction used as training data)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)

# + [markdown] tags=[]
# Hyperparameter Tuning
# ---------------------
#
# RandomizedSearchCV randomly search parameters within a range per hyperparameter. We define the hyperparameters to use and their ranges in the param_dist dictionary. In our case, we are using:
#
# - n_estimators: the number of decision trees in the forest. Increasing this hyperparameter prediction_scoperally improves the performance of the model but also increases the computational cost of training and predicting.
# - max_depth: the maximum depth of each decision tree in the forest. Setting a higher value for max_depth can lead to overfitting while setting it too low can lead to underfitting
# - random_state is a seed for reproducing results
#
# RandomizedSearchCV will train many models (defined by n_iter_ and save each one as variables, the code below creates a variable for the best model and prints the hyperparameters. In this case, we haven’t passed a scoring system to the function, so it defaults to accuracy. This function also uses cross validation, which means it splits the data into five equal-sized groups and uses 4 to train and 1 to test the result. It will loop through each group and give an accuracy score, which is averaged to find the best model.

# +
param_dist = {'n_estimators': randint(50,500),
              'max_depth': randint(1,20)}

## Create an instance of the Random Forest model, with the default parameters
rf = RandomForestClassifier(random_state=seed,
                            oob_score=True)

## Use random search to find the best hyperparameters
rand_search = RandomizedSearchCV(rf, 
                                 param_distributions = param_dist, 
                                 n_iter=5, 
                                 cv=5,
                                 random_state=seed)

## Fit the random search object to our training data. We pass both the features and the target variable, so the model can learn.
rand_search.fit(X_train, y_train)

## Create a variable for the best model
best_rf = rand_search.best_estimator_

## Print the best hyperparameters
#print('Best hyperparameters:', rand_search.best_params_)
# -

# At this point, we have a trained Random Forest model, but we need to find out whether it is making accurate predictions.

## prediction_scoperate predictions with the best model
y_pred = best_rf.predict(X_test)

# Evaluate this model using *out-of-bag* score; also can do accuracy, precision, and recall; we check the predictions against the actual values in the test set and count up how many the model got right

oob = best_rf.oob_score_

# The confusion matrix plots what the model predicted against what the correct prediction was. We can use this to understand the tradeoff between false positives (top right) and false negatives(bottom left)

# +
# Create the confusion matrix
#cm = confusion_matrix(y_test, y_test)

# plot it
#ConfusionMatrixDisplay(confusion_matrix=cm,display_labels=True).plot()
# -

# Plot the importance of each feature, using the model’s internal score to find the best way to split the data within each decision tree

# +
## Create a series containing feature importances from the model and feature names from the training data
feature_importances = pd.Series(best_rf.feature_importances_, index=X_train.columns).sort_values(ascending=False)

#feature_importances.plot.bar()
# -

## print output table
res = pd.DataFrame({'NMF_parameters': [NMF_parameters],
                    'test_size': [test_size],
                    'prediction_scope': [prediction_scope],
                    'alteration': [alteration],
                    'seed': [seed],
                    'oob': [oob],
                    'n_cases': [number_cases]})

# +
## write results table

#list(res.columns)

res.to_csv(f'res_{NMF_parameters}_{prediction_scope}_{alteration}_{seed}.csv', header=True, sep='\t', index=False)
# -

# Visualize the first 3 trees
# ---------------------------
#
# Each tree image is limited to only showing the first few nodes. These trees can get very large and difficult to visualize. The colors represent the majority class of each node (boxes). The colors get darker the closer the node gets to being fully a category. Each node also contains the following information:
#
# 1- The variable name and value used for splitting
#
# 2- The % of total samples in each split
#
# 3- The % split between classes in each split

'''
for i in range(3):
    tree = best_rf.estimators_[i]
    dot_data = export_graphviz(tree,
                               feature_names=X_train.columns,  
                               filled=True,  
                               max_depth=2, 
                               impurity=False, 
                               proportion=True)
    graph = graphviz.Source(dot_data)
    display(graph)
'''
