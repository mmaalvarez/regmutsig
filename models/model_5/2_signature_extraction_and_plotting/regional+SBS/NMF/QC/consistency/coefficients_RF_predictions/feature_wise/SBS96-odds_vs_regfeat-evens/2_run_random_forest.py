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
import sys
import argparse
import os
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir_path', type=str)
parser.add_argument('--predictor', type=str)
parser.add_argument('--topredict', type=str)

# if interactive, pass some values manually
if 'ipykernel' in sys.modules:
    args = parser.parse_args(['-i', '/home/jovyan/fsupek_data/users/malvarez/projects/RepDefSig/models/model_4/2_signature_extraction_and_plotting/regional+SBS/NMF/QC/consistency/coefficients_RF_predictions/rf_inputs/', 
                              '--predictor', '1st_half__regional_feature__original', 
                              '--topredict', '1st_half__regional_feature__resampled'])
else:
    args = parser.parse_args()
    
input_dir_path, predictor, topredict = args.input_dir_path, args.predictor, args.topredict

## load data
predictor_coeff = pd.read_csv(f'{input_dir_path}/predictor-{predictor}___topredict-{topredict}/predictor.tsv', sep="\t").drop(['id'], axis=1).fillna(0)
topredict_coeff = pd.read_csv(f'{input_dir_path}/predictor-{predictor}___topredict-{topredict}/topredict.tsv', sep="\t").drop(['id'], axis=1).fillna(0)

## fit model
model = LinearRegression()

model.fit(predictor_coeff, topredict_coeff)

## calculate R2
prediction = model.predict(predictor_coeff)

R2 = r2_score(topredict_coeff, prediction)

## print
res = pd.DataFrame({'predictor': [predictor],
                    'topredict': [topredict],
                    'R2': [R2]})

res.to_csv(f'res__predictor-{predictor}__topredict-{topredict}.tsv', header=True, sep='\t', index=False)
