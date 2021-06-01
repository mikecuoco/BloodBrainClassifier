#!/usr/bin/env python

# To use this module, add this import statement at the top of your code:
# import scripts.starter as data
# Once you do that, all of the objects defined here will be accessible as
# data.object_name (ex: data.tpm or data.split())

# This module creates the following variables:
# tpm - A gene expression matrix for AD, N, and C patients (rows = samples, columns = genes)
# donors - Info about each individual in the study
# genes - Info about each gene in the tpm matrix
# samples - Info about each sample in the tpm matrix
# data - The TPM matrix subsetted to AD and N patients. The last column will contain each patient's status as AD or N
# X - the data matrix, but without its last column
# y - the last column of the data matrix
# X_train - raw gene expression matrix for training AD vs N patients
# X_test - raw gene expression matrix for testing AD vs N patients
# y_train - labels (1 for AD and 0 for N) corresponding to rows in X_train
# y_test - labels (1 for AD and 0 for N) corresponding to rows in X_test

# In addition, a function called split() is provided that can be used to update
# the values of X, y, X_train, X_test, y_train, y_test

import pandas as pd
from sklearn.model_selection import train_test_split

# load data
tpm = pd.read_table("data/GSE136243.txt.gz")
donors = pd.read_csv("data/GSE136243_donors.csv", index_col='Donor ID')

# move gene information to a separate dataframe
genes = tpm[['gene_id', 'gene_name', 'gene_biotype']].set_index('gene_id')
# remove genes and transpose for sklearn
tpm = tpm.set_index('gene_id').drop(genes.columns, axis=1).T
tpm.index.name = 'sample'
# extract labels:
# AD: alzheimer's
# N: control
# C: converter
tpm['group'] = tpm.index.str.split('_').str[0]

# extract sample info from the sample IDs
# supplementary table S3 is too difficult to parse into a table from a PDF :(
samples = pd.DataFrame(tpm.index)
samples['donor'] = samples['sample'].str.rsplit('_', 2).str[0]
samples['year'] = 2000 + samples['sample'].str.split('_').str[2].astype('int')
samples['group'] = samples['sample'].str.split('_').str[0]
samples = samples.join(donors[donors['Group'] == 'converter']['Year diagnosis'].astype('int'), on='donor')
samples.set_index('sample', inplace=True)

# subset to only the AD samples
# if True, this will include samples from converters after their diagnosis
if True:
    data = tpm[(samples['year'] > samples['Year diagnosis']) | (tpm['group'] != 'C')]
else:
    # remove samples coming from converters before their 
    data = tpm[tpm['group'] != 'C']

# ignore pesky warnings for the following line
pd.options.mode.chained_assignment = None  # (the default='warn')
# convert to boolean
data['group'].replace({'AD': 1, 'N': 0, 'C': 1}, inplace=True)

# define this code as a function so it can be called from outside the module, if needed
# note that calling this function will change the values of X_train, X_test, y_train,
# and y_test globally!
def split(test_size=0.33, random_state=42):
    ''' create a test/train split '''
    global data, X, y, X_train, X_test, y_train, y_test
    X = data.loc[:, data.columns != 'group']
    y = data['group']
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )

# split dataset into training and testing
split()
