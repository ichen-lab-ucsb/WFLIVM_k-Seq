### This script generates csv files to be used as inputs for calculating promiscuity indexes using the online tool available at http://hetaira.herokuapp.com/

### import libraries
import pandas as pd
import math
import os

### read data
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r.csv').sort_values(by='seq')

### set values to include in analysis
value_list = ['seq','r_BWO','r_BFO','r_BLO','r_BIO','r_BVO','r_BMO']

### set file size
FileLength = 245    ### this number is based on the maximum number of items that can be processed by the calculator in a single file

### make directory for output files
if not os.path.exists('data/promiscuity_index_tables/'):
    os.mkdir('data/promiscuity_index_tables/')
if not os.path.exists('data/promiscuity_index_tables/WFLIVM-r_input_tables/'):
    os.mkdir('data/promiscuity_index_tables/WFLIVM-r_input_tables/')

### make files
TotalFiles = range(0,math.ceil(len(df)/FileLength))
rows = 0
for i in TotalFiles:
    df_working = df[value_list].iloc[rows:rows+FileLength].set_index('seq').transpose()
    filename = 'data/promiscuity_index_tables/WFLIVM-r_input_tables/df' + str(i) + '.csv'
    df_working.to_csv(filename, index=False)
    rows += FileLength
