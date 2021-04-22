### this script concatenates the output files from the promiscuity index calculator (http://hetaira.herokuapp.com/) and merges I values with the master results table
### input files must be in an accessible directory (eg 'data/promiscuity_index_tables/WFLIVM-r_results_tables/') and must be in the format 'results (i).csv'

### import libraries
import pandas as pd


### create dataframe and
PI_DF = pd.DataFrame(columns=['seq','I'])
FileRange = range(0,40) ### set file range based on number of input files
for i in FileRange:
    FileName = 'data/promiscuity_index_tables/WFLIVM-r_results_tables/results (' + str(i) + ').csv'
    data = pd.read_csv(FileName, header=None, index_col=False)
    data.columns = ['seq','I']
    data.drop(data.tail(1).index,inplace=True)
    PI_DF = pd.concat([PI_DF,data], ignore_index=True)


### merge I values to master file
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r.csv').sort_values(by='seq')
merged = df.merge(PI_DF, on='seq')
merged.to_csv('data/WFLIVM-k-seq_merged_+r+I.csv', index=False)