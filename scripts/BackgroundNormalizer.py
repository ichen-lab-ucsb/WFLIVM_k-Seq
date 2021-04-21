### import libraries
import pandas as pd

#read data
df = pd.read_csv('data/WFLIVM-k-seq_merged.csv').sort_values(by='seq')

### set background rates
BGDict = {'BWO':1.57,'BFO':1.37,'BLO':1.46,'BIO':0.87,'BVO':1.21,'BMO':2.14}

### normalize values to background rate
W_norm = []
F_norm = []
L_norm = []
I_norm = []
V_norm = []
M_norm = []
Y_norm = []

for i in df['bs_kA_50%_BWO']:
    norm = (i / BGDict['BWO'])
    W_norm.append(norm)
for i in df['bs_kA_50%_BFO']:
    norm = (i / BGDict['BFO'])
    F_norm.append(norm)
for i in df['bs_kA_50%_BLO']:
    norm = (i / BGDict['BLO'])
    L_norm.append(norm)
for i in df['bs_kA_50%_BIO']:
    norm = (i / BGDict['BIO'])
    I_norm.append(norm)
for i in df['bs_kA_50%_BVO']:
    norm = (i / BGDict['BVO'])
    V_norm.append(norm)
for i in df['bs_kA_50%_BMO']:
    norm = (i / BGDict['BMO'])
    M_norm.append(norm)

df['r_BWO'] = W_norm
df['r_BFO'] = F_norm
df['r_BLO'] = L_norm
df['r_BIO'] = I_norm
df['r_BVO'] = V_norm
df['r_BMO'] = M_norm

### average and sum r values and calculate ratios
aromatic_hmean = 2/((1/df['r_BWO'])+(1/df['r_BFO']))
nonaromatic_hmean = 4/((1/df['r_BLO'])+(1/df['r_BIO'])+(1/df['r_BVO'])+(1/df['r_BMO']))
all_hmean = 6/((1/df['r_BWO'])+(1/df['r_BFO'])+(1/df['r_BLO'])+(1/df['r_BIO'])+(1/df['r_BVO'])+(1/df['r_BMO']))

aromatic_sum = df['r_BWO']+df['r_BFO']
nonaromatic_sum = df['r_BLO']+df['r_BIO']+df['r_BVO']+df['r_BMO']
all_sum = df['r_BWO']+df['r_BFO']+df['r_BLO']+df['r_BIO']+df['r_BVO']+df['r_BMO']

ratio1 = (aromatic_hmean/nonaromatic_hmean)
ratio2 = (nonaromatic_hmean/aromatic_hmean)
ratio3 = (aromatic_hmean/all_hmean)
ratio4 = (aromatic_sum/nonaromatic_sum)
ratio5 = (nonaromatic_sum/aromatic_sum)
ratio6 = (aromatic_sum/all_sum)

df['aromatic_mean'] = aromatic_hmean
df['nonaromatic_mean'] = nonaromatic_hmean
df['all_mean'] = all_hmean
df['aromatic_sum'] = aromatic_sum
df['nonaromatic_sum'] = nonaromatic_sum
df['all_sum'] = all_sum
df['aromatic/non_mean'] = ratio1
df['non/aromatic_mean'] = ratio2
df['aromatic/all_mean'] = ratio3
df['aromatic/non_sum'] = ratio4
df['non/aromatic_sum'] = ratio5
df['aromatic/all_sum'] = ratio6

### save file
df.to_csv('data/WFLIVM-k-seq_merged_+r.csv', index=False)