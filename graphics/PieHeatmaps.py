### this script makes colored pie charts representing activity with each substrate for each sequence; pies are stored as images in a directory coded by sequence

### import libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from palettable.scientific.sequential import Bilbao_20

### set working family
fam = 'Fam1B.1'


### read and sort data
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv').sort_values(by=[fam,'all_sum'], ascending=[True,False])

### get log of catalytic enhancement values for color map
df['log10-r_BWO'] = np.log10(df['r_BWO'])
df['log10-r_BFO'] = np.log10(df['r_BFO'])
df['log10-r_BLO'] = np.log10(df['r_BLO'])
df['log10-r_BIO'] = np.log10(df['r_BIO'])
df['log10-r_BVO'] = np.log10(df['r_BVO'])
df['log10-r_BMO'] = np.log10(df['r_BMO'])

### mask values below a set value
bg = 1
df['log10-r_BWO'].mask(df['log10-r_BWO'] <= np.log10(bg), inplace=True)
df['log10-r_BFO'].mask(df['log10-r_BFO'] <= np.log10(bg), inplace=True)
df['log10-r_BLO'].mask(df['log10-r_BLO'] <= np.log10(bg), inplace=True)
df['log10-r_BIO'].mask(df['log10-r_BIO'] <= np.log10(bg), inplace=True)
df['log10-r_BVO'].mask(df['log10-r_BVO'] <= np.log10(bg), inplace=True)
df['log10-r_BMO'].mask(df['log10-r_BMO'] <= np.log10(bg), inplace=True)

### set distance and return only log activity values
Family = df[df['Family'] == fam]
working = Family[['log10-r_BWO','log10-r_BFO','log10-r_BLO','log10-r_BIO','log10-r_BVO','log10-r_BMO']]

### make list of sequences
seq_list = Family['seq'].tolist()

### determine maximum and minimum r values for all substrates; MaxMax will be used as the upper limit for colormapping
W_max = df['log10-r_BWO'].max()
F_max = df['log10-r_BFO'].max()
L_max = df['log10-r_BLO'].max()
I_max = df['log10-r_BIO'].max()
V_max = df['log10-r_BVO'].max()
M_max = df['log10-r_BMO'].max()
maximums = [W_max,F_max,L_max,I_max,V_max,M_max]
MaxMax = max(maximums)

W_min = df['log10-r_BWO'].min()
F_min = df['log10-r_BFO'].min()
L_min = df['log10-r_BLO'].min()
I_min = df['log10-r_BIO'].min()
V_min = df['log10-r_BVO'].min()
M_min = df['log10-r_BMO'].min()
minimums = [W_min,F_min,L_min,I_min,V_min,M_min]
MinMin = min(minimums)

print('Maximum value: ',MaxMax)
print('Minimum value: ',MinMin)

### set color map, labels, and slice sizes
cmap = Bilbao_20.mpl_colormap
cmap.set_bad(cmap(np.log10(1)/MaxMax))
labels = 'BWO','BFO','BLO','BIO','BVO','BMO'
sizes = [1,1,1,1,1,1]

### make save directories
dir = 'outputs/pies/' + fam + '/'
if not os.path.exists('outputs/pies/'):
    os.mkdir('outputs/pies/')
if not os.path.exists(dir):
    os.mkdir(dir)

### iterate rows and make pie charts
fig_num = 0
plt.figure(1)
for id, row in working.iterrows():
    values = np.array(list(row))
    fig_name = dir + seq_list[fig_num] + '.png'
    fig_num += 1
    colors = cmap(row/MaxMax)
    plot = plt.pie(sizes, colors=colors, counterclock=False, startangle=150, textprops={'fontsize': 24}, wedgeprops = {'edgecolor': 'black', 'linewidth': 1})#, labels=labels,)
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(6,6)
    plt.savefig(fig_name, transparent=True, dpi=300, format='png')