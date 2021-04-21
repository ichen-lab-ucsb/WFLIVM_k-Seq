### import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import numpy as np
import itertools
from pylab import *

### plot families together or separately
subplots = False ### if True, will plot each Family as separate plots; if False, will plot all Families together


### import csv dataset and get log values
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv').sort_values(by='all_sum', ascending=True)

### data subsets
f1A = df[df['Family'] == 'Fam1A.1']
f1B = df[df['Family'] == 'Fam1B.1']
f21 = df[df['Family'] == 'Fam2.1']
f22 = df[df['Family'] == 'Fam2.2']
f3 = df[df['Family'] == 'Fam3.1']


### Calculate correlations
fams = [f1A,f1B,f21,f22,f3]
FamNames = ['1A.1','1B.1','2.1','2.2','3.1']
r_I = ['all_sum','I']
r_table = pd.DataFrame(columns=['Family','comp1','comp2','Pearson_r','Pearson_p','Spearman_r','Spearman_p'])
f = 0
for fam in fams:
    FamName = FamNames[f]
    f += 1
    for i,j in list(itertools.combinations(r_I,2)):
        pearson = stats.pearsonr(fam[i], fam[j])
        spearman = stats.spearmanr(fam[i], fam[j])
        r_data = pd.DataFrame([[FamName,i,j,pearson[0],pearson[1],spearman[0],spearman[1]]], columns=['Family','comp1','comp2','Pearson_r','Pearson_p','Spearman_r','Spearman_p'])
        r_table = r_table.append(r_data, ignore_index=True)


### for legend
palette = ['#1E3888','#47A8BD','#66101F','#F32E2B','#D3BF0D']
leg21 = mpatches.Patch(color=palette[0], label='2.1', alpha=0.75)
leg22 = mpatches.Patch(color=palette[1], label='2.2', alpha=0.75)
leg1A = mpatches.Patch(color=palette[2], label='1A.1', alpha=0.75)
leg1B = mpatches.Patch(color=palette[3], label='1B.1', alpha=0.75)
leg3 = mpatches.Patch(color=palette[4], label='3.1', alpha=0.75)


### plot data
if subplots == True:
    n = 0
    for i in fams:
        n+=1
        plt.subplot(3,2,n)
        x = i['all_sum']
        y = i['I']
        df = i
        plt.scatter(x=x, y=y, marker='.', c=i['color_label'], alpha=0.5, edgecolors=[])
        ax = plt.gca()
        ax.set_xscale('log')
        plt.xlabel("$\sum_X r_{s}^{BXO}$", fontsize=16)
        plt.ylabel('Promiscuity Index ($I_s$)',fontsize=16)
        plt.xlim(2,2000)
        plt.ylim(0,1.1)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        fig = plt.gcf()
        plt.title('Family ' + FamNames[n - 1], fontsize=18)
        fig.set_size_inches(12, 12)

if subplots == False:
    x = df['all_sum']
    y = df['I']
    plt.scatter(x=x, y=y, marker='.', c=df['color_label'], alpha=0.5, edgecolors=[])
    ax = plt.gca()
    ax.set_xscale('log')
    plt.xlabel("$\sum_X r_{s}^{BXO}$", fontsize=16)
    plt.ylabel('Promiscuity Index ($I_s$)',fontsize=16)
    plt.xlim(2,2000)
    plt.ylim(0,1.1)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(6,4)
    plt.legend(handles=[leg1A, leg1B, leg21, leg22, leg3], prop={'size': 12}, loc='best', title='Family',title_fontsize=14)

plt.tight_layout()


### save data
if subplots == True:
    plt.savefig('outputs/WFLIVM_sum-r_v_WFLIVM-I_EachFam.png', dpi=500, format='png')
if subplots == False:
    plt.savefig('outputs/WFLIVM_sum-r_v_WFLIVM-I_AllFams.png', dpi=500, format='png')

r_table.to_csv('outputs/WFLIVM_sum-r_I_correlations.csv', index=False)