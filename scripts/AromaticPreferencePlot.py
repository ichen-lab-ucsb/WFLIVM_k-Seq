### import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
import itertools
from scipy import stats


### working parameters
Family = 'AllFams'  ### set working family to plot (eg 'Fam1A.1' or 'AllFams')
CompTable = ['all_sum','aromatic/all_sum']  ### choose values to compare; will need to also set axes labels on plot and save names
leg = True  ### will add legend to plot if True


### import csv dataset and define dataframes
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv').sort_values(by='all_sum', ascending=True)

Fam1A = df[df['Family'] == 'Fam1A.1']
Fam1B = df[df['Family'] == 'Fam1B.1']
Fam21 = df[df['Family'] == 'Fam2.1']
Fam22 = df[df['Family'] == 'Fam2.2']
Fam3 = df[df['Family'] == 'Fam3.1']

WorkingDict = {'AllFams':df,'Fam1A.1':Fam1A,'Fam1B.1':Fam1B,'Fam2.1':Fam21,'Fam2.2':Fam22,'Fam3.1':Fam3}
working = WorkingDict[Family]


### calculate Pearson and Spearman correlation coefficients and p-values
fams = [Fam1A,Fam1B,Fam21,Fam22,Fam3]
FamNames = ['Fam1A.1','Fam1B.1','Fam2.1','Fam2.2','Fam3.1']
r_table = pd.DataFrame(columns=['Family','comp1','comp2','Pearson_r','Pearson_p','Spearman_r','Spearman_p'])
f = 0
for fam in fams:
    FamName = FamNames[f]
    f += 1
    for i,j in list(itertools.combinations(CompTable,2)):
        pearson = stats.pearsonr(fam[i], fam[j])
        spearman = stats.spearmanr(fam[i], fam[j])
        r_data = pd.DataFrame([[FamName,i,j,pearson[0],pearson[1],spearman[0],spearman[1]]], columns=['Family','comp1','comp2','Pearson_r','Pearson_p','Spearman_r','Spearman_p'])
        r_table = r_table.append(r_data, ignore_index=True)


### plot data
plt.scatter(x=working[CompTable[0]], y=(working[CompTable[1]]), marker='.', c=working['color_label'], alpha=0.5, edgecolors=[], zorder=4)

plt.xscale('log')
plt.xlabel("$\sum_X r_{s}^{BXO}$", fontsize=16) ### change if necessary
plt.ylabel("Aromatic Preference", fontsize=16)  ### change if necessary
plt.xlim(2,2000)
plt.ylim(0,1.1)
plt.locator_params(axis='y', nbins=5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig = plt.gcf()
fig.set_size_inches(6,4)
plt.tight_layout()

### plot legend
if leg == True:
    palette = ['#1E3888','#47A8BD','#66101F','#F32E2B','#D3BF0D']
    leg21 = mpatches.Patch(color=palette[0], label='2.1', alpha=0.75)
    leg22 = mpatches.Patch(color=palette[1], label='2.2', alpha=0.75)
    leg1A = mpatches.Patch(color=palette[2], label='1A.1', alpha=0.75)
    leg1B = mpatches.Patch(color=palette[3], label='1B.1', alpha=0.75)
    leg3 = mpatches.Patch(color=palette[4], label='3.1', alpha=0.75)
    plt.legend(handles=[leg1A, leg1B, leg21, leg22, leg3], prop={'size': 12}, loc='best', title='Family', title_fontsize=14)


### save outputs
plt.savefig('outputs/' + Family + '_AllSum,AroPref_plot.png', dpi=500, format='png')
r_table.to_csv('outputs/AllFams_AllSum,AroPref_pearson,spearman_table.csv', index=False)