### import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *

### import csv dataset as dataframe
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv').sort_values(by='seq')

palette = ['#1E3888','#47A8BD','#66101F','#F32E2B','#D3BF0D']
Families = ['Fam1A.1','Fam1B.1','Fam2.1','Fam2.2','Fam3.1']
subs = ['BWO','BFO','BLO','BIO','BVO','BMO']


### calculate precision and plot
n = 0
for i in subs:
    PrecisionList = []
    for index,row in df.iterrows():
        Precision = (row['bs_kA_97.5%_' + i] - row['bs_kA_2.5%_' + i]) / row['bs_kA_50%_' + i]
        PrecisionList.append(Precision)

    df['Precision_' + i] = PrecisionList
    plt.subplot(3,2,n+1)
    plt.scatter(x=df['bs_kA_50%_' + i], y=df['Precision_' + i], marker='.', c=df['color_label'], alpha=0.5, edgecolors=[])

    plt.xscale('log')
    plt.yscale('log')
    plt.title(i, fontsize=20)
    plt.xlim(0.2,1500)
    plt.ylim(0.08,200)
    plt.xlabel('$k_{s}A_{s}^{' + i + '}$ (Bootstrap median)', fontsize=14)
    plt.ylabel('95% CI range / median', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(11,11)

    if n == 0:
        leg21 = mpatches.Patch(color=palette[0], label='2.1', alpha=0.75)
        leg22 = mpatches.Patch(color=palette[1], label='2.2', alpha=0.75)
        leg1A = mpatches.Patch(color=palette[2], label='1A.1', alpha=0.75)
        leg1B = mpatches.Patch(color=palette[3], label='1B.1', alpha=0.75)
        leg3 = mpatches.Patch(color=palette[4], label='3.1', alpha=0.75)
        plt.legend(handles=[leg1A, leg1B, leg21, leg22, leg3], prop={'size': 12}, title='Family', title_fontsize=14)

    n += 1

plt.tight_layout()

plt.savefig('outputs/AllSubs_median-v-95CI;median_legend.png', format='png', dpi=300)