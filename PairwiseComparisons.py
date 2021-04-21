### import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import *
import itertools
import numpy as np


### set working family (eg 'Fam1A.1' or 'AllFams')
family = 'AllFams'

linreg = False  ### set True to plot linear regression
plotlegend = True   ### set True to plot legend


### import csv dataset as dataframe
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv').sort_values(by='all_sum', ascending=False)

working = df[df['Family'] == family]
if family == 'AllFams':
    working = df


### define layout of subplots
FigNum = [21,16,11,6,1,17,12,7,2,13,8,3,9,4,5]
fn=0


SubNames = ['r_BWO','r_BFO','r_BLO','r_BIO','r_BVO','r_BMO']

AllFitsDF = pd.DataFrame(columns=['Sub1','Sub2','slope','intercept'])

### plot data
for i,j in list(itertools.combinations(SubNames,2)):
    plt.subplot(5,5,FigNum[fn])
    plt.scatter(x=working[i], y=working[j], marker='.', c=working['color_label'], alpha=0.5, edgecolors=[], zorder=4)

    plt.xlim(0.15,1500)
    plt.ylim(0.15,1500)
    plt.xscale('log')
    plt.yscale('log')

    plt.tick_params(axis='both', which='both', direction='in', labelbottom=False, labelleft=False)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.plot([10000,0.0001], [10000, 0.0001], color='k', linestyle='--', linewidth=1.5, alpha=0.25)
    if FigNum[fn] in [1,6,11,16,21]:
        plt.tick_params(labelleft=True)
        plt.yticks(fontsize=16)
    if FigNum[fn] in [21,17,13,9,5]:
        plt.tick_params(labelbottom=True)
        plt.xticks(fontsize=16)
    LabelSize = 32
    LabelPad = 20
    if FigNum[fn] == 1:
        plt.xlabel("$r_{s}^{BWO}$", fontsize=LabelSize, labelpad=LabelPad)
        ax.xaxis.set_label_position('top')
    if FigNum[fn] == 2:
        plt.xlabel("$r_{s}^{BFO}$", fontsize=LabelSize, labelpad=LabelPad)
        ax.xaxis.set_label_position('top')
    if FigNum[fn] == 3:
        plt.xlabel("$r_{s}^{BLO}$", fontsize=LabelSize, labelpad=LabelPad)
        ax.xaxis.set_label_position('top')
    if FigNum[fn] == 4:
        plt.xlabel("$r_{s}^{BIO}$", fontsize=LabelSize, labelpad=LabelPad)
        ax.xaxis.set_label_position('top')
    if FigNum[fn] == 5:
        plt.xlabel("$r_{s}^{BVO}$", fontsize=LabelSize, labelpad=LabelPad)
        ax.xaxis.set_label_position('top')
    if FigNum[fn] == 21:
        plt.ylabel("$r_{s}^{BFO}$", fontsize=LabelSize)
    if FigNum[fn] == 16:
        plt.ylabel("$r_{s}^{BLO}$", fontsize=LabelSize)
    if FigNum[fn] == 11:
        plt.ylabel("$r_{s}^{BIO}$", fontsize=LabelSize)
    if FigNum[fn] == 6:
        plt.ylabel("$r_{s}^{BVO}$", fontsize=LabelSize)
    if FigNum[fn] == 1:
        plt.ylabel("$r_{s}^{BMO}$", fontsize=LabelSize)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    fn += 1

### linear regression fitting
    if linreg == True:
        m,b = np.polyfit(working[i], working[j],1)
        plt.plot(np.arange(0.1,2000,0.1), m * np.arange(0.1,2000,0.1) + b, color='k', linewidth=2, zorder=10)
        FitDF = pd.DataFrame([[i,j,m,b]], columns=['Sub1','Sub2','slope','intercept'])
        AllFitsDF = AllFitsDF.append(FitDF)

plt.subplot(5,5,19) ### set position of legend
plt.axis('off')
palette = ['#1E3888','#47A8BD','#66101F','#F32E2B','#D3BF0D']
leg21 = mpatches.Patch(color=palette[0], label='2.1', alpha=0.75)
leg22 = mpatches.Patch(color=palette[1], label='2.2', alpha=0.75)
leg1A = mpatches.Patch(color=palette[2], label='1A.1', alpha=0.75)
leg1B = mpatches.Patch(color=palette[3], label='1B.1', alpha=0.75)
leg3 = mpatches.Patch(color=palette[4], label='3.1', alpha=0.75)
if plotlegend == True:
    plt.legend(handles=[leg1A, leg1B, leg21, leg22, leg3], prop={'size': 28}, loc='best', title='Family', title_fontsize=32)

fig = plt.gcf()
fig.set_size_inches(20,20)

if linreg == True:
    AllFitsDF.to_csv('outputs/' + family + '_r-LinRegFits.csv')

plt.savefig('outputs/' + family + '_r-id_leg_TL.png', dpi=200)