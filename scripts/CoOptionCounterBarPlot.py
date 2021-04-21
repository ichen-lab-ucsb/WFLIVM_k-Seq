### this script counts the number of sequences (Occurrences) with r values above a threshold value in x number (ActiveCount) of substrate columns

### import libraries
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


### set r threshold
rt = 5


### read data
df = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv')

### mask values below r threshold
df['r_BWO'].mask(df['r_BWO'] <= rt, inplace=True)
df['r_BFO'].mask(df['r_BFO'] <= rt, inplace=True)
df['r_BLO'].mask(df['r_BLO'] <= rt, inplace=True)
df['r_BIO'].mask(df['r_BIO'] <= rt, inplace=True)
df['r_BVO'].mask(df['r_BVO'] <= rt, inplace=True)
df['r_BMO'].mask(df['r_BMO'] <= rt, inplace=True)

### make family dataframes
f1A = df.loc[df['Family'] == 'Fam1A.1']
f1B = df.loc[df['Family'] == 'Fam1B.1']
f21 = df.loc[df['Family'] == 'Fam2.1']
f22 = df.loc[df['Family'] == 'Fam2.2']
f3 = df.loc[df['Family'] == 'Fam3.1']
Families = [f1A, f1B, f21, f22, f3]


### plot data
palette = ['#66101F','#F32E2B','#1E3888','#47A8BD','#D3BF0D']   ### color palette
w = 0.15    ### width of bars
PlotList = [-2*w,-w,0,w,2*w]    ### bar arrangement

for i in range(0,5):
    working = Families[i][['r_BWO', 'r_BFO', 'r_BLO', 'r_BIO', 'r_BVO', 'r_BMO']]
    ActiveCount = []
    for index,row in working.iterrows():
        ActiveCount.append(row.count())
    working['ActiveCount'] = ActiveCount
    ActiveCountCounts = working['ActiveCount'].value_counts()
    CountsDF = pd.DataFrame(data=list(zip(ActiveCountCounts.index.to_list(),ActiveCountCounts.to_list())), columns=['ActiveCount','Occurrences']).sort_values(by='ActiveCount',ascending=True)
    print(CountsDF)
    plt.bar(CountsDF['ActiveCount']+PlotList[i], (CountsDF['Occurrences']/1954), color=palette[i], width=w)
    plt.xlabel('#Substrates with $r_s$ > $r_t$', fontsize=18)
    plt.ylabel('Fraction of Family', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylim(0,1)
    plt.xlim(-0.5,6.5)
    plt.locator_params(axis='y', nbins=6)
    plt.locator_params(axis='x', nbins=7)
    plt.tight_layout()
    fig=plt.gcf()
    fig.set_size_inches(8,6)

### plot legend
leg1A = mpatches.Patch(color=palette[1], label='1A.1')
leg1B = mpatches.Patch(color=palette[0], label='1B.1')
leg21 = mpatches.Patch(color=palette[2], label='2.1')
leg22 = mpatches.Patch(color=palette[3], label='2.2')
leg3 = mpatches.Patch(color=palette[4], label='3.1')
plt.legend(handles=[leg1A, leg1B, leg21, leg22, leg3], prop={'size': 14}, loc='best', title='Family', title_fontsize=16)


### save outputs
plt.savefig('outputs/CoOptionCounterBarPlot.png', dpi=500, format='png')