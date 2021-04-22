import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import networkx as nx

### set parameters for identifying nodes
Family = 'Fam1B.1'  ### family to plot
MaxSeqs = 6     ### maximum number of sequences for each category in SortList (some sequences may be redundant)
MaxEdges = 1954    ### maximum number of edges for a single node
rt = 30        ### all_sum activity threshold
SortList = ['r_BWO','r_BFO','r_BLO','r_BIO','r_BVO','r_BMO']    ### parameters to include in identifying relevant sequences
ascend = False  ### if False, will sort parameters in descending order
RemoveIsos = True   ### if True, remove isolated nodes from final plot


### read data
df1 = pd.read_csv('data/WFLIVM-k-seq_merged_+r+I.csv')
df2 = df1[df1['Family'] == Family].sort_values(by=Family)
df3 = df2[df2['all_sum']>rt]
AllSeqs = df3['seq'].to_list()

### define function for identifying hamming distances from one sequence to a list of sequences
### inputs: a reference sequence string and a list of sequence strings
### outputs: a list of hamming distances (substitutions) from reference sequence to each sequence in list of sequences
def hamming(refseq,seqlist):
    HamList = []
    for seq in seqlist:
        ham = 0
        zipped = zip(refseq, seq)
        for i,j in zipped:
            if i != j:
                ham +=1
        HamList.append(ham)
    return(HamList)

### define function for identifying all intervening single substitution mutants between two double-substitution mutants
### inputs: two sequence strings
### outputs: a list of all direct single mutants (strings) between two double mutants
def FindMuts(start,end):
    zipped = list(zip(start,end))
    MutList = []
    pos=0
    MutPos=[]
    for i,j in zipped:
        if i!=j:
            MutPos.append(pos)
        pos+=1
    pos=0
    for m in range(len(MutPos)):
        mutant = list(start)
        mutant[MutPos[pos]] = zipped[MutPos[pos]][1]
        mutant = ''.join(mutant)
        MutList.append(mutant)
        pos+=1
    return(MutList)


### makes list of all sequences (strings) containing the top MaxSeqs values in SortList categories
### prints wildtype sequence and sequence with highest value for all items in SortList
SortedSeqs = [AllSeqs[0]]
print('wt: ',AllSeqs[0])
for i in range(MaxSeqs):
    for sub in SortList:
        seq = df3.sort_values(by=sub, ascending=ascend).iloc[i]
        SortedSeqs.append(seq['seq'])
        if i == 0:
            print(sub + ': ',seq['seq'])
SeqList = list(set(SortedSeqs))


### adds all intervening direct single mutant sequences between double mutant sequences in SeqList
### removes all sequences absent from AllSeqs list (ie intervening sequences absent in dataset)
n=0
AddSeqs = []
for seq in SeqList:
    ham = hamming(seq, SeqList)
    if 2 not in ham:
        print(seq + ' has no single or double mutants in the list')
    n+=1
    s=0
    for h in ham:
        if h == 2:
            muts = FindMuts(seq,SeqList[s])
            for i in muts:
                AddSeqs.append(i)
        s+=1

AddSeqs = list(set(AddSeqs))
SeqList = SeqList + AddSeqs
SeqList = list(set(SeqList))

DelSeqs = []
for seq in SeqList:
    if seq not in AllSeqs:
        DelSeqs.append(seq)
for seq in DelSeqs:
    SeqList.remove(seq)


### due to the low activity of many Family 3.1 sequences, the article figure combines these sequence lists: the main peak and orphan sequences
if Family == 'Fam3.1':
    Fam31Seqs = ['AAGTCCGCTAATAGTCGCAAG', 'AGGTTTGCTAATAGTCGCAAG', 'AGGTCTGCTAATAGTCGCAAG', 'GAGTTTGCTAATAGTCGCAAG', 'AAGTTTGCTAATAGTCGTAAG', 'AAGTCTGCTAATAGTCGTAAG', 'AAGCCTGCTAATAGTCGCAAG', 'GAGTCTGCTAATAGTCGCAAG', 'AAGTTTGCTAATAGTCGCAAG', 'TAGTCTGCTAATAGTCGCAAG', 'AAGTCTGCTAATAGTCGCAAG']
    Fam31orphs = ['AACTTCGCTAATAGTCGCAAG', 'AAGCTTGCTAATAGTCGCGAG', 'AACTTTGCTAAGAGTCGCAAG', 'AACTTTGCTAATAGTCGCGAG', 'AAGTTTGCTAATAGTCGCGCG', 'AAGTTTGCTAATAGTCGCGGG', 'AAGTTTGCTAATAGTCGCCGG']
    SeqList = Fam31Seqs


### make network plot
G = nx.Graph()

### add edges to nodes with hamming distance = 1
spos = 0
for seq in SeqList:
    file = 'outputs/pies/' + Family + '/' + seq + '.png'
    img = mpimg.imread(file)
    dists = hamming(seq,SeqList)
    G.add_node(spos,image=img)
    mpos = 0
    for h in dists:
        if h == 1:
            G.add_edge(spos,mpos)
        mpos += 1
        if mpos == MaxEdges:
            break
    spos += 1

### remove isolated nodes
if len(list(nx.isolates(G))) != 0:
    print('List of isolated nodes: ',list(nx.isolates(G)))
    if RemoveIsos == True:
        G.remove_nodes_from(list(nx.isolates(G)))
        print('Isolated nodes removed from final plot')

### define node layout
pos = nx.fruchterman_reingold_layout(G, k=1)
pos = nx.kamada_kawai_layout(G)

### plot network and map images to nodes
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
nx.draw_networkx_edges(G,pos,ax=ax)

plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)

trans = ax.transData.transform
trans2 = fig.transFigure.inverted().transform

piesize=0.05    ### image size
p2=piesize/2.0
for n in G:
    xx,yy=trans(pos[n])
    xa,ya=trans2((xx,yy))
    a = plt.axes([xa-p2,ya-p2, piesize, piesize])
    a.set_aspect('equal')
    a.imshow(G.nodes[n]['image'])
    a.axis('off')
ax.axis('off')


### make table of sequences in plot with SortList values and X,Y coordinates from plot
SeqDF = pd.DataFrame()
SeqDF['seq'] = SeqList
SeqDF = SeqDF.merge(df1[['seq']+SortList], on='seq')
PosDF = pd.DataFrame()
PosDF = PosDF.from_dict(pos, orient='index')
PosDF = PosDF.rename(mapper={0:'X',1:'Y'}, axis=1)
SeqDF = SeqDF.merge(PosDF, left_index=True, right_index=True)


### save output files
SeqDF.to_csv('outputs/pies/' + Family + '_PieNetwork_v.csv', index=False)
plt.savefig('outputs/pies/'  + Family + '_PieNetwork_v.png', dpi=500, transparent=True, format='png')