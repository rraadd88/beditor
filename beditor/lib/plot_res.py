import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

def plot_nt_composition(dintseqguidesflt,pos_min,pos_max,plotp=None):
    poss=np.arange(pos_min,pos_max+1)
    plt.figure(figsize=[8,2*len(poss)])
    for posi,pos in enumerate(poss):
        seqs=list(dintseqguidesflt.loc[:,[c for c in dintseqguidesflt if ('guide sequence+PAM' in c) and ('mutation position={}'.format(pos) in c)]].unstack().dropna())
        # [c for c in dintseqguidesflt if 'dintseqguidesflt codon position=-16' in c]
        instances=[Seq.Seq(s) for s in seqs]
        motif=motifs.create(instances)
        d=pd.DataFrame(motif.counts)
        d.index=d.index-20
        d.index.name='Position'
        ax=plt.subplot(len(poss),1,posi+1)
        ax=d.plot(grid=True,xticks=d.index,ax=ax,lw=2)
        ax.set_ylabel('count')
        ax.axvspan(0,2,lw=4,color='lime',zorder=-1, label='PAM')    
        ax.axvline(pos,ax.get_ylim()[0],ax.get_ylim()[1],lw=3,color='k',zorder=-1, label='Mutation at')
        ax.axvspan(pos_min,pos_max,lw=4,color='lightgray',zorder=-1,label='Activity window')
        ax.legend(loc=2,bbox_to_anchor=[1,1])
    #     break
    plt.tight_layout()
    if not plotp is None:
        plt.savefig(plotp)
