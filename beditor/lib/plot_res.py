import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

def plot_nt_composition(dintseqguidesflt,pos_muts,plotp=None):
    fig=plt.figure(figsize=[16,10])
    for mi,method in enumerate(pos_muts.index):
        pos_min=pos_muts.loc[method,'Position of mutation from PAM: minimum']
        pos_max=pos_muts.loc[method,'Position of mutation from PAM: maximum']
        poss=np.arange(pos_min,pos_max+1)
        for posi,pos in enumerate(poss):
            seqs=list(dintseqguidesflt.loc[:,[c for c in dintseqguidesflt if ('guide sequence+PAM' in c) \
                                              and ('mutation position={}'.format(pos) in c)\
                                             and (' {}: '.format(method) in c)]].unstack().dropna())
            if len(seqs)>2:
                # [c for c in dintseqguidesflt if 'dintseqguidesflt codon position=-16' in c]
                instances=[Seq.Seq(s) for s in seqs]
                motif=motifs.create(instances)
                d=pd.DataFrame(motif.counts)
                d.index=d.index-20
                d.index.name='Position'
                ax=plt.subplot(5,len(pos_muts.index),posi*2+1+mi) #012 135 
                ax=d.plot(grid=True,xticks=d.index,ax=ax,lw=3)
                ax.set_ylabel('count')
                ax.axvspan(0,2,lw=4,color='lime',zorder=-1, label='PAM')    
                ax.axvline(pos,ax.get_ylim()[0],ax.get_ylim()[1],lw=3,color='k',zorder=-1, label='Mutation at')
                ax.axvspan(pos_min,pos_max,lw=4,color='lightgray',zorder=-1,label='Activity window')
                ax.legend(loc=2,bbox_to_anchor=[1,1])
    plt.suptitle(list(pos_muts.index))
    plt.subplots_adjust(wspace = 0.5)
            #     break
#     plt.tight_layout()
    if not plotp is None:
        plt.savefig(plotp)