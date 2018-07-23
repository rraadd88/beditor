import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

def get_df4features(dguideslin,id,types=['guide+pam']):
    features=[]
    df=dguideslin.loc[(dguideslin['id']==id),:]
    df.loc[:,'sense']=df.apply(lambda x : f"{x['strand']}1",axis=1)
    if 'codon' in types:
        #codon
        df=pd.DataFrame({'start':[21],
                        'end':[23],
                        'codon':['codon'],
                        'sense':['0']})
        feature_codon=df2features(df)
        features+=feature_codon
    if 'pam' in types:
        #pam
        df.loc[:,'PAM feature name']=''#df.apply(lambda x : f"{x['strand']}PAM",axis=1)
        features_pam=df2features(df.loc[:,['position of PAM ini','position of PAM end','PAM feature name','sense']].drop_duplicates())
        features+=features_pam
    if 'guide' in types:
        #guide
        df.loc[:,'sense_']=0
        # df.loc[:,'guide ini']=df.apply(lambda x : x['position of PAM end'] if x['position']=='up' else x['position of PAM ini']-x['guide sequence length'],axis=1)
        # df.loc[:,'guide end']=df.apply(lambda x : x['position of PAM ini']-1 if x['position']=='down' else x['position of PAM ini']+x['guide sequence length'],axis=1)
        df.loc[:,'guide ini']=df.apply(lambda x : x['position of PAM end'] if x['position']=='up' else x['position of PAM end']-x['guide sequence length']-1,axis=1)
        df.loc[:,'guide end']=df.apply(lambda x : x['position of PAM ini']-1 if x['position']=='down' else x['position of PAM ini']+x['guide sequence length'],axis=1)
        features_guide=df2features(df.loc[:,['guide ini','guide end','strategy','sense_']].drop_duplicates())
        features+=features_guide
    if 'guide+pam' in types:
        #guide+pam
        df.loc[:,'guide+pam ini']=df.apply(lambda x : x['position of PAM ini'] if x['position']=='up' else x['position of PAM ini']-x['guide sequence length']-1,axis=1)
        df.loc[:,'guide+pam end']=df.apply(lambda x : x['position of PAM end'] if x['position']=='down' else x['position of PAM end']+x['guide sequence length'],axis=1)
        features_guidepam=df2features(df.loc[:,['guide+pam ini','guide+pam end','strategy','sense']].drop_duplicates())
        features+=features_guidepam
    return features

def df2features(df):
    """
    cols= ini, end, name,sense
    """
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from beditor.lib.io_dfs import set_index
    colini,colend,colname,colsense=df.columns
    df=set_index(df,colname)
    features=[]
    df=df.reset_index()
    for name in df.index:
#         print(int(df.loc[name,colini]),int(df.loc[name,colend]))
        features.append(SeqFeature(FeatureLocation(start=int(df.loc[name,colini]), 
                                                   end=int(df.loc[name,colend])+1,
                                                   strand=int(df.loc[name,colsense]),                                                   
                                                  ), 
                                   type=df.loc[name,colname],
                                   
                                  ))
    return features

def make_gb(sequence_string,features=[],
            seqid='seqid',seqname='seqname',seqdesc='seqdesc',
            gbp=None,
           ):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    # thanks to https://caretdashcaret.com/2016/06/15/how-to-create-genbank-files-with-biopython/

    # Create a sequence
#     sequence_string = "ggggaaaattttaaaaccccaaaa"
    sequence_object = Seq(sequence_string, IUPAC.unambiguous_dna)

    # Create a record
    record = SeqRecord(sequence_object,
                       id=seqid, # random accession number
                       name=seqname,
                       description=seqdesc)

    # Add annotation
    for feature in features: 
        record.features.append(feature)

    if gbp is not None:
        # Save as GenBank file
        with open(gbp, 'w') as output_file:
            SeqIO.write(record, output_file, 'genbank')
    return record

def plot_seq(record,annot_residuei=8,
             title='',
             xlabel='',
             plotp=None):
    from dna_features_viewer import BiopythonTranslator
    # graphic_record = BiopythonTranslator().translate_record("seqname.gb")
    graphic_record = BiopythonTranslator().translate_record(record)
    ax, _ = graphic_record.plot(figure_width=12.5,
                               annotate_inline=True, 
                                level_offset=0.5,
                               )
    graphic_record.plot_sequence(ax=ax)
    graphic_record.plot_translation(ax=ax,location=[0,45])    
    ax.plot([annot_residuei*3-3.5,annot_residuei*3-0.5],[-2,-2],lw=5,color='r')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
#     ax.plot([21,23],[-2,-2])
    if not plotp is None:
        ax.figure.savefig(plotp,format='png')

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