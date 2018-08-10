def is_interactive():
    # thanks to https://stackoverflow.com/a/22424821/3521099
    import __main__ as main
    return not hasattr(main, '__file__')

if not is_interactive():
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

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

def get_nt_composition(seqs):
    instances=[Seq.Seq(s) for s in seqs]
    motif=motifs.create(instances)
    d=pd.DataFrame(motif.counts)
    d.index=range(1,d.shape[0]+1)
    return d

def plot_ntcompos(dp,
                  legend=True,title=None,
                  ax=None,fig=None):
    if fig is None:
        fig=plt.figure(figsize=[15,4])
    if ax is None:
        ax=plt.subplot(111)
    ax=dp.loc[:,list('ATGC')].plot(grid=True,xticks=dp.index,ax=ax,lw=3,
           color=['r','b','orange','g'],
           label=None,
    #                cmap='Spectral'
           legend=False,
          )

    ax.set_ylabel('count')
    #pam
    ax.axvspan(dp['PAM'].min(),dp['PAM'].max(),lw=4,color='lime',zorder=-1, label='PAM')    
    #mutation
    ax.axvline(dp['Mutation at'].min(),ax.get_ylim()[0],ax.get_ylim()[1],lw=3,color='k',zorder=-1, label='Mutation at')
    #activity window
    ax.axvspan(dp['Activity window'].min(),dp['Activity window'].max(),lw=4,color='lightgray',zorder=-1,label='Activity window')
    for pampos in dp['PAM'].dropna():
        ax.text(pampos,0,dp.loc[pampos,'PAM nt'],
                ha='center',va='center')
    if not title is None:
        ax.set_title(title)
    if legend:
        ax.legend(loc=2,bbox_to_anchor=[1,1])
    return ax

def get_dntcompos(dguideslin_sub123,dpam,pos,pam):
    dp=get_nt_composition(dguideslin_sub123['guide+PAM sequence'])
    paml=len(dguideslin_sub123['PAM'].unique()[0])

    if dpam.loc[pam,'position']=='down':
        dp.index=(dp.index[::-1]-paml)*-1
        pos=pos
        pamposs=dp.index.tolist()[-1*paml:]
        windowmin=dguideslin_sub123['Position of mutation from PAM: minimum'].unique()[0]
        windowmax=dguideslin_sub123['Position of mutation from PAM: maximum'].unique()[0]        
    else:
        dp.index=dp.index-paml
        pos=pos
        pamposs=dp.index.tolist()[:paml]        
        windowmin=dguideslin_sub123['Position of mutation from PAM: minimum'].unique()[0]*-1
        windowmax=dguideslin_sub123['Position of mutation from PAM: maximum'].unique()[0]*-1        

    dp.loc[dguideslin_sub123['Position of mutation from PAM: minimum'].unique()[0],'Activity window']=windowmin
    dp.loc[dguideslin_sub123['Position of mutation from PAM: maximum'].unique()[0],'Activity window']=windowmax
    dp.loc[pos,'Mutation at']=pos
    for pamnti,pamnt in enumerate(list(pam)):
        dp.loc[pamposs[pamnti],'PAM']=pamposs[pamnti]
        dp.loc[pamposs[pamnti],'PAM nt']=pamnt
    return dp

def plot_bar_dguides(dstep,plotp):
    cols=['method','PAM','strand']
    fig,axes=plt.subplots(nrows=len(cols),figsize=[3,3],sharex=True)
    for i,col in enumerate(cols):
        dstep[col].value_counts().plot.barh(ax=axes[i])
        axes[i].set_ylabel(col)
        if i==len(cols)-1:
            axes[i].set_xlabel('count of guides designed')
    plt.tight_layout()
    plt.savefig(plotp)

def plot_dist_dguides(dguideslin,dpam,plotpf):
    #method fig
    for met in dguideslin['method'].unique():
        dps_met={}
        dguideslin_sub1=dguideslin.loc[((dguideslin['method']==met)),:]
        #pam
        dps_pam={}
        for pam in np.sort(dguideslin_sub1['PAM'].unique()):
            dguideslin_sub12=dguideslin_sub1.loc[(dguideslin_sub1['PAM']==pam),:]
            #position
            dps_pos={}
            for pos in np.sort(dguideslin_sub12['distance of mutation in codon from PAM'].unique()):
                dguideslin_sub123=dguideslin_sub12.loc[(dguideslin_sub12['distance of mutation in codon from PAM']==pos),:]
                #combine strands
                dps_pos[pos]=get_dntcompos(dguideslin_sub123,dpam,pos,pam)
            dps_pam[pam]=dps_pos
    #             break
    #         break
        dps_met[met]=dps_pam
        pamc=np.max([len(dps_met[met]) for pam in dps_met])
        posc=np.max([len(dps_met[met][pam]) for pam in dps_met[met]])
        fig,axes=plt.subplots(ncols=pamc,
                              nrows=posc,
                              figsize=[8*(pamc),2*(posc)],
    #                          sharey=True,
                             )
    #     for met in dps_met.keys():
        for pami,pam in enumerate(dps_met[met].keys()):
            for posi,pos in enumerate(dps_met[met][pam].keys()):
                ax=axes[posi,pami]
                plot_ntcompos(dp=dps_met[met][pam][pos],
                      legend=False,title=f"PAM={pam}; position of mutation={pos}",
                      ax=ax,fig=fig)
                if posi==0 and pami==len(dps_met[met].keys())-1:
                    ax.legend(loc=2,bbox_to_anchor=[1,1.15])

    #             break
    #         break
        plt.tight_layout()
        plt.savefig(plotpf.format(method=met))
#         break                      

def plot_dist_dofftargets(dofftargets,plotp)
    dofftargets=dofftargets.replace([np.inf, -np.inf], np.nan)
    plt.figure(figsize=[3,3])
    ax=plt.subplot(111)
    _=dofftargets['beditor score (log10)'].hist(bins=100,ax=ax)
    ax.set_xlabel('beditor score (log10 scale)')
    ax.set_ylabel('number of guides')
    plt.tight_layout()
    plt.savefig(plotp)