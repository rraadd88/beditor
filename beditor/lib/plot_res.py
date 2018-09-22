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
plt.style.use('ggplot')
from os.path import exists
from os import makedirs

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

import logging

from beditor.lib.io_strs import make_pathable_string
from beditor.lib.io_dfs import df2info
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from beditor.lib.io_dfs import df2info

def data2sub_matrix(data_fit,
                    values_col,
                    index_col,
                    col_ref='mut',
                    aggfunc='mean',test=False): 
    """
    This creates substitition matrix from input data (frequncies `data_lbl` or `data_fit`).
    
    :param data_fit: fitness data (`data_fit`).
    :param values_col: column name with values (str). 
    :param index_col: column name with index (str).
    :param type_form: type of values ["aas" : amino acid | "cds" : codons].
    """ 
    if aggfunc=='sum':
        data_fit['count']=1
    if aggfunc=="mean":
        data_sub_matrix=pd.pivot_table(data_fit,values=values_col,index=index_col,columns=col_ref)
    else:
        data_sub_matrix=pd.pivot_table(data_fit,values=values_col,index=index_col,columns=col_ref,aggfunc=aggfunc)            
    return data_sub_matrix
    
def plot_submap_possibilities(dmutagenesis,plotpf,test=False):
    import seaborn as sns
#     dmutagenesis=d.dropna(how='all')
    dmutagenesis=dmutagenesis.dropna(how='all')
    muttype2c={'All mutations':0,
               'Single nucleotide mutations':1,
               'Double nucleotide mutations':2,
               'Triple nucleotide mutations':3}

    aa2grp=pd.DataFrame({'aa':['A','G','I','L','P','V','M','C','N','Q','S','T','H','K','R','D','E','F','W','Y','*'],
    'grp':['non polar','non polar','non polar','non polar','non polar','non polar','neutral','neutral polar','neutral polar','neutral polar','neutral polar','neutral polar','positive','positive','positive','negative','negative','aromatic','aromatic','aromatic','*'],
    }).set_index('aa')

    grp2clr=dict(zip(['non polar', 'neutral', 'neutral polar', 'positive', 'negative','aromatic', '*'],
                    ['r', 'g', 'b', 'orange', 'm','c', 'k']))
    aa2grp['c']=[grp2clr[g] for g in aa2grp['grp']]

    aa2grp['rank']=range(len(aa2grp))
    aa2grp.index.name='amino acid'
    
    from Bio.Alphabet import IUPAC
    aas=IUPAC.IUPACData.protein_letters+'*'
    from beditor.lib.get_mutations import get_codon_table
    dcodontable=get_codon_table(aa=list(aas))
    
    cd2grp=dcodontable.set_index('amino acid').join(aa2grp).sort_values(by='rank').reset_index().set_index('codon')    
#     df2info(dcodontable)
    methods=['All']+list(dmutagenesis['method'].unique())
    for muttype1 in ['amino acid','codon']:
        if muttype1=='codon':
            figsize=[60,60]
        else:
            figsize=[30,30]        
        fig,axes=plt.subplots(figsize=figsize,nrows=4, ncols=len(methods)+1,)
        for methodi,method in enumerate(methods):
            for muttypei,muttype in enumerate(muttype2c.keys()):
                if not 'All' in method:
                    dplot=dmutagenesis.loc[(dmutagenesis['method']==method),:]
                else:
                    dplot=dmutagenesis.copy()
                    
                if not 'All' in muttype:
                    dplot=dplot.loc[(dplot['nucleotide mutation: count']==muttype2c[muttype]),:]
                ax=axes[muttypei][methodi]
                if len(dplot)!=0:
#                     print(method,muttype)
#                     df2info(dplot)
                    dplot=data2sub_matrix(dplot,values_col='count',col_ref=muttype1,
                                                   index_col=muttype1+' mutation',
                                                   aggfunc='sum',
#                                                   cmap=,
                                                 )                    
                    if muttype1=='amino acid':
                        annot=False
                        muttype12grp=aa2grp
                        rotation=90
                    else:
                        annot=False
                        muttype12grp=cd2grp
                        rotation=90

                    dplot=dplot.loc[muttype12grp.index,muttype12grp.index]
                    dplot[pd.notna(dplot)]=1
                    ax=sns.heatmap(dplot,ax=ax,#fmt="d",
                                  vmin=0, vmax=1,
                                  cbar=False,
                                  annot=annot,
                                  square=True,
                                  linewidths=0.1, linecolor='gray',
                                  cmap='YlOrRd')
                    ax.set_ylabel('Mutation to',labelpad=20,color='gray')
                    ax.set_xlabel('Wild type',labelpad=20,color='gray')
                    for ticklabely,aa in enumerate(muttype12grp.index.tolist()):
                        if not test:
                            ax.set_yticklabels([])
                            ticklabelx=-1
                        else:
                            ticklabelx=-2
                        ax.text(ticklabelx+0.5,ticklabely+0.75,aa,color=muttype12grp.loc[aa,'c'],fontsize=10,
                            ha='right',)
                        if test:
                            ticklabelx+=2
                        else:
                            ax.set_xticklabels([])
                        ax.text(ticklabely+0.5,ticklabelx+len(muttype12grp.index)+1.5,aa,color=muttype12grp.loc[aa,'c'],fontsize=10,
                            ha='center',va='top',rotation=rotation)
                    ax.set_title('method={}\n{}'.format(method,muttype))
        #                 plt.savefig('plots/heatmap_{}_{}_{}.svg'.format(muttype1.replace(' ','_'),method,muttype))
        #                 break
        #             break
        # legend
                else:
                    fig.delaxes(ax)
        for muttypei,muttype in enumerate(muttype2c.keys()):
            ax=axes[muttypei][len(methods)]
            if muttypei==0:
                for aa in aa2grp.index:
                #     ax.scatter(1,cd2grp.loc[aa,'rank'],color=cd2grp.loc[aa,'c'])
                    ax.scatter(1,aa2grp.loc[aa,'rank'],color=aa2grp.loc[aa,'c'])
                    ax.text(0.98,aa2grp.loc[aa,'rank']+0.15,aa,color=aa2grp.loc[aa,'c'],
                           ha='center',va='center')        
                    ax.text(1.01,aa2grp.loc[aa,'rank']+0.15,cd2grp.loc[(cd2grp['amino acid']==aa),:].index.tolist(),
                            color=aa2grp.loc[aa,'c'],
                           ha='left',va='center')    
                for c in aa2grp['c'].unique():
                    ax.text(0.96,aa2grp.loc[aa2grp['c']==c,:]['rank'].mean()+0.15,
                            aa2grp.loc[aa2grp['c']==c,:]['grp'].unique()[0],
                            color=c,
                           ha='right',va='center')    
                ax.invert_yaxis()
                ax.set_xlim([0.8,1.5])
                ax.set_facecolor('w')
                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)
            #                 plt.axis('off')
            else:
                fig.delaxes(ax)
        if not '{mutation_type}' in plotpf:
            plotpf=plotpf+'_{mutation_type}.png'
        plt.savefig(plotpf.format(mutation_type=muttype1.replace(' ','_')))
    #     break

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
        df.loc[:,'sense_']=0        
        #pam
        df.loc[:,'PAM feature name']=df.apply(lambda x : f"{x['PAM']}",axis=1)
        features_pam=df2features(df.loc[:,['position of PAM ini','position of PAM end','PAM feature name','sense_']].drop_duplicates())
        features+=features_pam
    if 'guide' in types:
        #guide
        df.loc[:,'sense_']=0
        features_guide=df2features(df.loc[:,['position of guide ini','position of guide end','strategy','sense_']].drop_duplicates())
        features+=features_guide
    if 'guide+pam' in types:
        df.loc[:,'sense_']=0
        #guide+pam
        df.loc[:,'guide+pam ini']=df.apply(lambda x : np.min([x['position of PAM ini'],x['position of guide ini']]),axis=1)
        df.loc[:,'guide+pam end']=df.apply(lambda x : np.max([x['position of PAM end'],x['position of guide end']]),axis=1)   
        features_guidepam=df2features(df.loc[:,['guide+pam ini','guide+pam end','strategy','sense_']].drop_duplicates())
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
        features.append(SeqFeature(FeatureLocation(start=int(df.loc[name,colini]),
                                                   end=int(df.loc[name,colend])+1,
                                                   strand=int(df.loc[name,colsense]),), 
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
    graphic_record.plot_sequence(ax=ax,)
    graphic_record.plot_translation(ax=ax,location=[0,45])    
    ax.plot([annot_residuei*3-3.5,annot_residuei*3-0.5],[-2,-2],lw=5,color='r')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
#     ax.plot([21,23],[-2,-2])
    if not plotp is None:
        plt.tight_layout()
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
    ax.axvspan(dp['Activity window'].min(),dp['Activity window'].max(),lw=4,color='khaki',zorder=-1,label='Activity window')
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
        pamposs=dp.index.tolist()[-1*paml:]
        windowmin=dguideslin_sub123['distance of mutation from PAM: minimum'].unique()[0]*-1
        windowmax=dguideslin_sub123['distance of mutation from PAM: maximum'].unique()[0]*-1        
    else:
        dp.index=dp.index-paml
        pamposs=dp.index.tolist()[:paml]        
        windowmin=dguideslin_sub123['distance of mutation from PAM: minimum'].unique()[0]
        windowmax=dguideslin_sub123['distance of mutation from PAM: maximum'].unique()[0]        

    dp.loc[dguideslin_sub123['distance of mutation from PAM: minimum'].unique()[0],'Activity window']=windowmin
    dp.loc[dguideslin_sub123['distance of mutation from PAM: maximum'].unique()[0],'Activity window']=windowmax
    dp.loc[pos,'Mutation at']=pos
    for pamnti,pamnt in enumerate(list(pam)):
        dp.loc[pamposs[pamnti],'PAM']=pamposs[pamnti]
        dp.loc[pamposs[pamnti],'PAM nt']=pamnt
    return dp

def plot_bar_dguides(dstep,plotp,figsize=None,figfmt='png',):
    cols=['method','PAM','strand']
    fig_ht=np.max([len(dstep[c].unique()) for c in cols])
    if figsize is None:
        figsize=[5,3+fig_ht*0.2]
    fig,axes=plt.subplots(nrows=len(cols),figsize=figsize,sharex=True)
    for i,col in enumerate(cols):
        dplot=pd.DataFrame(dstep[col].value_counts()).T
        if dplot.shape[0]!=0:
            if dplot.shape[1]!=0:
                dplot.plot.barh(ax=axes[i],stacked=True)
                axes[i].legend(bbox_to_anchor=[1,1])
                if i==len(cols)-1:
                    axes[i].set_xlabel('# of guides designed')            
    plt.tight_layout()
    plt.savefig(plotp,format=figfmt)
    
def plot_dist_dguides(dguideslin,dpam,plotpf=None):
    #method fig
    dps=[]
    for met in dguideslin['method'].unique():
        dps_met={}
        dguideslin_sub1=dguideslin.loc[((dguideslin['method']==met)),:]
        #pam
        dps_pam={}
        for pam in np.sort(dguideslin_sub1['PAM'].unique()):
            dguideslin_sub12=dguideslin_sub1.loc[(dguideslin_sub1['PAM']==pam),:]
            #position
            dps_pos={}
            poss=dguideslin_sub12['position of mutation from PAM'].unique()
            for pos in np.sort(poss):
                dguideslin_sub123=dguideslin_sub12.loc[(dguideslin_sub12['position of mutation from PAM']==pos),:]
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
                              figsize=[8.25*(pamc),2*(posc)],
    #                          sharey=True,
                             )
        for pami,pam in enumerate(dps_met[met].keys()):
            for posi,pos in enumerate(dps_met[met][pam].keys()):
                if len(np.shape(axes))==2:
                    ax=axes[posi,pami]
                else:
                    try:
                        ax=axes[posi]
                    except:
                        ax=axes
                plot_ntcompos(dp=dps_met[met][pam][pos],
                      legend=False,title=f"PAM={pam}; position of mutation={pos}",
                      ax=ax,fig=fig)
                if posi==0 and pami==len(dps_met[met].keys())-1:
                    ax.legend(loc=2,bbox_to_anchor=[1,1.15])

    #             break
    #         break
        plt.tight_layout()
        dps.append(dps_met) 
        if not plotpf is None:
            if not '{method}' in plotpf:
                plotpf=plotpf+'_{method}.png'
            plt.savefig(plotpf.format(method=met))
    return dps
#         break                      
def plot_dna_features_view(cfg,dsequences,dguides,plotd,more=False,guideids=None):
    if guideids is None:
        guideids=dguides['id'].unique()
        if not more:
            np.random.seed(seed=8888)            
    guidesc=len(guideids)
    if guidesc!=0:
        makedirs(plotd,exist_ok=True)
        ids=[guideids[i] for i in np.sort(np.random.choice(range(guidesc), 10 if guidesc>10 else guidesc-1, replace=False))]
        for id in ids:
            gbp=f"{plotd}/{make_pathable_string(id)}.gb"
            plotp=f"{gbp}.png"        
            if id in dsequences['id']:
                seq_transcript=dsequences.loc[(dsequences['id']==id),'transcript: sequence'].tolist()[0]
            else:
                logging.info('mutation_format=nucleotide detected.')
                if cfg['mutation_format']=='nucleotide':
                    seq_transcript=dsequences.loc[(dsequences['id']==id.split('|')[0]),'transcript: sequence'].tolist()[0] 
                elif cfg['mutation_format']=='aminoacid':
                    seq_transcript=dsequences.loc[(dsequences['id']==id),'transcript: sequence'].tolist()[0] 
            record=make_gb(sequence_string=seq_transcript,
                    features=get_df4features(dguides,id,types=['guide+pam']),gbp=gbp,
                       )
            plot_seq(record,
                     annot_residuei=8,
                     xlabel=f'Sequence flanking target site\n bed id: {id}',
                     plotp=plotp
                    )

def plot_dist_dofftargets(dofftargets,plotp=None):
    dofftargets=dofftargets.replace([np.inf, -np.inf], np.nan)
    plt.figure(figsize=[3,3])
    ax=plt.subplot(111)
    _=dofftargets['beditor score (log10)'].hist(bins=100,ax=ax)
    ax.set_xlim(-50,0)
    ax.set_xlabel('beditor score (log10 scale)')
    ax.set_ylabel('number of guides')
    plt.tight_layout()
    plt.savefig(plotp)
    
from glob import glob
from beditor.lib.io_dfs import del_Unnamed,set_index
from os.path import exists,splitext,dirname,splitext,basename,realpath
def plot_vizbysteps(cfg):  
    prjd=cfg['prjd']
    #make one output table and stepwise plots
    datad=f"{prjd}/05_output"
                               
    # step2 # make submap
    stepi=2
    plotp=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_substitution_map"
    plotps=glob(plotp+'*')
    if len(plotps)==0 or cfg['force']:
        plotpf=plotp+"_{mutation_type}.png"
        dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
        if exists(dstepp):
            dstep=del_Unnamed(pd.read_table(dstepp)).drop_duplicates()
            if len(dstep)<1000:
                logging.info('plot_submap_possibilities')
                plot_submap_possibilities(dmutagenesis=dstep,
                                          plotpf=plotpf,test=False)
            else:
                logging.warning(f'skipped: plot_submap_possibilities')
        else:
            logging.warning(f'not found: {dstepp}')

    # step3 
    # stats by strategies
    stepi=3
    plotp=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_stats_by_strategies.png"
    if not exists(plotp) or cfg['force']:                               
        dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
        if exists(dstepp):
            dstep=del_Unnamed(pd.read_table(dstepp)).drop_duplicates()
            logging.info('plot_bar_dguides')
            plot_bar_dguides(dstep,plotp)
        else:
            logging.warning(f'not found: {dstepp}')

    # make nt_composition plot
    stepi=3
    plotp=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_nt_compositions"
    plotps=glob(plotp+'*')
    if len(plotps)==0 or cfg['force']:
        plotpf=plotp+"_{method}.png"
        dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
        if exists(dstepp):
            dstep=del_Unnamed(pd.read_table(dstepp)).drop_duplicates()
            dpam=pd.read_table('{}/../data/dpam.tsv'.format(dirname(realpath(__file__))))
            dpam=set_index(dpam,'PAM')
            logging.info('plot_dist_dguides')
            plot_dist_dguides(dstep,dpam,plotpf)
        else:
            logging.warning(f'not found: {dstepp}')

    # make plot_dna_features_view
    stepi=3
    plotd=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_dna_features_view"
    plotps=glob(plotd+'/*')
    if len(plotps)==0 or cfg['force']:
        dguidesp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
        dsequencesp=f"{cfg[stepi-2]}/d{cfg[stepi-2].replace('/','').split('_')[-1]}.tsv"
        if exists(dguidesp):
            logging.info('plot_dna_features_view')
            plot_dna_features_view(cfg,dsequences=del_Unnamed(pd.read_table(dsequencesp)).drop_duplicates(),
                       dguides=del_Unnamed(pd.read_table(dguidesp)).drop_duplicates(),
                       plotd=plotd,more=False)
        else:
            logging.warning(f'not found: {dstepp}')
            
#     # step2 # make submap #FIXME get all the columns used for plotting in the dguides.
#     stepi=3
#     plotp=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_submap_used_for_mutagenesis"
#     plotps=glob(plotp+'*')
#     if len(plotps)==0 or cfg['force']:
#         plotpf=plotp+"_{mutation_type}.png"
#         dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
#         dstep=del_Unnamed(pd.read_table(dstepp)).drop_duplicates()
#         logging.info('plot_submap_possibilities')
#         plot_submap_possibilities(dmutagenesis=dstep,
#                                   plotpf=plotpf,test=False)

    # step4 offtargets correlations  
    stepi=4
    plotp=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_dist_beditor_score.png"
    if not exists(plotp) or cfg['force']:
        dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
        if exists(dstepp):
            dstep=del_Unnamed(pd.read_table(dstepp)).drop_duplicates()
            logging.info('plot_dist_dofftargets')
            plot_dist_dofftargets(dstep,plotp)
        else:
            logging.warning(f'not found: {dstepp}')

#     # step5
#     stepi=4
#     plotp=f"{datad}/plot_d{cfg[stepi].replace('/','').split('_')[-1]}_dist_beditor_score.png"
#     if not exists(plotp) or cfg['force']:                               
#         dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
#         dstep=del_Unnamed(pd.read_table(dstepp)).drop_duplicates()
#         logging.info('plot_dist_dofftargets')
#         plot_dist_dofftargets(dstep,plotp)                
