import pandas as pd
import numpy as np
import json
import sys
from os.path import exists,realpath,dirname
from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO
import logging
from tqdm import trange
import pandas as pd
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO
import logging

from beditor.lib.io_dfs import *
from beditor.lib.global_vars import multint2reg,multint2regcomplement,nt2complement
from beditor.lib.io_seqs import reverse_complement_multintseq,reverse_complement_multintseqreg,str2seq
from beditor.lib.io_strs import s2re
from beditor.lib.global_vars import stepi2cols    


def get_pam_searches(dpam,seq,pos_codon,
                    test=False):
    """
    Search PAM occurance

    :param dpam: dataframe with PAM sequences
    :param seq: target sequence
    :param pos_codon: reading frame
    :param test: debug mode on
    :returns dpam_searches: dataframe with positions of pams
    """
    import regex as re
    def get_guide_pam(match,pam_stream,guidel,pos_codon):  
        if pam_stream=='down':
            # IIIIIIINGG
            # 0123456789
            seq_guidepam=seq[match.span()[0]-guidel:match.span()[1]]
            seq_guide=seq[match.span()[0]-guidel:match.span()[0]]
            seq_pam=seq[match.span()[0]:match.span()[1]]
            dist_codon=pos_codon-match.span()[0]
#             if test:
#                 print(match.span()[0]-pos_codon)
        elif pam_stream=='up':
            # TTNIIIIIII
            # 0123456789
            seq_guidepam=seq[match.span()[0]:match.span()[1]+guidel]
            seq_guide=seq[match.span()[1]:match.span()[1]+guidel]
            seq_pam=seq[match.span()[0]:match.span()[1]]
            dist_codon=pos_codon-match.span()[1]
        if seq_pam!=match.group():
            logging.error(f'indexing is wrong:seq_guidepam: {seq_guidepam}, seq_guide: {seq_guide}, seq_pam: {seq_pam},match.group(): {match.group()}')               
        return seq_guidepam,seq_guide,seq_pam,abs(dist_codon)
    pams=dpam.index.tolist()
    dpamposs=pd.DataFrame(columns=['guide+PAM sequence','guide sequence','PAM sequence'])
    pamposi=0
    for pam in pams:
        matchesiter=re.finditer(dpam.loc[pam,'rPAM'], seq, overlapped=True)
        for match in matchesiter:
            dpamposs.loc[pamposi,'position of PAM ini'],dpamposs.loc[pamposi,'position of PAM end']=match.span() 
            dpamposs.loc[pamposi,'position of PAM end']=dpamposs.loc[pamposi,'position of PAM end']-1                
            dpamposs.loc[pamposi,'guide+PAM sequence'],dpamposs.loc[pamposi,'guide sequence'],dpamposs.loc[pamposi,'PAM sequence'],dpamposs.loc[pamposi,'distance of codon from PAM']\
            =get_guide_pam(match,dpam.loc[pam,'position'],dpam.loc[pam,'guide length'],pos_codon)
            dpamposs.loc[pamposi,'PAM']=pam
            pamposi+=1
    dpamposs['codon: from pam search']=seq[pos_codon:pos_codon+3]
    dpamposs['guide sequence']=dpamposs['guide sequence'].fillna('')
    if len(dpamposs)==0:
        return None
    dpamposs['guide sequence length']=dpamposs.apply(lambda x : len(x['guide sequence']),axis=1)

    dpamposs['franking sequence']=seq
    dpamposs=dpam.join(set_index(dpamposs,'PAM'),how='right')
    dpamposs.index.name='PAM'
    return dpamposs

def guide2dpositions(x,dbug=False): 
    """
    Get positions of guides relative to the target site and PAM sequence
    Note:
    Index and flank sequence based indexing are 0-based
    Distances and positions from pam are 1-based

    :param x: lambda section of dguides dataframe  
    """
    dpositions=pd.DataFrame(index=range(45),
                           columns=['guide+PAM sequence'])
    dpositions.index.name='position'

    dpositions.loc[x['position of PAM ini']:x['position of PAM end'],'location PAM']=True
    if x['position']=='up':
        dpositions.loc[x['position of PAM end']+1:x['position of PAM end']+x['guide sequence length'],'location guide']=True
        dpositions.loc[x['position of PAM end']+x['distance of mutation from PAM: minimum']:x['position of PAM end']+x['distance of mutation from PAM: maximum'],'location window']=True        
    elif x['position']=='down':
        dpositions.loc[x['position of PAM ini']-x['guide sequence length']:x['position of PAM ini']-1,'location guide']=True
        dpositions.loc[x['position of PAM ini']-x['distance of mutation from PAM: maximum']:x['position of PAM ini']-x['distance of mutation from PAM: minimum'],'location window']=True        

    dpositions.loc[21:23,'location codon']=True
    dpositions.loc[21-1+x['position of mutation in codon'],'location mutation']=True
    dpositions[[c for c in dpositions if 'location' in c]]=dpositions[[c for c in dpositions if 'location' in c]].fillna(False)
    dpositions['location guide+PAM']=dpositions['location guide'] | dpositions['location PAM']
    if x['original position']=='down' and x['strand']=='+':
        dpositions['position from PAM']=np.array(dpositions.index.tolist())-dpositions.loc[dpositions['location PAM'],:].index.tolist()[0]
    elif x['original position']=='down' and  x['strand']=='-':
        dpositions['position from PAM']=np.array(dpositions.index.tolist())-dpositions.loc[dpositions['location PAM'],:].index.tolist()[-1]
        dpositions['position from PAM']=dpositions['position from PAM']*-1    
    elif x['original position']=='up' and x['strand']=='+':
        dpositions['position from PAM']=dpositions.loc[dpositions['location PAM'],:].index.tolist()[-1]-np.array(dpositions.index.tolist())
        dpositions['position from PAM']=dpositions['position from PAM']*-1    
    elif x['original position']=='up' and  x['strand']=='-':
        dpositions['position from PAM']=np.array(dpositions.index.tolist())-dpositions.loc[dpositions['location PAM'],:].index.tolist()[0]
        dpositions['position from PAM']=dpositions['position from PAM']*-1    

    dpositions.loc[dpositions['location mutation'],'nucleotide wild-type']=x['nucleotide']
    dpositions.loc[dpositions['location mutation'],'nucleotide mutation']=x['nucleotide mutation']

    dpositions.loc[dpositions['location codon'],'codon wild-type']=list(x['codon: wild-type'])
    dpositions.loc[dpositions['location codon'],'codon mutation']=list(x['codon mutation'])
    if x['strand']=='+':
        dpositions.loc[dpositions['location guide+PAM'],'guide+PAM sequence']=list(x['guide+PAM sequence'])
        activity_sequence=''.join(dpositions.loc[dpositions['location window'],'guide+PAM sequence'].tolist())
    elif x['strand']=='-':
        dpositions.loc[dpositions['location guide+PAM'],'guide+PAM sequence']=list(x['guide+PAM sequence'])[::-1]
        activity_sequence=''.join(dpositions.loc[dpositions['location window'],'guide+PAM sequence'].tolist())[::-1]
    posmut=dpositions.loc[dpositions['location mutation'],:].index[0]
    posmutfrompam=int(dpositions.loc[dpositions['location mutation'],'position from PAM'])
    distmutfrompam=abs(posmutfrompam)
    posguideini=dpositions.loc[dpositions['location guide'],:].index.min()
    posguideend=dpositions.loc[dpositions['location guide'],:].index.max()
    if dbug:
        print(x[['strategy','strand','distance of mutation from PAM']+[s for s in x.index if ('sequence' in s) or ('length' in s)]])
        print({'posmut':posmut,
               'posmutfrompam':posmutfrompam,
               'distmutfrompam':distmutfrompam,
               'posguideini':posguideini,
               'posguideend':posguideend,
              'activity_sequence':activity_sequence})
        return dpositions
    else:
        return posmut,posmutfrompam,distmutfrompam,posguideini,posguideend,activity_sequence
    
def make_guides(cfg,dseq,dmutagenesis,dpam,
               test=False,
               dbug=False):
    """
    Wrapper around submodules that design guides by
    1. searching all PAM sequences on 'both' the strands,
    2. filtering guides by all possible strategies (given in dmutagenesis) e.g. activity window,
    Finally generates a table.

    :param cfg: configuration dict
    :param dseq: dsequences dataframe
    :param dmutagenesis: dmutagenesis dataframe
    :param dpam: dpam dataframe
    :param test: debug mode on
    :param dbug: more verbose
    """
    from beditor.lib.io_strs import s2re
    flankaas=7#FIXME if flank length changes

    dseq=dseq.reset_index()
    dseq.index=range(len(dseq))
    dpam=set_index(dpam,'PAM')                
#     for gi in trange(len(dseq),desc='designing guides'):
    gierrfltmutpos=[]
    gierrdenan=[]
    gierrfltguidel=[]
    gierrpamnotfound=[]
    gierrcannotmutate=[]
    dseq_cols=['transcript: id','aminoacid: position','aminoacid: wild-type','codon: wild-type','id',]
    for gi in dseq.index:
        if cfg['mutations']=='mutations':
            dseqi=pd.DataFrame(dseq.loc[gi,dseq_cols+['amino acid mutation']]).T
            dmutagenesis_gi=pd.merge(dseqi,
                dmutagenesis,
                how='inner',
                left_on=['aminoacid: wild-type','codon: wild-type','amino acid mutation'],
                right_on=['amino acid','codon','amino acid mutation'])                    
        else:
            dseqi=pd.DataFrame(dseq.loc[gi,dseq_cols]).T
            dmutagenesis_gi=pd.merge(dseqi,
                dmutagenesis,
                how='inner',
                left_on=['aminoacid: wild-type','codon: wild-type'],
                right_on=['amino acid','codon'])        
        if len(dmutagenesis_gi)!=0:
            logging.info(f"working on {dseq.loc[gi,'id']}")
#             codon=dseq.loc[gi,'codon: wild-type']
            pos_codon=(flankaas)*3
            dpamsearches=get_pam_searches(dpam=dpam,
                 seq=dseq.loc[gi,'transcript: sequence'],
                 pos_codon=pos_codon,
                 test=test)
            # print(dpamsearches.columns) #RMME
            if dpamsearches is None:
                continue
            if len(dpamsearches)!=0:
                # filter by guide length
                dpamsearchesflt=dpamsearches.loc[dpamsearches['guide length']==dpamsearches['guide sequence length'],:]
                if len(dpamsearchesflt)!=0:
                    dpamsearches_strategy=pd.merge(dpamsearchesflt.reset_index(),dmutagenesis_gi.reset_index(),
                             how='inner',
                             on=['strand'],suffixes=['',': dmutagenesis_gi'])
                    if len(dpamsearches_strategy)!=0:                                 
                        if not 'dguides' in locals():
                            dguides=dpamsearches_strategy.copy()
                        else:
                            dguides=dguides.append(dpamsearches_strategy)
                        del dpamsearches_strategy
                    else:
                        gierrdenan.append(gi)
                        if dbug:
                            print('empty after removing nan seqs')
                else:
                    gierrfltguidel.append(gi)
                    if dbug:
                        print(f"empty after filtering by guide length. {dpamsearches['guide sequence length'].tolist()}")
            else:
                gierrpamnotfound.append(gi)
                if dbug:
                    print(f"no pam among {dpam.index.tolist()} were found {dseq.loc[gi,'transcript: sequence']}")
        else:
            gierrcannotmutate.append(gi)
            if dbug:
                print(f"can not mutate {dseqi['codon: wild-type'].tolist()}. its not in {dmutagenesis['codon'].tolist()}")

    gierrfltmutpos=[]
    gierrdenan=[]
    gierrfltguidel=[]
    gierrpamnotfound=[]
    gierrcannotmutate=[]

    err2idxs={'gierrfltmutpos':gierrfltmutpos,
              'gierrdenan':gierrdenan,
              'gierrfltguidel':gierrfltguidel,
              'gierrpamnotfound':gierrpamnotfound,
              'gierrcannotmutate':gierrcannotmutate,
             }

    if 'dguides' in locals():        
        # 0-based indexing 'position of guide ini/end', 'position of PAM ini/end'
        # 1-based indexing 'position of mutation in codon'
    
        logging.info('#reverse complement guides on negative strand sequences')
        dguides.loc[:,'PAM']=dguides.apply(lambda x : reverse_complement_multintseq(x['PAM'],nt2complement) if x['is a reverse complement'] else x['PAM'],axis=1)
        for colseq in ['guide+PAM sequence','guide sequence','PAM sequence']:
            dguides.loc[:,colseq]=dguides.apply(lambda x : str(str2seq(x[colseq]).reverse_complement()) if x['is a reverse complement'] else x[colseq],axis=1)

        logging.info('get dposition')
        dpositions=dguides.apply(lambda x: guide2dpositions(x),axis=1).apply(pd.Series)
#         posmut,posmutfrompam,distmutfrompam,posguideini,posguideend,activity_sequence
        dpositions.columns=['position of mutation','position of mutation from PAM','distance of mutation from PAM',
                           'position of guide ini','position of guide end','activity sequence']
#         dguides.to_csv('test_dguides.csv',sep='\t')
#         dpositions.to_csv('test_dposition.csv',sep='\t')
        for col in dpositions:
            dguides[col]=dpositions[col]
        
        logging.info('filter by # of editable nts in activity seq')
        logging.info(dguides.shape)
        dguides_noflt=dguides.copy()
        dguides=dguides.loc[(dguides.apply(lambda x : np.sum([x['activity sequence'].count(nt) for nt in  x['nucleotide']])==len(x['nucleotide']),axis=1)),:]
        logging.info(dguides.shape)
        if len(dguides)!=0:
            # if dbug:
            #     dguides.to_csv('test.tsv',sep='\t')
            logging.info('#filter by location of mutation within guide')
            dguides=dguides.loc[dguides.apply(lambda x : True if (x['distance of mutation from PAM: minimum']<=abs(x['distance of mutation from PAM'])<=x['distance of mutation from PAM: maximum']) else False,axis=1),:]
            if len(dguides)!=0:

                dguides.loc[:,'strategy']=dguides.apply(lambda x: f"{x['method']};{x['strand']};@{int(x['distance of mutation from PAM'])};{x['PAM']};{x['codon']}:{x['codon mutation']};{x['amino acid']}:{x['amino acid mutation']};",axis=1)
                dguides.loc[:,'guide: id']=dguides.apply(lambda x: f"{x['id']}|{int(x['aminoacid: position']) if not pd.isnull(x['aminoacid: position']) else 'nucleotide'}|({x['strategy']})",axis=1)
                dguides.loc[:,'guide+PAM length']=dguides.apply(lambda x: len(x['guide+PAM sequence']),axis=1)
                dguides=dguides.drop_duplicates(subset=['guide: id'])
                return dguides,dguides_noflt,err2idxs
            else:
                return None,dguides_noflt,None        
        else:
            return None,dguides_noflt,None        
    else:
        return None,None,None        

def dpam2dpam_strands(dpam,pams):
    """
    Duplicates dpam dataframe to be compatible for searching PAMs on - strand

    :param dpam: dataframe with pam information
    :param pams: pams to be used for actual designing of guides.
    """
    
    dpam=del_Unnamed(dpam)
    dpam['rPAM']=dpam.apply(lambda x : s2re(x['PAM'],multint2reg) ,axis=1)
    dpam=set_index(dpam,'PAM')
    dpam['strand']='+'
    dpamr=pd.DataFrame(columns=dpam.columns)
    dpam.loc[:,'reverse complement']=np.nan
    dpam.loc[:,'original']=np.nan
    for pam in dpam.index:
        pamr=reverse_complement_multintseq(pam,nt2complement)
        dpam.loc[pam,'reverse complement']=pamr
        dpam.loc[pam,'original']=pam
        dpamr.loc[pamr,'original']=pam
        dpam.loc[pam,'original position']=dpam.loc[pam,'position']
        dpamr.loc[pamr,'original position']=dpam.loc[pam,'position']
        dpamr.loc[pamr,['position','guide length','Description']]=dpam.loc[pam,['position','guide length','Description']]
        dpamr.loc[pamr,['rPAM']]=reverse_complement_multintseqreg(pam,multint2regcomplement,nt2complement)    
    dpamr['position']= dpamr.apply(lambda x: 'up' if x['position']=='down' else 'down',axis=1)
    dpamr['strand']='-'
    dpam_strands=dpam.append(dpamr,sort=True)
    dpam_strands.index.name='PAM'
    dpam_strands.loc[:,'is a reverse complement']=pd.isnull(dpam_strands.loc[:,'reverse complement'])
    pams_strands=pams+dpam_strands.loc[pams,'reverse complement'].dropna().tolist()
    dpam_strands=dpam_strands.loc[pams_strands,:]
    return dpam_strands

def dinnucleotide2dsequencesproper(dsequences,dmutagenesis,dbug=False):
    """
    Makes dseqeunces dataframe of nucleotide mutation format compatible to guide design modules

    :param dsequences: dsequences dataframe
    :param dmutagenesis: dmutagenesis dataframe
    """
    dmutagenesis=dmutagenesis.loc[(dmutagenesis['position of mutation in codon']==2),:]
    dsequences=pd.merge(dsequences,dmutagenesis,
             left_on=['nucleotide wild-type','nucleotide mutation','codon: wild-type','codon: mutation'],
             right_on=['nucleotide: wild-type','nucleotide: mutation','codon','codon mutation'],
             suffixes=['',': dmutagenesis'])
    if len(dsequences)!=0:
        dsequences['transcript: id']=dsequences['genome coordinate']
        dsequences['aminoacid mutation']=dsequences['amino acid mutation']
        dsequences['aminoacid: wild-type']=dsequences['amino acid']
        if dbug:
            dsequences['codon: wild-type']=dsequences['codon']
            df2info(dsequences,'nucle')
            df2info(dmutagenesis,'nucle')
        dsequences['id']=dsequences.apply(lambda x: f"{x['genome coordinate']}|{x['method']}|{x['mutation on strand'].replace(' strand','')}|{x['nucleotide wild-type']}:{x['nucleotide mutation']}|{x['codon: wild-type']}:{x['codon mutation']}",axis=1)
    else:
        logging.warning('empty dsequences after merging with dmutagenesis')
    cols_missing=[c for c in stepi2cols[1] if not c in dsequences]
    for c in cols_missing:
        dsequences[c]=np.nan
    return dsequences,dmutagenesis

def dseq2dguides(cfg):
    """
    Wrapper around make guides function.
    
    :param cfg: configuration dict.    
    """
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    dguideslinp=f"{cfg['datad']}/dguides.tsv"
    dguides_nofltp=f"{cfg['datad']}/dguides_noflt.tsv"
    dmutagenesisp=f"{cfg['datad']}/dmutagenesis.tsv"
    dpam_strandsp=f"{cfg['datad']}/dpam_strands.csv"
    if not exists(dguideslinp) or cfg['force']:
        dmutagenesis=pd.read_csv(f"{cfg[cfg['step']-1]}/dmutagenesis.tsv",sep='\t')
        if cfg['mutation_format']=='nucleotide':
            dsequences=pd.read_csv(f"{cfg[cfg['step']-2]}/dsequences.tsv",sep='\t') #FIXME if numbering of steps is changed, this is gonna blow
            dsequences,dmutagenesis=dinnucleotide2dsequencesproper(dsequences,dmutagenesis)
        elif cfg['mutation_format']=='aminoacid':
            dsequences=pd.read_csv(f"{cfg[cfg['step']-2]}/dsequences.tsv",sep='\t') #FIXME if numbering of steps is changed, this is gonna blow
            if 'reverse_mutations' in cfg:
                if cfg['reverse_mutations']: 
#                     from beditor.lib.global_vars import stepi2cols
                    cols_dsequences=dsequences.columns.tolist()
                    dsequences=pd.merge(dsequences,
                        dmutagenesis,
                        how='inner',
                        left_on=['aminoacid: wild-type','codon: mutation','amino acid mutation'],
                        right_on=['amino acid','codon mutation','amino acid mutation'],
                                       suffixes=['', ': dmutagenesis']) 
                    dsequences['codon: wild-type']=dsequences['codon']
                    dsequences=dsequences.loc[:,cols_dsequences]
                    
        dsequences.to_csv(f"{cfg[cfg['step']]}/dsequences.tsv",sep='\t') #FIXME if numbering of steps is changed, this is gonna blow
        # make pam table
        dpam=pd.read_table(f'{dirname(realpath(__file__))}/../data/dpam.tsv')
        if sum(dpam['PAM'].isin(cfg['pams']))!=len(cfg['pams']):
            logging.error(f"PAM/s {cfg['pams']} are not supported {dpam['PAM'].tolist()}")
            sys.exit(1)
        dpam_strands=dpam2dpam_strands(dpam,pams=cfg['pams'])
        dpam_strands.to_csv(dpam_strandsp,sep='\t')

        dmutagenesis['strand']=dmutagenesis.apply(lambda x : x['mutation on strand'].replace(' strand',''),axis=1)        
        dmutagenesis.to_csv(dmutagenesisp,sep='\t')
        
        dguideslin,dguides_noflt,err2idxs=make_guides(cfg,dsequences,
                    dmutagenesis,
                    dpam=dpam_strands,
                       test=cfg['test'],
                       # dbug=True,
                     )
        if not dguides_noflt is None:
            dguides_noflt.to_csv(dguides_nofltp,sep='\t')
        if not ((dguideslin is None) and (err2idxs is None)):
            dguideslin.to_csv(dguideslinp,sep='\t')
            if cfg['test']:
                logging.info(err2idxs)
            with open(dguideslinp+'.err.json', 'w') as f:
                json.dump(err2idxs, f)
        else:
            from beditor.lib.global_vars import saveemptytable
            logging.warning('no guides designed; saving an empty table.')
            saveemptytable(cfg,dguideslinp)
        import gc
        gc.collect()