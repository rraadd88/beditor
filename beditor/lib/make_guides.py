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


def get_pam_searches(dpam,seq,pos_codon,
                    test=False):
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
    dpamposs['codon']=seq[pos_codon:pos_codon+3]
    dpamposs['guide sequence']=dpamposs['guide sequence'].fillna('')
    if len(dpamposs)==0:
        return None
    dpamposs['guide sequence length']=dpamposs.apply(lambda x : len(x['guide sequence']),axis=1)

    
    dpamposs=dpam.join(set_index(dpamposs,'PAM'),how='right')
    dpamposs.index.name='PAM'
    return dpamposs


def get_activity_seq(guide_seq,pam_pos,pam_dist_min,pam_dist_max,dbug=False):
    if dbug:
        print(guide_seq,pam_pos,pam_dist_min,pam_dist_max)
    if pam_pos=='up':
        seq=guide_seq[pam_dist_min-1:pam_dist_max]
    elif pam_pos=='down':
        seq=guide_seq[::-1][pam_dist_min-1:pam_dist_max][::-1]
    return seq

def get_pos_mut_from_pam(x,
    pos_codon_ini=22, #FIXME if fragment size if changed
    pos_codon_end=24, #FIXME if fragment size if changed
                        ):
    # For debugging
        # x=dguides.loc[0,:]
    strand=x['strand'] #FIXME annotate the strand wrt to fragment(BE/dPAM) and genome
    loc_pam=x['position']
#     #keep all 1-based indexed
#     pos_guide_ini=x['position of guide ini']+1
#     pos_guide_end=x['position of guide end']+1
    pos_pam_ini=x['position of PAM ini']+1
    pos_pam_end=x['position of PAM end']+1
    pos_mut_incodon=x['position of mutation in codon']

    pos_mut_infragment=pos_mut_incodon+pos_codon_ini-1

    if loc_pam=='down' and strand=='+':
        dist_mut_frompam=pos_mut_infragment-pos_pam_ini
    elif loc_pam=='down' and strand=='-':
        dist_mut_frompam=pos_mut_infragment-pos_pam_end
        dist_mut_frompam=dist_mut_frompam*-1
    elif loc_pam=='up' and strand=='+':
        dist_mut_frompam=pos_mut_infragment-pos_pam_end        
    elif loc_pam=='up' and strand=='-':
        dist_mut_frompam=pos_mut_infragment-pos_pam_end        
        dist_mut_frompam=dist_mut_frompam*-1
    else:
        raise(ValueError(f'loc_pam or strand are invalid ({loc_pam},{strand})'))
    return dist_mut_frompam

def make_guides(cfg,dseq,dmutagenesis,
                dpam,
               test=False,
               dbug=False):
    """
    Makes guides by
    1. searching all PAM sequences on 'both' the strands,
    2. filtering guides by all possible strategies (given in dmutagenesis) e.g. activity window,
    Finally generates a table.
    0-based indexing
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
                             on=['codon','strand'])
                    if len(dpamsearches_strategy)!=0:                                 
                        # print(dpamsearches_strategy.columns) #RMME
                        if not 'dguides' in locals():
                            dguides=dpamsearches_strategy.copy()
                        else:
                            dguides=dguides.append(dpamsearches_strategy)
#                             if dbug:
#                                 return dpamsearches_strategy        
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
        dguides['distance of mutation in codon from PAM']=dguides.apply(lambda x : get_pos_mut_from_pam(x),axis=1)
        dguides.loc[:,'position of guide ini']=dguides.apply(lambda x : x['position of PAM end'] if x['position']=='up' else x['position of PAM end']-x['guide sequence length']-1,axis=1)
        dguides.loc[:,'position of guide end']=dguides.apply(lambda x : x['position of PAM ini']-1 if x['position']=='down' else x['position of PAM ini']+x['guide sequence length'],axis=1)
        
        # print(dguides.columns) #RMME
        logging.info('#reverse complement guides on negative strand sequences')
        dguides.loc[:,'PAM']=dguides.apply(lambda x : reverse_complement_multintseq(x['PAM'],nt2complement) if x['is a reverse complement'] else x['PAM'],axis=1)
        for colseq in ['guide+PAM sequence','guide sequence','PAM sequence']:
            dguides.loc[:,colseq]=dguides.apply(lambda x : str(str2seq(x[colseq]).reverse_complement()) if x['is a reverse complement'] else x[colseq],axis=1)

        logging.info('get activity sequence')
        dguides['activity sequence']=dguides.apply(lambda x: get_activity_seq(x['guide sequence'],x['original position'],
                     int(x['distance of mutation from PAM: minimum']),
                     int(x['distance of mutation from PAM: maximum'])),axis=1)   

        logging.info('filter by # of editable nts in activity seq')
        logging.info(dguides.shape)
        dguides=dguides.loc[(dguides.apply(lambda x : np.sum([x['activity sequence'].count(nt) for nt in  x['nucleotide']])==len(x['nucleotide']),axis=1)),:]
        logging.info(dguides.shape)
        dguides_noflt=dguides.copy()
        if len(dguides)!=0:
            # if dbug:
            #     dguides.to_csv('test.tsv',sep='\t')
            logging.info('#filter by location of mutation within guide')
            dguides=dguides.loc[dguides.apply(lambda x : True if (x['distance of mutation from PAM: minimum']<=abs(x['distance of mutation in codon from PAM'])<=x['distance of mutation from PAM: maximum']) else False,axis=1),:]
            if len(dguides)!=0:

                dguides.loc[:,'strategy']=dguides.apply(lambda x: f"{x['method']};{x['strand']};@{int(x['distance of mutation in codon from PAM'])};{x['PAM']};{x['codon']}:{x['codon mutation']};{x['amino acid']}:{x['amino acid mutation']};",axis=1)
                dguides.loc[:,'guide: id']=dguides.apply(lambda x: f"{x['id']}|{int(x['aminoacid: position'])}|({x['strategy']})",axis=1)
                dguides.loc[:,'guide+PAM length']=dguides.apply(lambda x: len(x['guide+PAM sequence']),axis=1)
                dguides=dguides.drop_duplicates(subset=['guide: id'])
                return dguides,dguides_noflt,err2idxs
            else:
                return None,None,None        
        else:
            return None,None,None        
    else:
        return None,None,None        

def dpam2dpam_strands(dpam,pams):
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
    
def dseq2dguides(cfg):
    """
    Wrapper around make guides function.
    :param cfg: conffguration settings given in yml file.    
    """
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    dguideslinp='{}/dguides.tsv'.format(cfg['datad'])
    dguides_nofltp='{}/dguides_noflt.tsv'.format(cfg['datad'])
    dmutagenesisp='{}/dmutagenesis.tsv'.format(cfg['datad'])
    dpam_strandsp='{}/dpam_strands.csv'.format(cfg['datad'])
    if not exists(dguideslinp) or cfg['force']:
        dseq=pd.read_csv('{}/dsequences.tsv'.format(cfg[cfg['step']-2]),sep='\t') #FIXME if numbering of steps is changed, this is gonna blow
        dmutagenesis=pd.read_csv(f"{cfg[cfg['step']-1]}/dmutagenesis.tsv",sep='\t')
        # make pam table
        dpam=pd.read_table('{}/../data/dpam.tsv'.format(dirname(realpath(__file__))))
        if sum(dpam['PAM'].isin(cfg['pams']))!=len(cfg['pams']):
            logging.error(f"PAM/s {cfg['pams']} are not supported {dpam['PAM'].tolist()}")
            sys.exit(1)
        dpam_strands=dpam2dpam_strands(dpam,pams=cfg['pams'])
        dpam_strands.to_csv(dpam_strandsp,sep='\t')

        for col in ['Position of mutation from PAM: ','Position of mutation from PAM: ',
                    'Position of codon start from PAM: ','Position of codon start from PAM: ']:    
            mum1,mum2=('minimum','maximum')
        #             for mum1,mum2 in zip(['minimum','maximum'],['minimum','maximum']):
            dmutagenesis[col.replace('Position','distance')+mum1]=dmutagenesis.apply(lambda x : np.min([abs(x[col+'minimum']),abs(x[col+'maximum'])]),axis=1)
            dmutagenesis[col.replace('Position','distance')+mum2]=dmutagenesis.apply(lambda x : np.max([abs(x[col+'minimum']),abs(x[col+'maximum'])]),axis=1)

        dmutagenesis['strand']=dmutagenesis.apply(lambda x : x['mutation on strand'].replace(' strand',''),axis=1)        

        dmutagenesis.to_csv(dmutagenesisp,sep='\t')
#         sys.exist(1)
        
        dguideslin,dguides_noflt,err2idxs=make_guides(cfg,dseq,
                    dmutagenesis,
                    dpam=dpam_strands,
                       test=cfg['test'],
                       # dbug=True,
                     )
        if not ((dguideslin is None) and (err2idxs is None)):
            dguideslin.to_csv(dguideslinp,sep='\t')
            dguides_noflt.to_csv(dguides_nofltp,sep='\t')
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