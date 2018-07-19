import pandas as pd
import numpy as np

from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

import logging
from tqdm import trange

def get_guide_seq(dguides,method):
    seq_activity=seq_target[20+dmutagenesis.loc[muti,'Position of mutation from PAM: minimum']:20+1+dmutagenesis.loc[muti,'Position of mutation from PAM: maximum']]
    if seq_activity.count(dmutagenesis.loc[muti,'nucleotide'])==1:
        if (pos_codon_from_PAM>=dmutagenesis.loc[muti,'Position of codon start from PAM: minimum']) and (pos_codon_from_PAM<=dmutagenesis.loc[muti,'Position of codon start from PAM: maximum']):
            if strand=='+ strand':
                pos_mut_from_PAM=pos_codon_from_PAM-1+dmutagenesis.loc[muti,'position of mutation in codon']                    
            elif strand=='- strand':
                pos_mut_from_PAM=pos_codon_from_PAM-1+4-dmutagenesis.loc[muti,'position of mutation in codon']                    
            if (pos_mut_from_PAM>=dmutagenesis.loc[muti,'Position of mutation from PAM: minimum']) and (pos_mut_from_PAM<=dmutagenesis.loc[muti,'Position of mutation from PAM: maximum']):
                strategy='{}; {}: {} to {}; {} to {}, codon position={}; mutation position={};'.format(dmutagenesis.loc[muti,'mutation on strand'],
                                                             method,
                                                            dmutagenesis.loc[muti,'codon'],
                                                            dmutagenesis.loc[muti,'codon mutation'],
                                                            dmutagenesis.loc[muti,'amino acid'],
                                                            dmutagenesis.loc[muti,'amino acid mutation'],
                                                            pos_codon_from_PAM,
                                                            pos_mut_from_PAM
                                                            )        
                codon_mut=dmutagenesis.loc[muti,'codon mutation']
                if strand=='- strand':
                    codon_mut=str(Seq.Seq(codon_mut,Alphabet.generic_dna).reverse_complement())
#                                     dguides.loc[genei,'guide sequence ({0})'.format(strategy)]=seq_target
                dguides.loc[genei,'guide sequence+PAM({0})'.format(strategy)]=seq_target+seq[pampos-1:pampos+2]
                if test:
                    print('{}:pos_mut_from_PAM={};pos_codon_from_PAM={};seq_activity={};{}'.format(gene,pos_mut_from_PAM,pos_codon_from_PAM,seq_activity,strategy))
                if dbug:
                    print(pampos)
                    print(strand)
                    print(codon)
                    print(seq)
                    print(seq_target)
                    print(seq_activity)
                    print(pos_codon_from_PAM)
                    print(pos_mut_from_PAM)
#                                     print(seq_guide)
                    print(sdf)
                    break
    return dguides

def get_seq_codon(seq,codon,pos_codon):            
    if strand=='- strand':
        # print(seq)
        seq=str(Seq.Seq(seq,Alphabet.generic_dna).reverse_complement())
        codon=str(Seq.Seq(codon,Alphabet.generic_dna).reverse_complement())
        pos_codon=len(seq)-(pos_codon)-3
    if test:
        if codon!=seq[(pos_codon):pos_codon+3]:
            print('pos_codon is wrong')
            break
    return seq,codon,pos_codon

def get_pam_searches(pam,dpam,seq,pos_codon,
                    test=False):
    if dpam.loc[pam,'position']=='down':
        pamposs=[match.start() for match in re.finditer(pam, seq)]
    elif dpam.loc[pam,'position']=='up':
        pamposs=[match.end() for match in re.finditer(pam, seq)]                    
    for pamposi,pampos in enumerate(pamposs):
        if dpam.loc[pam,'position']=='down':
            seq_target=seq[pampos-dpam.loc[pam,'guide lendth']+1:pampos-1]
            pos_codon_from_PAM=pos_codon-(pampos)+1
        elif dpam.loc[pam,'position']=='up':
            seq_target=seq[pampos+1:pampos+dpam.loc[pam,'guide lendth']+1]
            pos_codon_from_PAM=pos_codon-(pampos)+1                        
    return dpamposs

def make_guides(dseq,dmutagenesis,
                pams,dpam,
               test=False,
               dbug=False):
    """
    Makes guides by
    1. searching all PAM sequences on 'both' the strands,
    2. filtering guides by all possible strategies (given in dmutagenesis) e.g. activity window,
    Finally generates a table.
    0-based indexing
    """
    flankaas=7#FIXME if flank length changes
    dguides=dseq.copy()
    dguides=dguides.reset_index()
    dguides.index=range(len(dguides))
    dpam=set_index(dpam,'PAM')                
    for pam in pams:
        logging.info(f'working on {pam} PAM')
        for strand in dmutagenesis.loc[:,'mutation on strand'].unique():
            logging.info('working on {strand}')
            for genei in range(len(dguides)):
                gene=dguides.loc[genei,'gene: name']
                logging.info('working on {gene}')
                seq,codon,pos_codon=get_seq_codon(seq=dguides.loc[genei,'transcript: sequence'],
                                    codon=dguides.loc[genei,'codon: wild-type'],
                                    pos_codon=(flankaas)*3,test=test)
                dpamposs=get_pam_searches(pam,dpam,seq,pos_codon,test=test)
                for pamposi in dpamposs.index:
                    for muti in dmutagenesis[(dmutagenesis['codon']==codon) & (dmutagenesis['mutation on strand']==strand)].index:
                        method=dmutagenesis.loc[muti,'method']
                        
    return dguides

def dseq2dguides(cfg):
    """
    Wrapper around make guides function.
    :param cfg: conffguration settings given in yml file.    
    """
    from beditor.lib.global_vars import BEs,pos_muts
    from beditor.lib.io_dfs import df2unstack
    from beditor.lib.plot_res import plot_nt_composition

    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    dguidesp='{}/dguides.csv'.format(cfg['datad'])
    dguideslinp='{}/dguideslin.csv'.format(cfg['datad'])
    if not exists(dguidesp) or cfg['force']:
        dseq=pd.read_csv('{}/dseq.csv'.format(cfg[cfg['step']-2])) #FIXME if numbering of steps is changed, this is gonna blow
        dmutagenesis=pd.read_csv('{}/dmutagenesis.csv'.format(cfg[cfg['step']-1]))

        dguides=make_guides(dseq,dmutagenesis)
        
        colns_guideseq=[c for c in dguides if c.startswith('guide sequence+PAM')]
        genes_dguides=dguides.loc[:,colns_guideseq].dropna(how='all')
        dguides=dguides.loc[genes_dguides.index,:]
        print('Out of total {} sites, for {} sites,\nguides were designed.'.format(len(dseq),len(dguides)))
        dguides.to_csv(dguidesp)


        dguideslin=dguides.set_index('id').loc[:,colns_guideseq]
        dguideslin.index.name='guide: id'

        dguideslin=df2unstack(dguideslin,coln='strategy',col='guide sequence+PAM').dropna(axis=0,how='any')
        dguideslin['guide: id']=dguideslin.apply(lambda x : x['guide: id']+'|'+x['strategy'],axis=1)
        dguideslin.to_csv(dguideslinp)

        plot_nt_composition(dguides,
                        pos_muts,
                        plotp='{}/nt_compositions_of_guides.png'.format(cfg['plotd']))
