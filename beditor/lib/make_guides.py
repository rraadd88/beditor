import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

import logging
from tqdm import trange

def make_guides(dseq,dmutagenesis,
               test=False,
               dbug=False):
    """
    Makes guides by
    1. searching all PAM sequences on 'both' the strands,
    2. filtering guides by all possible strategies (given in dmutagenesis) e.g. activity window,
    Finally generates a table.
    """
    dguides=dseq.copy()
    dguides=dguides.reset_index()
    dguides.index=range(len(dguides))
    for strand in dmutagenesis.loc[:,'mutation on strand'].unique():
        logging.info('working on {}'.format(strand))
        for genei in trange(len(dguides)-1,desc='mutations'):
            gene=dguides.loc[genei,'gene: name']
            seq=dguides.loc[genei,'transcript: sequence']
    #             pos_codon=(int(dguides.loc[genei,'aminoacid: position'])-1)*3 # 0based index
            pos_codon=(7)*3 # 0based index
            codon=dguides.loc[genei,'codon: wild-type']
            if test:
                if codon!=seq[(pos_codon):pos_codon+3]:
                    print('pos_codon is wrong')
                    break
            if strand=='- strand':
                print(seq)
                seq=str(Seq.Seq(seq,Alphabet.generic_dna).reverse_complement())
                codon=str(Seq.Seq(codon,Alphabet.generic_dna).reverse_complement())
                pos_codon=len(seq)-(pos_codon)-3
            posGGs=[i for i in range(len(seq)) if seq.startswith('GG',i)]
            if len(posGGs)!=0:
                for posGGi,posGG in enumerate(posGGs):
                    seq_target=seq[posGG-21:posGG-1]
                    pos_codon_from_PAM=pos_codon-(posGG)+1
                    for muti in dmutagenesis[(dmutagenesis['codon']==codon) & (dmutagenesis['mutation on strand']==strand)].index:
                        method=dmutagenesis.loc[muti,'method']
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
                                    dguides.loc[genei,'guide sequence+PAM({0})'.format(strategy)]=seq_target+seq[posGG-1:posGG+2]
                                    if test:
                                        print('{}:pos_mut_from_PAM={};pos_codon_from_PAM={};seq_activity={};{}'.format(gene,pos_mut_from_PAM,pos_codon_from_PAM,seq_activity,strategy))
                                    if dbug:
                                        print(posGG)
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
    #                     if posGG==45:
    #                         break
    #                 break
    #             break
    #         break
    #     break
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
