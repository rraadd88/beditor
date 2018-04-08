import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

def make_guides(dseq,dmutagenesis,
               test=False,
               dbug=False):
    dguides=dseq.copy()
    for strand in dmutagenesis.loc[:,'mutation on strand'].unique():
        for subi,sub in zip(dguides.index,dguides['gene: name'].tolist()):
            seq=dguides.loc[subi,'transcript: sequence']
            pos_codon=(int(dguides.loc[subi,'aminoacid: position'])-1)*3 # 0based index
            codon=dguides.loc[subi,'codon: wild-type']
            if test:
                if codon!=seq[(pos_codon):pos_codon+3]:
                    print('pos_codon is wrong')
                    break
            if strand=='- strand':
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
#                                     dguides.loc[subi,'guide sequence ({0})'.format(strategy)]=seq_target
                                    dguides.loc[subi,'guide sequence+PAM({0})'.format(strategy)]=seq_target+seq[posGG-1:posGG+2]
                                    if test:
                                        print('{}:pos_mut_from_PAM={};pos_codon_from_PAM={};seq_activity={};{}'.format(sub,pos_mut_from_PAM,pos_codon_from_PAM,seq_activity,strategy))
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

from global_vars import BEs,pos_muts
def dseq2dguides(cfg):
    dseq=pd.read_csv('{}/dseq.csv'.format(cfg['datad']))
    dmutagenesis=pd.read_csv('{}/dmutagenesis.csv'.format(cfg['datad']))

    dguides=make_guides(dseq,dmutagenesis)

    genes_dguides=dguides.loc[:,[c for c in dguides if c.startswith('guide sequence+PAM')]].dropna(how='all')
    dguidesflt=dguides.loc[genes_dguides.index,:]
    print('Out of total {} sites, for {} sites,\nguides were designed.'.format(len(dseq),len(dguidesflt)))

    genes_dguides=dguidesflt.loc[:,[c for c in dguides if (c.startswith('guide sequence+PAM')) \
                                 and ((' to A' in c) or (' to L' in c) or (' to H' in c)) ]].dropna(how='all')
    dguidesflt1=dguidesflt.loc[genes_dguides.index,:]
    print('Out of total {} sites, for {} sites,\nguides were designed for mutations to A, L or H.'.format(len(dguidesflt),len(dguidesflt1)))
    

    dguides.to_csv('{}/dguides.csv'.format(cfg['datad']))
    dguidesflt.to_csv('{}/dguidesflt.csv'.format(cfg['datad']))
    dguidesflt1.to_csv('{}/dguidesflt1.csv'.format(cfg['datad']))
    
    from beditor.lib.plot_res import plot_nt_composition
    plot_nt_composition(dguidesflt1,
                    pos_muts,
                    plotp='{}/nt_compositions_of_guides.png'.format(plotd))
