import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists,abspath,dirname

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

#lib modules
import logging
from beditor.lib.global_vars import hosts

def translate(dnaseq,host='human',fmtout=str):
    if isinstance(dnaseq,str): 
        dnaseq=Seq.Seq(dnaseq,Alphabet.generic_dna)
    prtseq=dnaseq.translate(table=hosts[host])
    if fmtout is str:
        return str(prtseq)
    else:
        return prtseq

def get_seq_yeast(dseq,orfs_fastap,
            host,
           test=False):
        
    orfs=SeqIO.to_dict(SeqIO.parse(orfs_fastap,'fasta'))
    for subi,sub in enumerate(dseq['gene: name'].tolist()):
        if sub in orfs:
            seq=orfs[sub]
            if test:
                print(seq.id)
    #             break
        else:
            print('{}: not found in fasta'.format(sub))
            break
        dnaseq=str(seq.seq)
#         print(host)
        prtseq=translate(seq.seq,host=host,fmtout=str)
#             prtseq=str(seq.seq.translate(table=hosts[host]))
        
        dseq.loc[subi,'DNA sequence']=dnaseq
        dseq.loc[subi,'Protein sequence']=prtseq
#         if and (not '|' in dseq.loc[subi,'aminoacid: position']):
        dseq.loc[subi,'gene: name']=sub
#             dseq.loc[subi,'aminoacid: position']=int(dseq.loc[subi,'aminoacid: position'])
        dseq.loc[subi,'aminoacid: wild-type']=prtseq[int(dseq.loc[subi,'aminoacid: position'])-1]
        if dseq.loc[subi,'aminoacid: wild-type']=='Y':
            dseq.loc[subi,'codon: wild-type']=dnaseq[(int(dseq.loc[subi,'aminoacid: position'])-1)*3:(int(dseq.loc[subi,'aminoacid: position'])-1)*3+3]
            dseq.loc[subi,'transcript: sequence']=dnaseq
            dseq.loc[subi,'protein: sequence']=prtseq

#             psites=[int(i) for i in dseq.iloc[subi,:]['aminoacid: position'].split('|')]
#             for psitei,psite in enumerate(psites):
#                 if prtseq[int(psite)-1] == 'Y':
#                     dseq.loc[subi,'Psite{0:02d}'.format(psitei+1)]=psite                
#                 else:
#                     if test:
#                         print('{}: p site not found; found {}; seq-10+10: {}'.format(sub,prtseq[int(psite)-1],prtseq[int(dint.iloc[subi,:]['aminoacid: position'])-10:int(dint.iloc[subi,:]['aminoacid: position'])+10]))
#                     if prtseq[int(psite)] == 'Y':
#                         dseq.loc[subi,'Psite{0:02d}'.format(psitei+1)]=psite+1
#                         if test:
#                             print('{}: p site corrected'.format(sub))
#                     else:
#                         if test:
#                             print('{}: p site not found; found {}; seq-10+10: {}'.format(sub,prtseq[int(psite)-1],prtseq[int(dint.iloc[subi,:]['aminoacid: position'])-10:int(dint.iloc[subi,:]['aminoacid: position'])+10]))
#                         break
#     #         print(dseq.shape)
    dseq=dseq.dropna(axis=0,how='any')
#         print(dseq.shape)
    return dseq

def get_codon_seq(dintseqflt01,test=False):
    for subi,sub in zip(dintseqflt01.index,dintseqflt01['gene: name'].tolist()):
        psite=int(dintseqflt01.loc[subi,'Psite01'])
        dintseqflt01.loc[subi,'P']=dintseqflt01.loc[subi,'Protein sequence'][(psite-1)]
        dintseqflt01.loc[subi,'P-codon']=dintseqflt01.loc[subi,'DNA sequence'][(psite-1)*3:(psite-1)*3+3]    
        ini=(psite-1)-10
        end=(psite-1)+10+1
        if ini<0:
            ini=0
        if end>len(dintseqflt01.loc[subi,'Protein sequence']):
            end=len(dintseqflt01.loc[subi,'Protein sequence'])        
        dintseqflt01.loc[subi,'10[P]10']=dintseqflt01.loc[subi,'Protein sequence'][ini:end]
        dintseqflt01.loc[subi,'10[P]10: P position']=(psite-1)-ini
        ini=(psite-1)*3-(10*3)
        end=(psite-1)*3+((10+1)*3)
        if ini<0:
            ini=0
        if end>len(dintseqflt01.loc[subi,'DNA sequence']):
            end=len(dintseqflt01.loc[subi,'DNA sequence'])        
        dintseqflt01.loc[subi,'30[P-codon]30']=dintseqflt01.loc[subi,'DNA sequence'][ini:end]
        dintseqflt01.loc[subi,'30[P-codon]30: P-codon position']=(psite-1)*3-ini
        if test:
            print('{}:{}'.format(dintseqflt01.loc[subi,'P'],dintseqflt01.loc[subi,'P-codon']))
    return dintseqflt01

def din2dseq(cfg):
    # get dna and protein sequences 
    dseqp='{}/dseq.csv'.format(cfg['datad'])
    if not exists(dseqp) or cfg['force']:
        din=pd.read_csv(cfg['dinp'])    
        if cfg['host']=='human':
            import pyensembl
            #import ensembl object that would fetch genes 
            ensembl = pyensembl.EnsemblRelease(release=cfg['release'])
            dseq=din.copy()
            dseq.index=range(len(dseq.index))
            for rowi in dseq.index:
                gene_id=dseq.loc[rowi,'gene: id']
                transcript=ensembl.transcript_by_id(dseq.loc[rowi,'transcript: id'])
                dnaseq=transcript.coding_sequence
                prtseq=transcript.protein_sequence
                dseq.loc[rowi,'aminoacid: wild-type']=prtseq[dseq.loc[rowi,'aminoacid: position']-1]
                dseq.loc[rowi,'codon: wild-type']=dnaseq[(dseq.loc[rowi,'aminoacid: position']-1)*3:(dseq.loc[rowi,'aminoacid: position']-1)*3+3]
                dseq.loc[rowi,'transcript: sequence']=dnaseq
                dseq.loc[rowi,'protein: sequence']=prtseq
        else:
            dseq=get_seq_yeast(din,
                          orfs_fastap='{}/../data/yeast/orf_coding_all.fasta'.format(abspath(dirname(__file__))),
                          host=cfg['host'],
                          test=cfg['test'])
        din.to_csv('{}/din.csv'.format(cfg['datad']))
        dseq.to_csv(dseqp)
        logging.info('Counts of amino acids to mutate:')
        logging.info(dseq['aminoacid: wild-type'].value_counts())