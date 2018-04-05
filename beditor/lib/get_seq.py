import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO


def get_seq(dintseq,orfs_fastap,
            host,
           test=False):
        
    orfs=SeqIO.to_dict(SeqIO.parse(orfs_fastap,'fasta'))
    for subi,sub in enumerate(dintseq['Substrate'].tolist()):
        if sub in orfs:
            seq=orfs[sub]
            if test:
                print(seq.id)
    #             break
        else:
            print('{}: not found in fasta'.format(sub))
            break
        dnaseq=str(seq.seq)
        prtseq=str(seq.seq.translate(table=hosts[host]))
        dintseq.loc[subi,'DNA sequence']=dnaseq
        dintseq.loc[subi,'Protein sequence']=prtseq        
        psites=[int(i) for i in dintseq.iloc[subi,:]['Ysite'].split('|')]
        for psitei,psite in enumerate(psites):
            if prtseq[int(psite)-1] == 'Y':
                dintseq.loc[subi,'Psite{0:02d}'.format(psitei+1)]=psite
            else:
                if test:
                    print('{}: p site not found; found {}; seq-10+10: {}'.format(sub,prtseq[int(psite)-1],prtseq[int(dint.iloc[subi,:]['Ysite'])-10:int(dint.iloc[subi,:]['Ysite'])+10]))
                if prtseq[int(psite)] == 'Y':
                    dintseq.loc[subi,'Psite{0:02d}'.format(psitei+1)]=psite+1
                    if test:
                        print('{}: p site corrected'.format(sub))
                else:
                    if test:
                        print('{}: p site not found; found {}; seq-10+10: {}'.format(sub,prtseq[int(psite)-1],prtseq[int(dint.iloc[subi,:]['Ysite'])-10:int(dint.iloc[subi,:]['Ysite'])+10]))
                    break
    return dintseq

def get_codon_seq(dintseqflt01,test=False):
    for subi,sub in zip(dintseqflt01.index,dintseqflt01['Substrate'].tolist()):
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

