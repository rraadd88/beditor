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
#       prtseq=str(seq.seq.translate(table=hosts[host]))
#       if not isinstance(dseq.loc[subi,'aminoacid: position'],str):
        dseq.loc[subi,'transcript: sequence']=dnaseq
        dseq.loc[subi,'Protein sequence']=prtseq
        dseq.loc[subi,'gene: name']=sub
    #   dseq.loc[subi,'aminoacid: position']=int(dseq.loc[subi,'aminoacid: position 
        try:
            dseq.loc[subi,'aminoacid: position']=int(float(dseq.loc[subi,'aminoacid: position']))
            dseq.loc[subi,'aminoacid: wild-type']=prtseq[int(dseq.loc[subi,'aminoacid: position'])-1]
        #   if dseq.loc[subi,'aminoacid: wild-type']=='Y':
            dseq.loc[subi,'codon: wild-type']=dnaseq[(int(dseq.loc[subi,'aminoacid: position'])-1)*3:(int(dseq.loc[subi,'aminoacid: position'])-1)*3+3]
            dseq.loc[subi,'transcript: sequence']=dnaseq
            dseq.loc[subi,'protein: sequence']=prtseq
        except:
            dseq.loc[subi,'aminoacid: position']=dseq.loc[subi,'aminoacid: position']
            if test:
                print(print(seq.id,dseq.loc[subi,'aminoacid: position']))                        
    dseq=dseq.dropna(axis=0,how='any')
#   print(dseq.shape)
    return dseq

def get_codon_seq(dintseqflt01,test=False):
    for subi,sub in zip(dintseqflt01.index,dintseqflt01['gene: name'].tolist()):
        psite=int(dintseqflt01.loc[subi,'Psite01'])
        dintseqflt01.loc[subi,'P']=dintseqflt01.loc[subi,'Protein sequence'][(psite-1)]
        dintseqflt01.loc[subi,'P-codon']=dintseqflt01.loc[subi,'transcript: sequence'][(psite-1)*3:(psite-1)*3+3]    
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
        if end>len(dintseqflt01.loc[subi,'transcript: sequence']):
            end=len(dintseqflt01.loc[subi,'transcript: sequence'])        
        dintseqflt01.loc[subi,'30[P-codon]30']=dintseqflt01.loc[subi,'transcript: sequence'][ini:end]
        dintseqflt01.loc[subi,'30[P-codon]30: P-codon position']=(psite-1)*3-ini
        if test:
            print('{}:{}'.format(dintseqflt01.loc[subi,'P'],dintseqflt01.loc[subi,'P-codon']))
    return dintseqflt01

def din2dseq(cfg):
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    # get dna and protein sequences 
    dseqp='{}/dseq.csv'.format(cfg['datad'])
    if not exists(dseqp) or cfg['force']:
        din=pd.read_csv(cfg['dinp'])    
        if cfg['host']=='human':
            import pyensembl
            #import ensembl object that would fetch genes 
            ensembl = pyensembl.EnsemblRelease(release=cfg['release'])
            
            din.index=range(len(din))
            dbedp='{}/dbedflank.bed'.format(cfg['datad'])
            dbed=pd.DataFrame(columns=bed_colns)
            terrpositions=[]
            bedrowi=0
            for i in tnrange(len(din)-1,desc='get positions for bedtools'):
                t=ensembl.transcript_by_id(din.loc[i,'transcript: id'])
                coding_sequence_positions=tboundaries2positions(t)
                if len(coding_sequence_positions)==len(t.coding_sequence):
                #TODO     need to check if the seq made from coding_sequence_positions is same as t.coding_seqeunce
                    dbed.loc[bedrowi,'chromosome']=t.contig
                    dcoding=t2pmapper(t,coding_sequence_positions)
                    dcodingmutpos=dcoding.loc[(dcoding['protein index']==din.loc[i,'aminoacid: position']),:]
                    if t.strand=='+':
                        dbed.loc[bedrowi,'codon start'],dbed.loc[bedrowi,'codon end']=dcodingmutpos['coding sequence positions'].tolist()[0],dcodingmutpos['coding sequence positions'].tolist()[1]
                    elif t.strand=='-':
                        dbed.loc[bedrowi,'codon end'],dbed.loc[bedrowi,'codon start']=dcodingmutpos['coding sequence positions'].tolist()[0],dcodingmutpos['coding sequence positions'].tolist()[1]
                    dbed.loc[bedrowi,'start'],dbed.loc[bedrowi,'end']=dbed.loc[bedrowi,'codon start']-20,dbed.loc[bedrowi,'codon end']+20
                    dbed.loc[bedrowi,'reference residue']=dcodingmutpos['protein sequence'].tolist()[0]
                    dbed.loc[bedrowi,'reference codon']=''.join(dcodingmutpos['coding sequence'].tolist())
                    dbed.loc[bedrowi,'strand']=t.strand
                    dbed.loc[bedrowi,'id']='{}|{}|{}|{}|{}'.format(din.loc[i,'transcript: id'],
                                                                  dbed.loc[bedrowi,'chromosome'],
                             dbed.loc[bedrowi,'strand'],dbed.loc[bedrowi,'start'],dbed.loc[bedrowi,'end'])        
                    dbed.loc[bedrowi,'gene: id']=t.gene_id
                    dbed.loc[bedrowi,'gene: name']=t.gene.name
                    dbed.loc[bedrowi,'protein: id']=t.protein_id
                    dbed.loc[bedrowi,'aminoacid: position']=din.loc[i,'aminoacid: position']
            #         break
                    bedrowi+=1
                else:
                    terrpositions.append(t.id)
            dbed.loc[:,'start']=dbed.loc[:,'start'].astype(int)
            dbed.loc[:,'end']=dbed.loc[:,'end'].astype(int)

            print(terrpositions)
            dbed.loc[:,bed_colns].to_csv(dbedp,sep='\t',
                            header=False,index=False)

            bedp='{}/dbedflank.bed'.format(cfg['datad'])
            fastap='{}/dbedflank.fa'.format(cfg['datad'])
            cmd='bedtools getfasta -s -name -fi {} -bed {} -fo {}'.format(cfg['genomep'],bedp,fastap)
            runbashcmd(cmd)

            dflankfa=fa2df(fastap,ids2cols=True)
            dflankfa.loc[:,'sequence']=dflankfa.loc[:,'sequence'].apply(lambda x : x.upper())
            dflankfa.loc[:,'sequence: length']=[len(s) for s in dflankfa['sequence']]
            dseq=dbed.set_index('id').join(dflankfa,rsuffix='.1')

            dseq2compatible={'aminoacid: position':'aminoacid: position',
             'gene: id':'gene: id',
             'gene: name':'gene: name',
             'protein: id':'protein: id',
             'transcript: id':'seqid',
             'aminoacid: wild-type':'reference residue',
             'codon: wild-type':'reference codon',
             'contig':'contig',
             'strand':'strand',
             'start':'start',
             'end':'end',
             'codon start':'codon start',
             'codon end':'codon end',
            }
            dseq=dseq[list(dseq2compatible.values())]
            dseq.columns=list(dseq2compatible.keys())
#             dseq.to_csv('data/dseq.csv')            
            
        else:
            dseq=get_seq_yeast(din,
                          orfs_fastap='{}/../data/yeast/orf_coding_all.fasta'.format(abspath(dirname(__file__))),
                          host=cfg['host'],
                          test=cfg['test'])
        din.to_csv('{}/din.csv'.format(cfg['datad']))
        dseq.to_csv(dseqp)
        logging.info('Counts of amino acids to mutate:')
        logging.info(dseq['aminoacid: wild-type'].value_counts())
