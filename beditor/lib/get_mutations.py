import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists,abspath,dirname
import itertools

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

import logging

# local scripts
from beditor.lib.global_vars  import hosts

def get_codon_table(aa, host):
    # get codon table
    codontable=Data.CodonTable.unambiguous_dna_by_id[hosts[host]]

    dcodontable=pd.DataFrame(pd.Series(codontable.forward_table))

    dcodontable.index.name='codon'
    dcodontable.columns=['amino acid']

    for cdn in codontable.stop_codons:
        dcodontable.loc[cdn,'amino acid']='*'        

    dcodontable=dcodontable.reset_index()
    rows=[]
    if isinstance(aa,list):
        for s in dcodontable['amino acid'].tolist():
            if s in aa:
                rows.append(True)
            else:
                rows.append(False)
    else:
        rows=dcodontable['amino acid']==aa
#     print(sum(rows))
    dcodontable=dcodontable.loc[rows,:].set_index('codon').reset_index()
    return dcodontable

def get_codon_usage(cuspp):
    # get codon usage stats
    dcodonusage=pd.read_csv(cuspp,sep='\t',header=5)
    cols=''.join(dcodonusage.columns.tolist()).split(' ')
    dcodonusage.columns=[cols[-1]]
    dcodonusage.index.names=cols[:-1]

    dcodonusage=dcodonusage.reset_index().set_index('Codon')
    dcodonusage['amino acid']=[SeqUtils.seq1(s) for s in dcodonusage['#AA']]
    return dcodonusage


def get_possible_mutagenesis(dcodontable,dcodonusage,
                             BEs,pos_muts,
                             host,
                            ): 
    def write_dmutagenesis(cdni,posi,codon,codonmut,ntwt,ntmut,aa,aamut,method):
        dmutagenesis.loc[cdni,'codon']=codon
        dmutagenesis.loc[cdni,'position of mutation in codon']=int(posi)
        dmutagenesis.loc[cdni,'codon mutation']=codonmut
        dmutagenesis.loc[cdni,'nucleotide']=ntwt
        dmutagenesis.loc[cdni,'nucleotide mutation']=ntmut
        dmutagenesis.loc[cdni,'amino acid']=aa
        dmutagenesis.loc[cdni,'amino acid mutation']=aamut
        dmutagenesis.loc[cdni,'mutation on strand']=method.split(' on ')[1]
        dmutagenesis.loc[cdni,'method']=method.split(' on ')[0]                        
        dmutagenesis.loc[cdni,'codon mutation usage Fraction']=dcodonusage.loc[codonmut,'Fraction']
        dmutagenesis.loc[cdni,'codon mutation usage Frequency']=dcodonusage.loc[codonmut,'Frequency']
        return dmutagenesis

    positions={0:'@1st position',1:'@2nd position',2:'@3rd position'}
    #double nucleotide mutations
    positions_dm=[(i,j)  for i in positions.keys() for j in positions.keys() if i<j]
    #double nucleotide mutations
    positions_tm=[[0,1,2]]

    dmutagenesis=dcodontable.copy()
    # test=True
    test=False
    for cdni in dmutagenesis.index:
        codon=dmutagenesis.loc[cdni,'codon']
        aa=dmutagenesis.loc[cdni,'amino acid']
        muti=0
        if test:
            print(codon)
        #single nucleuotide mutations
        for method in BEs:
            for posi in positions: 
                if BEs[method][0]==codon[posi]:
                    for ntmut in BEs[method][1]:
                        if posi==0:
                            codonmut='{}{}{}'.format(ntmut,codon[1],codon[2])
                        elif posi==1:
                            codonmut='{}{}{}'.format(codon[0],ntmut,codon[2])
                        elif posi==2:
                            codonmut='{}{}{}'.format(codon[0],codon[1],ntmut)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=hosts[host]))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt=BEs[method][0]
                        if '-' in method.split(' on ')[1]:
                            ntwt=str(Seq.Seq(ntwt,Alphabet.generic_dna).reverse_complement())
                            ntmut=str(Seq.Seq(ntmut,Alphabet.generic_dna).reverse_complement())
                        dmutagenesis=write_dmutagenesis(**{'cdni':cdni,
                        'posi':posi+1,
                        'codon':codon,
                        'codonmut':codonmut,
                        'ntwt':ntwt,
                        'ntmut':ntmut,
                        'aa':aa,
                        'aamut':aamut,
                        'method':method})
        #double nucleotide mutations
        for method in BEs:
            for posi1,posi2 in positions_dm: 
                if (BEs[method][0]==codon[posi1]) and (BEs[method][0]==codon[posi2]):
                    for ntmut1,ntmut2 in itertools.product(''.join(BEs[method][1]), repeat=2):
                        if (posi1==0) and (posi2==1):
                            codonmut='{}{}{}'.format(ntmut1,ntmut2,codon[2])
                        elif (posi1==1) and (posi2==2):
                            codonmut='{}{}{}'.format(codon[0],ntmut1,ntmut2)
                        elif (posi1==0) and (posi2==2):
                            codonmut='{}{}{}'.format(ntmut1,codon[1],ntmut2)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=hosts[host]))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt=BEs[method][0]
                        ntmut='{}{}'.format(ntmut1,ntmut2)
                        if '-' in method.split(' on ')[1]:
                            ntwt=str(Seq.Seq(ntwt,Alphabet.generic_dna).reverse_complement())
                            ntmut=str(Seq.Seq(ntmut,Alphabet.generic_dna).reverse_complement())
                        dmutagenesis=write_dmutagenesis(
                        **{'cdni':cdni,
                        'posi':'{}{}'.format(posi1,posi2),
                        'codon':codon,
                        'codonmut':codonmut,
                        'ntwt':ntwt,
                        'ntmut':ntmut,
                        'aa':aa,
                        'aamut':aamut,
                        'method':method})
        #triple nucleotide mutations
        for method in BEs:
            for posi1,posi2,posi3 in positions_tm:
                if (BEs[method][0]==codon[posi1]) and (BEs[method][0]==codon[posi2]) and (BEs[method][0]==codon[posi3]):
                    for ntmut1,ntmut2,ntmut3 in itertools.product(''.join(BEs[method][1]), repeat=3):
                        codonmut='{}{}{}'.format(ntmut1,ntmut2,ntmut3)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=hosts[host]))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt=BEs[method][0]
                        ntmut='{}{}{}'.format(ntmut1,ntmut2,ntmut3)
                        if '-' in method.split(' on ')[1]:
                            ntwt=str(Seq.Seq(ntwt,Alphabet.generic_dna).reverse_complement())
                            ntmut=str(Seq.Seq(ntmut,Alphabet.generic_dna).reverse_complement())
                        dmutagenesis=write_dmutagenesis(
                        **{'cdni':cdni,
                        'posi':'123',
                        'codon':codon,
                        'codonmut':codonmut,
                        'ntwt':ntwt,
                        'ntmut':ntmut,
                        'aa':aa,
                        'aamut':aamut,
                        'method':method})
                        
    dmutagenesis['nucleotide mutation: count']=[len(s) for s in dmutagenesis['nucleotide mutation']]
    dmutagenesis=dmutagenesis.sort_values('codon')  
    # Adding information of Allowed activity window
    dmutagenesis=dmutagenesis.set_index('method').join(pos_muts)
    dmutagenesis=dmutagenesis.reset_index()
    return dmutagenesis


from beditor.lib.global_vars import BEs,pos_muts
def dseq2dmutagenesis(cfg):
    dmutagenesisp='{}/dmutagenesis.csv'.format(cfg['datad'])
    if not exists(dmutagenesisp) or cfg['force']:
        dseq=pd.read_csv('{}/dseq.csv'.format(cfg['datad']))
        aas=list(dseq['aminoacid: wild-type'].unique())#['S','T','Y']

        dcodontable=get_codon_table(aa=aas, host=cfg['host'])

        dcodonusage=get_codon_usage(cuspp='{}/../data/64_1_1_all_nuclear.cusp.txt'.format(abspath(dirname(__file__))))

        dmutagenesis=get_possible_mutagenesis(dcodontable,dcodonusage,
        #                          aa=aas,
                                              BEs=BEs,pos_muts=pos_muts,
                                     host=cfg['host'])
        
        dmutagenesis.to_csv(dmutagenesisp)

        print('Possible 1 nucleotide mutations of the phospho-sites:')
        print(dmutagenesis.set_index('amino acid')[['amino acid mutation','method','codon','codon mutation',
        #               'position of mutation in codon','mutation on strand',
        #               'nucleotide','nucleotide mutation',
                     ]])
        for aa in aas:
            print(aa+' can be mutated to:')
            print(list(dmutagenesis.loc[dmutagenesis.loc[:,'amino acid']==aa,:].loc[:,'amino acid mutation'].unique()))