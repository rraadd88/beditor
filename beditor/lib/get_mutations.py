import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

# local scripts
from global_vars import hosts

def get_codon_table(aa, host):
    # get codon table
    codontable=Data.CodonTable.unambiguous_dna_by_id[hosts[host]]

    dcodontable=pd.DataFrame(pd.Series(codontable.forward_table))

    dcodontable.index.name='codon'
    dcodontable.columns=['amino acid']

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

    positions={0:'@1st position',1:'@2nd position',2:'@3rd position'}
    dmutagenesis=dcodontable.copy()
    # test=True
    test=False
    for cdni in dmutagenesis.index:
        codon=dmutagenesis.loc[cdni,'codon']
        aa=dmutagenesis.loc[cdni,'amino acid']
        muti=0
        if test:
            print(codon)
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
                        if (aamut!='*') and (aamut!=aa):
                            if muti==0:
                                cdni=cdni
                            else:
                                cdni=len(dmutagenesis)+1
                            muti+=1
                            ntwt=BEs[method][0]
                            if '-' in method.split(' on ')[1]:
                                ntwt=str(Seq.Seq(ntwt,Alphabet.generic_dna).reverse_complement())
                                ntmut=str(Seq.Seq(ntmut,Alphabet.generic_dna).reverse_complement())
                            dmutagenesis.loc[cdni,'codon']=codon
                            dmutagenesis.loc[cdni,'position of mutation in codon']=posi+1
                            dmutagenesis.loc[cdni,'codon mutation']=codonmut
                            dmutagenesis.loc[cdni,'nucleotide']=ntwt
                            dmutagenesis.loc[cdni,'nucleotide mutation']=ntmut
                            dmutagenesis.loc[cdni,'amino acid']=aa
                            dmutagenesis.loc[cdni,'amino acid mutation']=aamut
                            dmutagenesis.loc[cdni,'mutation on strand']=method.split(' on ')[1]
                            dmutagenesis.loc[cdni,'method']=method.split(' on ')[0]                        
                            dmutagenesis.loc[cdni,'codon mutation usage Fraction']=dcodonusage.loc[codonmut,'Fraction']
                            dmutagenesis.loc[cdni,'codon mutation usage Frequency']=dcodonusage.loc[codonmut,'Frequency']
    dmutagenesis=dmutagenesis.sort_values('codon')  
    # Adding information of Allowed activity window
    dmutagenesis=dmutagenesis.set_index('method').join(pos_muts)
    dmutagenesis=dmutagenesis.reset_index()
    return dmutagenesis
