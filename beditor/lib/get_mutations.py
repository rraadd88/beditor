import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os.path import exists,realpath,dirname
import itertools

from Bio import SeqIO, Alphabet, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

import logging

# local scripts

def get_codon_table(aa, tax_id=None):
    """
    Gets host specific codon table.

    :param aa: list of amino acids
    :param host: name of host
    :returns: codon table (pandas dataframe)
    """
    from Bio import Data
    # get codon table
    if tax_id is None:
        codontable=Data.CodonTable.unambiguous_dna_by_name["Standard"]
    else:
        codontable=Data.CodonTable.unambiguous_dna_by_id[tax_id]

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
    """
    Creates codon usage table.

    :param cuspp: path to cusp generated file
    :returns: codon usage table (pandas dataframe)
    """
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
    """
    Assesses possible mutagenesis strategies, given the set of Base editors and positions of mutations.

    :param dcodontable: Codon table
    :param dcodonusage: Codon usage table
    :param BEs: Base editors (dict), see global_vars.py
    :param pos_muts: positions of mutations
    :param host: host organism
    :returns: possible mutagenesis strategies as a pandas dataframe
    """
    def write_dmutagenesis(cdni,posi,codon,codonmut,ntwt,ntmut,aa,aamut,method):
        """
        Write dmutagenesis table for each iteraction in get_possible_mutagenesis.
        """
        dmutagenesis.loc[cdni,'codon']=codon
        dmutagenesis.loc[cdni,'position of mutation in codon']=int(posi)
        dmutagenesis.loc[cdni,'codon mutation']=codonmut
        dmutagenesis.loc[cdni,'nucleotide']=ntwt
        dmutagenesis.loc[cdni,'nucleotide mutation']=ntmut
        dmutagenesis.loc[cdni,'amino acid']=aa
        dmutagenesis.loc[cdni,'amino acid mutation']=aamut
        dmutagenesis.loc[cdni,'mutation on strand']=method.split(' on ')[1]
        dmutagenesis.loc[cdni,'strand: mutation']=method.split(' on ')[1].replace(' strand','')
        dmutagenesis.loc[cdni,'method']=method.split(' on ')[0]                        
        dmutagenesis.loc[cdni,'codon mutation usage Fraction']=dcodonusage.loc[codonmut,'Fraction']
        dmutagenesis.loc[cdni,'codon mutation usage Frequency']=dcodonusage.loc[codonmut,'Frequency']
        return dmutagenesis

    def get_sm(dmutagenesis,BEs,positions,codon,muti,cdni):
        """
        Fetches single nucleotide mutagenesis strategies.
        """
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
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=1))
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
        return dmutagenesis,muti
    def get_dm(dmutagenesis,BEs,positions_dm,codon,muti,cdni):
        """
        Fetches double nucleotide mutagenesis strategies.
        """
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
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=1))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt='{}{}'.format(BEs[method][0],BEs[method][0])
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
        return dmutagenesis,muti

    def get_tm(dmutagenesis,BEs,positions_tm,codon,muti,cdni):
        """
        Fetches triple nucleotide mutagenesis strategies.
        """
        for method in BEs:
            for posi1,posi2,posi3 in positions_tm:
                if (BEs[method][0]==codon[posi1]) and (BEs[method][0]==codon[posi2]) and (BEs[method][0]==codon[posi3]):
                    for ntmut1,ntmut2,ntmut3 in itertools.product(''.join(BEs[method][1]), repeat=3):
                        codonmut='{}{}{}'.format(ntmut1,ntmut2,ntmut3)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=1))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt='{}{}{}'.format(BEs[method][0],BEs[method][0],BEs[method][0])
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
        return dmutagenesis,muti

    def get_dm_combo(dmutagenesis,BEs,positions_dm,codon,muti,cdni,method):
        """
        Fetches double nucleotide mutagenesis strategies utilising 2 different base editors simultaneously.
        """
        methods=[m for m in itertools.product(BEs.keys(),repeat=2) if ((m[0].split('on')[1]==m[1].split('on')[1])) and (m[0]!=m[1])]
        for method1,method2 in methods:
            for posi1,posi2 in positions_dm: 
                if (BEs[method1][0]==codon[posi1]) and (BEs[method2][0]==codon[posi2]):
                    ntmuts=[(n1,n2) for n1 in ''.join(BEs[method1][1]) for n2 in ''.join(BEs[method2][1])]
                    for ntmut1,ntmut2 in ntmuts:
                        if (posi1==0) and (posi2==1):
                            codonmut='{}{}{}'.format(ntmut1,ntmut2,codon[2])
                        elif (posi1==1) and (posi2==2):
                            codonmut='{}{}{}'.format(codon[0],ntmut1,ntmut2)
                        elif (posi1==0) and (posi2==2):
                            codonmut='{}{}{}'.format(ntmut1,codon[1],ntmut2)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=1))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt='{}{}'.format(BEs[method1][0],BEs[method2][0])
                        ntmut='{}{}'.format(ntmut1,ntmut2)
                        if '-' in method1.split(' on ')[1]:
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
                        'method':method+' on '+method1.split('on')[1]})
        return dmutagenesis,muti

    def get_tm_combo(dmutagenesis,BEs,positions_tm,codon,muti,cdni,method):
        """
        Fetches triple nucleotide mutagenesis strategies utilising 2 different base editors simultaneously.
        """
        methods=[m for m in itertools.product(BEs.keys(),repeat=3) if ((m[0].split('on')[1]==m[1].split('on')[1]==m[2].split('on')[1])) and not (m[0]==m[1]==m[2])]
        for method1,method2,method3 in methods:
            for posi1,posi2,posi3 in positions_tm:
                if (BEs[method1][0]==codon[posi1]) and (BEs[method2][0]==codon[posi2]) and (BEs[method3][0]==codon[posi3]):
                    ntmuts=[(n1,n2,n3) for n1 in ''.join(BEs[method1][1]) for n2 in ''.join(BEs[method2][1]) for n3 in ''.join(BEs[method3][1])]
                    for ntmut1,ntmut2,ntmut3 in ntmuts:
                        codonmut='{}{}{}'.format(ntmut1,ntmut2,ntmut3)
                        aamut=str(Seq.Seq(codonmut,Alphabet.generic_dna).translate(table=1))
                        # if (aamut!='*') and (aamut!=aa): #  nonsence and synonymous
                        if muti==0:
                            cdni=cdni
                        else:
                            cdni=len(dmutagenesis)+1
                        muti+=1
                        ntwt='{}{}{}'.format(BEs[method1][0],BEs[method2][0],BEs[method3][0])
                        ntmut='{}{}{}'.format(ntmut1,ntmut2,ntmut3)
                        if '-' in method1.split(' on ')[1]:
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
                        'method':method+' on '+method1.split('on')[1]})
        return dmutagenesis,muti

    #double nucleotide mutations
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
        dmutagenesis,muti=get_sm(dmutagenesis,BEs,positions,codon,muti,cdni)
        #double nucleotide mutations
        dmutagenesis,muti=get_dm(dmutagenesis,BEs,positions_dm,codon,muti,cdni)
        #triple nucleotide mutations
        dmutagenesis,muti=get_tm(dmutagenesis,BEs,positions_tm,codon,muti,cdni)
        # #double nucleotide mutations combinations
        # dmutagenesis,muti=get_dm_combo(dmutagenesis,BEs,positions_dm,codon,muti,cdni,method='undefined')
        # #triple nucleotide mutations combinations
        # dmutagenesis,muti=get_tm_combo(dmutagenesis,BEs,positions_tm,codon,muti,cdni,method='undefined')

    dmutagenesis['nucleotide mutation: count']=[len(s) for s in dmutagenesis['nucleotide mutation']]
    dmutagenesis=dmutagenesis.sort_values('codon')  
    # Adding information of Allowed activity window
    dmutagenesis=dmutagenesis.set_index('method').join(pos_muts)
    dmutagenesis=dmutagenesis.reset_index()
    
    from beditor.lib.io_seqs import reverse_complement_multintseq
    from beditor.lib.global_vars import nt2complement
    dmutagenesis['nucleotide: wild-type']=dmutagenesis.apply(lambda x : x['nucleotide'] if x['strand: mutation']=='+' else reverse_complement_multintseq(x['nucleotide'],nt2complement),axis=1) 
    dmutagenesis['nucleotide: mutation']=dmutagenesis.apply(lambda x : x['nucleotide mutation'] if x['strand: mutation']=='+' else reverse_complement_multintseq(x['nucleotide mutation'],nt2complement),axis=1)    
    
    return dmutagenesis

from beditor.lib.io_dfs import df2unstack
from os.path import abspath,dirname
def get_submap(cfg):
    """
    Fetches mimetic substitution map that would be used to filter mutagenesis strategies.

    :param cfg: configurations from yml file.
    """
    mimetism_levels={'high': 1,
                     'medium': 5,
                     'low': 10}
    if not cfg['host'] in ['saccharomyces_cerevisiae','homo_sapiens']:
        host='saccharomyces_cerevisiae'
        logging.warning("for mimetic substitutions, substitution matrix of {host} is used")
    else:
        host=cfg['host']
    try:
        dsubmap=pd.read_csv(f'{dirname(realpath(__file__))}/../data/dsubmap_{host}.csv').set_index('AA1')
    except:
        if cfg['test']:
            print(f"{dirname(realpath(__file__))}/data/dsubmap_{host}.csv")
        dsubmap=pd.read_csv(f'data/dsubmap_{host}.csv').set_index('AA1')
        
    dsubmap.index.name='amino acid'
    dsubmap.columns.name='amino acid mutation'
    dsubmap=dsubmap.T

    dsubmaptop=pd.DataFrame(columns=dsubmap.columns,index=dsubmap.index)
    dsubmaptop=dsubmaptop.fillna(False)
    for i in dsubmap.index:
        for c in dsubmap:
            if i==c:
                dsubmaptop.loc[i,c]=True

    dsubmaptop_=dsubmaptop.copy()
    for c in dsubmap:
        dsubmaptop.loc[dsubmap.nlargest(mimetism_levels[cfg['mimetism_level']],c).index,c]=True

    dsubmaptop=df2unstack(dsubmaptop,col='mimetic')
    dsubmaptop=dsubmaptop[dsubmaptop['mimetic']]
    return dsubmaptop

def filterdmutagenesis(dmutagenesis,cfg):
    """
    Filters the mutagenesis strategies by multiple options provided in configuration file (.yml).

    :param dmutagenesis: mutagenesis strategies (pd.DataFrame)
    :param cfg: configurations from yml file
    """
    logging.info('filtering: dmutagenesis.shape: '+str(dmutagenesis.shape))
    # filter by mutation_type
    if 'mutation_type' in cfg:
        if not cfg['mutation_type'] is None:
            if cfg['mutation_type']=='S':
                dmutagenesis=dmutagenesis.loc[(dmutagenesis['amino acid']==dmutagenesis['amino acid mutation'])]
            elif cfg['mutation_type']=='N':
                dmutagenesis=dmutagenesis.loc[(dmutagenesis['amino acid']!=dmutagenesis['amino acid mutation'])]
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    
    # filter by nonsense
    if 'keep_mutation_nonsense' in cfg:
        if not cfg['keep_mutation_nonsense'] is None:
            if not cfg['keep_mutation_nonsense']:
                dmutagenesis=dmutagenesis.loc[(dmutagenesis['amino acid mutation']!='*'),:]
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    

    # filter by mutation per codon
    if 'max_subs_per_codon' in cfg:
        if not cfg['max_subs_per_codon'] is None:
            dmutagenesis=dmutagenesis.loc[(dmutagenesis['nucleotide mutation: count']==cfg['max_subs_per_codon']),:]
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    

    # filter by method
    if 'BEs' in cfg:
        if not cfg['BEs'] is None:
            dmutagenesis=dmutagenesis.loc[dmutagenesis['method'].isin(cfg['BEs']),:]
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    

    # filter by submap
    if 'mutations' in cfg:
        if (cfg['mutations']=='mimetic') or (cfg['mutations']=='substitutions'):
            if cfg['mutations']=='mimetic':                
                dsubmap=get_submap(cfg)
            elif cfg['mutations']=='substitutions':
                dsubmap=pd.read_csv(cfg['dsubmap_preferred_path'],sep='\t') # has two cols: amino acid and amino acid mutation

            import seaborn as sns
            dsubmap.to_csv(f"{cfg['datad']}/dsubmap.tsv",sep='\t')
            dmutagenesis=pd.merge(dsubmap,dmutagenesis,on=['amino acid','amino acid mutation'],how='inner')

            dsubmap['count']=1
            sns.heatmap(dsubmap.pivot_table(columns='amino acid',index='amino acid mutation',values='count'),square=True)
            plt.xlabel('wild-type amino acid')
            plt.ylabel('mutated amino acid')
            plt.savefig(f"{cfg['datad']}/heatmap_submap.svg")

            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    

    # filter non interchageables
    if 'non_intermutables' in cfg:
        if not cfg['non_intermutables'] is None:
            non_intermutables=list(itertools.permutations(''.join(cfg['non_intermutables']),2))
            dmutagenesis.loc[(dmutagenesis.apply(lambda row: not (row['amino acid'], 
                                                  row['amino acid mutation']) in non_intermutables, axis=1)),:]    
            logging.info('dmutagenesis.shape: '+str(dmutagenesis.shape))    
    return dmutagenesis

def dseq2dmutagenesis(cfg):
    # from beditor.lib.global_vars import BEs2mutations,pos_muts
    """
    Generates mutagenesis strategies from identities of reference and mutated codons (from dseq).
    
    :param cfg: configurations from yml file  
    """
    from .global_vars import stepi2cols

    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    dmutagenesisp=f"{cfg['datad']}/dmutagenesis.tsv"
    dmutagenesisallp=f"{cfg['datad']}/dmutagenesis_all.tsv"
    if not exists(dmutagenesisp) or cfg['force']:
        dseq=pd.read_csv('{}/dsequences.tsv'.format(cfg[cfg['step']-1]),sep='\t')
        if cfg['mutation_format']=='nucleotide':
            from .global_vars import aminoacids as aas
        elif cfg['mutation_format']=='aminoacid':
            if 'reverse_mutations' in cfg:
                if cfg['reverse_mutations']:
                    from .global_vars import aminoacids as aas
                else:
                    aas=list(dseq['aminoacid: wild-type'].unique())#['S','T','Y']
            else:
                aas=list(dseq['aminoacid: wild-type'].unique())#['S','T','Y']
        else:
            aas=list(dseq['aminoacid: wild-type'].unique())#['S','T','Y']

        dcodontable=get_codon_table(aa=aas)
        dcodonusage=get_codon_usage(cuspp='{}/../data/64_1_1_all_nuclear.cusp.txt'.format(abspath(dirname(__file__)))) #FIXME if prokaryote is used ?

        #create BEs and pos_muts for back-compatibility
        dBEs=pd.read_table(f"{dirname(realpath(__file__))}/../data/dBEs.tsv")
        dBEs=dBEs.loc[dBEs['method'].isin(cfg['BEs']),:]
        
        BEs2mutations={}
        for method in dBEs['method'].unique():
            for strand in dBEs['strand'].unique():
                dBEsi=dBEs.loc[(dBEs['method']==method) & (dBEs['strand']==strand),:]
                BEs2mutations[f"{method} on {strand} strand"]=[dBEsi['nucleotide'].unique().tolist()[0],
                                                               dBEsi['nucleotide mutation'].unique().tolist()]
        pos_muts=dBEs.loc[:,['method']+['distance of mutation from PAM: minimum',
         'distance of mutation from PAM: maximum',
         'distance of codon start from PAM: minimum',
         'distance of codon start from PAM: maximum']].drop_duplicates().set_index('method')

        dmutagenesis=get_possible_mutagenesis(dcodontable=dcodontable,dcodonusage=dcodonusage,
                                    BEs=BEs2mutations,pos_muts=pos_muts,
                                    host=cfg['host'])
        dmutagenesis.to_csv(dmutagenesisallp,sep='\t')
        
        dmutagenesis=filterdmutagenesis(dmutagenesis,cfg)            
        colns_pos=[c for c in dmutagenesis if ('position' in c) or ('Position' in c)]
        dmutagenesis.loc[:,colns_pos]=dmutagenesis.loc[:,colns_pos].astype(int)
        
        dmutagenesis.to_csv(dmutagenesisp,sep='\t')

        logging.info('Possible 1 nucleotide mutations:')
        logging.info(dmutagenesis.set_index('amino acid')[['amino acid mutation','method','codon','codon mutation',
        #               'position of mutation in codon','mutation on strand',
        #               'nucleotide','nucleotide mutation',
                     ]])
        for aa in aas:
            logging.info(aa+' can be mutated to:')
            logging.info(list(dmutagenesis.loc[dmutagenesis.loc[:,'amino acid']==aa,:].loc[:,'amino acid mutation'].unique()))

        import gc
        gc.collect()
