import pandas as pd
import numpy as np
from os.path import exists,abspath,dirname

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

#lib modules
import logging
from beditor.lib.io_dfs import set_index,df2info

import json

from beditor.lib.global_vars import bed_colns
from beditor.lib.io_seqs import fa2df 
from beditor.lib.io_sys import runbashcmd
from .io_dfs import del_Unnamed
from .global_vars import mutation_format2cols,stepi2cols


def tboundaries2positions(t):
    """
    Fetches positions from transcript boundaries.
    """
    coding_sequence_positions=[]
    for ini,end in t.coding_sequence_position_ranges:
        if t.strand == '+':
            coding_sequence_positions+=np.arange(ini, end+1,1).tolist()
        if t.strand == '-':
            coding_sequence_positions+=np.arange(end, ini-1,-1).tolist()            
    coding_sequence_positions+=np.arange(coding_sequence_positions[-1]+1, coding_sequence_positions[-1]+1+3,1).tolist()
    return coding_sequence_positions

def t2pmapper(t,coding_sequence_positions):
    """
    Maps transcript id with protein id. 
    """
    dcoding=pd.DataFrame(columns=['coding sequence positions','coding sequence'])
    # dcoding.index.name='transcript index'
    dcoding['coding sequence positions']=coding_sequence_positions
    dcoding['coding sequence']=list(t.coding_sequence)
    dcoding['protein sequence']=np.ravel(list(zip(list(t.protein_sequence+'*'),list(t.protein_sequence+'*'),list(t.protein_sequence+'*'))))
    prtidx=np.arange(1,len(t.protein_sequence+'*')+1)
    dcoding['protein index']=np.ravel(list(zip(prtidx,prtidx,prtidx)))
    dcoding['transcript index']=dcoding.index.values+1
    return dcoding.sort_values(by='transcript index',ascending=True)#.set_index('transcript index')

def get_seq_aminoacid(cfg,din):
    import pyensembl
    #import ensembl object that would fetch genes 
    # ensembl = pyensembl.EnsemblRelease(release=cfg['genomerelease'])
    ensembl = pyensembl.EnsemblRelease(species=pyensembl.species.Species.register(
    latin_name=cfg['host'],
    synonyms=[cfg['host']],
    reference_assemblies={
    cfg['genomeassembly']: (cfg['genomerelease'], cfg['genomerelease']),
    }),release=cfg['genomerelease'])
    
    din.index=range(len(din))
    dbedp='{}/dbedflank.bed'.format(cfg['datad'])
    dbed=pd.DataFrame(columns=bed_colns)
    terrpositions=[]
    terrnotfound=[]
    terrnoncoding=[]
    bedrowi=0
#             for i in trange(len(din)-1,desc='get positions for bedtools'):                
    for i in din.index:                
        if din.loc[i,'transcript: id'] in ensembl.transcript_ids():
            t=ensembl.transcript_by_id(din.loc[i,'transcript: id'])
            if t.is_protein_coding and t.contains_start_codon and t.contains_stop_codon:
                coding_sequence_positions=tboundaries2positions(t)
                if len(coding_sequence_positions)==len(t.coding_sequence):
                #TODO     need to check if the seq made from coding_sequence_positions is same as t.coding_seqeunce
                    dcoding=t2pmapper(t,coding_sequence_positions)
                    dcodingmutpos=dcoding.loc[(dcoding['protein index']==din.loc[i,'aminoacid: position']),:]
                    codon_positions=dcodingmutpos['coding sequence positions'].tolist()
                    if len(codon_positions)!=0:
                        dbed.loc[bedrowi,'chromosome']=t.contig
                        if cfg['test']:
                            print(din.loc[i,'transcript: id'],codon_positions)
                        if t.strand=='+':
                            dbed.loc[bedrowi,'codon start']=codon_positions[0]
                            dbed.loc[bedrowi,'codon end']=codon_positions[2]
                        elif t.strand=='-':
                            dbed.loc[bedrowi,'codon start']=codon_positions[2]
                            dbed.loc[bedrowi,'codon end']=codon_positions[0]
                        dbed.loc[bedrowi,'start']=dbed.loc[bedrowi,'codon start']-22 #FIXME put flank in the yml
                        dbed.loc[bedrowi,'end']=dbed.loc[bedrowi,'codon end']+21 #FIXME put flank in the yml

                        dbed.loc[bedrowi,'reference residue']=dcodingmutpos['protein sequence'].tolist()[0]
                        dbed.loc[bedrowi,'reference codon']=''.join(dcodingmutpos['coding sequence'].tolist())
                        dbed.loc[bedrowi,'strand']=t.strand
                        dbed.loc[bedrowi,'id']='{}|{}|{}|{}|{}'.format(din.loc[i,'transcript: id'],
                                                                      dbed.loc[bedrowi,'chromosome'],
                                 dbed.loc[bedrowi,'strand'],int(dbed.loc[bedrowi,'start']),int(dbed.loc[bedrowi,'end']))        
                        dbed.loc[bedrowi,'gene: id']=t.gene_id
                        dbed.loc[bedrowi,'gene: name']=t.gene.name
                        dbed.loc[bedrowi,'protein: id']=t.protein_id
                        dbed.loc[bedrowi,'aminoacid: position']=din.loc[i,'aminoacid: position']
                #         break
                        bedrowi+=1
                    else:
                        terrpositions.append(t.id)                        
                else:
                    terrpositions.append(t.id)
            else:
                terrnoncoding.append(t.id)                    
        else:
            terrnotfound.append(din.loc[i,'transcript: id'])
            if cfg['test']:
                logging.error('not found: {}'.format(din.loc[i,'transcript: id']))
    dbed=dbed.loc[(dbed.apply(lambda x : x['end']-x['start']==45, axis=1)),:] #FIXME put flank in the yml

    dbed.loc[:,'start']=dbed.loc[:,'start'].astype(int)
    dbed.loc[:,'end']=dbed.loc[:,'end'].astype(int)
    
    dbed.loc[:,bed_colns].to_csv(dbedp,sep='\t',
                    header=False,index=False)
    err2tids={'terrpositions':terrpositions,
              'terrnotfound':terrnotfound,
              'terrnoncoding':terrnoncoding,
             }
    if cfg['test']:
        print(err2tids)            
    with open(dbedp+'.err.json', 'w') as outfile:
        json.dump(err2tids, outfile)

    bedp=f"{cfg['datad']}/dbedflank.bed"
    fastap=f"{cfg['datad']}/dbedflank.fa"
    cmd=f"{cfg['bedtools']} getfasta -s -name -fi {cfg['genomep']} -bed {bedp} -fo {fastap}"
    runbashcmd(cmd)

    dflankfa=fa2df(fastap,ids2cols=True)                
    dflankfa.loc[:,'sequence']=dflankfa.loc[:,'sequence'].apply(lambda x : x.upper())
    dflankfa.loc[:,'sequence: length']=[len(s) for s in dflankfa['sequence']]
    dflankfa.index=[idx.split('(')[0] for idx in dflankfa.index]
    dflankfa.index.name='id'
    dseq=set_index(dbed,'id').join(set_index(dflankfa,'id'),rsuffix='.1')
    dseq2compatible={'aminoacid: position':'aminoacid: position',
     'gene: id':'gene: id',
     'gene: name':'gene: name',
     'protein: id':'protein: id',
     'transcript: id':'seqid',
     'transcript: sequence':'sequence',
     'aminoacid: wild-type':'reference residue',
     'codon: wild-type':'reference codon',
     'contig':'contig',
     'strand':'strand',
     'start':'start',
     'end':'end',
     'codon start':'codon start',
     'codon end':'codon end',
    }
    if 'amino acid mutation' in dseq:
        dseq2compatible['amino acid mutation']='amino acid mutation'
    dseq.to_csv(cfg['dseqtmpp'],sep='\t')
    
    dseq=dseq[list(dseq2compatible.values())]
    dseq.columns=list(dseq2compatible.keys())
#             dseq.to_csv('data/dseq.csv')            

    logging.info(dseq.columns.tolist())
    logging.info(din.columns.tolist())
    dseq=pd.merge(dseq.reset_index(),din,on=['transcript: id','aminoacid: position'])
    logging.info(dseq.columns.tolist())
    set_index(dseq,'id')
    if 'reverse_mutations' in cfg:
        if cfg['reverse_mutations']:
            from beditor.lib.io_dfs import dfswapcols
            dseq=dfswapcols(dseq,['aminoacid: wild-type','amino acid mutation'])
            dseq['codon: mutation']=dseq['codon: wild-type'].copy()
            
    dseq.to_csv(f"{cfg['dsequencesp']}",sep='\t')
    del ensembl

from .io_seqs import genomeocoords2bed,fa2df
from .global_vars import flankntc
def get_seq_nucleotide(cfg,din):
    bedp=f"{cfg['datad']}/dbedntmuts.bed"
    fastap=f"{cfg['datad']}/dbedntmuts.fa"
    dbedntmutsp=f"{cfg['datad']}/dbedntmuts.tsv"
    if not exists(cfg['dsequencesp']) or cfg['force']:
        if not exists(bedp) or cfg['force']:            
            dbed=genomeocoords2bed(din,col_genomeocoord='genome coordinate')
            dbed['start']=dbed['start'].astype(int)-flankntc-1
            dbed['end']=dbed['end'].astype(int)+flankntc
            dbed.to_csv(bedp,sep='\t',header=False, index=False)
        if not exists(fastap) or cfg['force']:
            cmd=f"{cfg['bedtools']} getfasta -s -name -fi {cfg['genomep']} -bed {bedp} -fo {fastap}"
            runbashcmd(cmd)
        if not exists(dbedntmutsp) or cfg['force']:
            dbedntmuts=fa2df(fastap)
            dbedntmuts.columns=['transcript: sequence']
            dbedntmuts['transcript: sequence']=dbedntmuts.apply(lambda x: x['transcript: sequence'].upper(),axis=1)
            dbedntmuts=dbedntmuts.reset_index()
            dbedntmuts['genome coordinate']=dbedntmuts.apply(lambda x : x['id'].split('(')[0] ,axis=1)
            dbedntmuts.to_csv(dbedntmutsp,sep='\t')
        else:
            dbedntmuts=pd.read_table(dbedntmutsp)
    dsequences=pd.merge(din,dbedntmuts,
            on=['genome coordinate'],suffixes=('', ': dbedntmuts'))
    dsequences=del_Unnamed(dsequences)
#     print(dsequences[['codon: wild-type']].head())
    col_nt_wt='nucleotide wild-type' if not 'nucleotide wild-type' in dsequences else 'nucleotide wild-type: from flanking sequence'    
    col_nt_mt='nucleotide mutation' if not 'nucleotide mutation' in dsequences else 'nucleotide mutation: from flanking sequence'    
    col_cd_wt='codon: wild-type' if not 'codon: wild-type' in dsequences else 'codon: wild-type: from flanking sequence'
    col_cd_mt='codon: mutation' if not 'codon: mutation' in dsequences else 'codon: mutation: from flanking sequence'        
    dsequences[col_nt_wt]=dsequences.apply(lambda x: x['transcript: sequence'][flankntc],axis=1)        
    
#     print(dsequences[['codon: wild-type']].head())
    dsequences[col_cd_wt]=dsequences.apply(lambda x: x['transcript: sequence'][flankntc-1:flankntc+2],axis=1)
#     print(dsequences[['codon: wild-type']].head())
            
    dsequences[col_cd_mt]=dsequences.apply(lambda x: f"{x['codon: wild-type'][0]}{x['nucleotide mutation']}{x['codon: wild-type'][2]}",axis=1)
#     if 'reverse_mutations' in cfg:
#         if cfg['reverse_mutations']:
#             dsequences[col_cd_mt]=dsequences.apply(lambda x: f"{x['codon: wild-type'][0]}{x['nucleotide wild-type']}{x['codon: wild-type'][2]}",axis=1)
            
#     print(dsequences[['codon: wild-type', 'codon: mutation']].head())
    dsequences['transcript: id']=dsequences['genome coordinate']
    dsequences_bedcols=genomeocoords2bed(dsequences, col_genomeocoord='genome coordinate')
    for col in dsequences_bedcols:
        dsequences[col]=dsequences_bedcols[col]
    if 'reverse_mutations' in cfg:
        if cfg['reverse_mutations']:
            from beditor.lib.io_dfs import dfswapcols
            dseq=dfswapcols(dsequences,['nucleotide wild-type', 'nucleotide mutation'])
            dseq=dfswapcols(dsequences,['codon: wild-type', 'codon: mutation'])
#     print(dsequences[['codon: wild-type', 'codon: mutation']].head())
    dsequences.to_csv(f"{cfg['dsequencesp']}",sep='\t')
    # return dsequences

def din2dseq(cfg):
    """
    Wrapper for converting input data (transcript ids and positions of mutation) to seqeunces flanking the codon. 
    """
    from beditor.lib.global_vars import stepi2cols_nucleotide
    
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    # get dna and protein sequences 
    dsequencesp=f"{cfg['datad']}/dsequences.tsv"
    dseqtmpp=f"{cfg['datad']}/dseqtmp.tsv"
    cfg['dsequencesp']=dsequencesp
    cfg['dseqtmpp']=dseqtmpp
    if not exists(dsequencesp) or cfg['force']:
        cfg['dinp']=f"{cfg[cfg['step']-1]}/dinput.tsv"
        din=pd.read_table(cfg['dinp'])
        din=del_Unnamed(din)
        if cfg['mutation_format']=='aminoacid':        
            cols_dseq=stepi2cols[cfg['step']]
            if cfg['reverse_mutations']:
                if not 'codon: mutation' in cols_dseq:
                    cols_dseq+=['codon: mutation']
        elif cfg['mutation_format']=='nucleotide':
            cols_dseq=stepi2cols_nucleotide[cfg['step']]
            din=din.dropna(subset=['genome coordinate'])
            if len(din)==0:
                din=pd.DataFrame(columns=cols_dseq)
                logging.warning('no genome coordinates; saving an empty table.')                
        ifdinpisdseq=all([True if c in din else False for c in cols_dseq])
        if not ifdinpisdseq:
            if cfg['mutation_format']=='aminoacid':
                if all([True for c in mutation_format2cols['aminoacid'] if c in din]):
                    get_seq_aminoacid(cfg,din)
            elif cfg['mutation_format']=='nucleotide':
                if all([True for c in mutation_format2cols['nucleotide'] if c in din]):
                    get_seq_nucleotide(cfg,din)
            else:
                raise(ValueError(f"invalid value of cfg['mutation_format']: {cfg['mutation_format']}"))
            
            import gc
            gc.collect()
        else:
            din.to_csv(dsequencesp,sep='\t')