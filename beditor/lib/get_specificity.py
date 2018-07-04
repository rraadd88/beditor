#!usr/bin/python

import logging
import subprocess
import re
import sys
import logging 
import operator

from collections import defaultdict
from os.path import join, basename, dirname, abspath, exists
from os import makedirs

import pandas as pd
from Bio import SeqIO

import pysam
import numpy as np
from glob import glob

from tqdm import trange

from beditor.lib.io_sys import runbashcmd
from beditor.lib.io_seqs import fa2df 
from beditor.lib.io_dfs import set_index,del_Unnamed 

def str2num(x):
    """
    This extracts numbers from strings. eg. 114 from M114R.
    :param x: string
    """
    return int(''.join(ele for ele in x if ele.isdigit()))

def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
#     print(s1,s2)
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1.upper(), s2.upper()))
def align(s1,s2,test=False):
    """
    Creates pairwise local alignment between seqeunces.
    Get the visualization and alignment scores.
    :param s1: seqeunce 1
    :param s2: seqeunce 2    
    """
    from Bio import pairwise2
    alignments = pairwise2.align.localms(s1.upper(),s2.upper(),1,-1,-5,-5)
    if test:
        print(alignments)
    alignsymb=np.nan
    score=np.nan
    sorted_alignments = sorted(alignments, key=operator.itemgetter(2))
    for a in alignments:
        alignstr=pairwise2.format_alignment(*a)
        alignsymb=alignstr.split('\n')[1]
        score=a[2]
        if test:
            print(alignstr)
        break
    return alignsymb,score
def gffatributes2ids(s):
    """
    Deconvolutes ids from `attributes` column in GFF3 to seprate columns.
    :param s: attribute string.
    :returns: tuple of ids
    """
    Name,gene_id,transcript_id,protein_id,exon_id=np.nan,np.nan,np.nan,np.nan,np.nan
    if '=' in s:
        d=dict([i.split('=') for i in s.split(';')])
        if 'Parent' in d:
            d[d['Parent'].split(':')[0]+'_id']=d['Parent'].split(':')[1]
        Name,gene_id,transcript_id,protein_id,exon_id=np.nan,np.nan,np.nan,np.nan,np.nan
        if 'Name' in d:    
            Name=d['Name']
        if 'gene_id' in d:    
            gene_id=d['gene_id']
        if 'transcript_id' in d:    
            transcript_id=d['transcript_id']
        if 'protein_id' in d:    
            protein_id=d['protein_id']
        if 'exon_id' in d:    
            exon_id=d['exon_id']
    return Name,gene_id,transcript_id,protein_id,exon_id
#--------------------------------------------
from beditor.lib.io_dfs import lambda2cols

def dguides2offtargets(cfg):
    """
    All the processes in offtarget detection are here.
    :param cfg: Configuration settings provided in .yml file
    """
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    
    stepn='04_offtargets'
    logging.info(stepn)
    dguidesp='{}/dguideslin.csv'.format(cfg[cfg['step']-1])
    datatmpd='{}/tmp'.format(cfg['datad'],stepn)
    for dp in [datatmpd]:
        if not exists(dp):
            makedirs(dp)

    dguides=pd.read_csv(dguidesp)
#     dguides.to_csv('{}/{}'.format(cfg['datad'],basename(dguidesp)),sep='\t')
    dguides=dguides.set_index('guide: id')

    guidesfap = '{}/guides.fa'.format(datatmpd)
    logging.info(basename(guidesfap))
    if not exists(guidesfap) or cfg['force']:
        with open(guidesfap,'w') as f:
            for gi in dguides.index:
                f.write('>{}\n{}\n'.format(gi.replace(' ','_'),dguides.loc[gi,'guide sequence+PAM']))

    # BWA: allow up to X mismatches
    # cfg['cores']=8

    # maximum number of occurences in the genome to get flagged as repeats. 
    # This is used in bwa samse, when converting the sam file
    # and for warnings in the table output.
    MAXOCC = 60000

    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000/MAXOCC

    #FIXME prepare genome
    # bwa index pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa
    # samtools faidx pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa    
    # cut pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa faidx   

    genomep=cfg['genomep']
    genomed = dirname(genomep) # make var local, see below
    genomegffp=cfg['genomegffp']
    
    # increase MAXOCC if there is only a single query, but only in CGI mode
    bwaM = MFAC*MAXOCC # -m is queue size in bwa
    guidessap = '{}/guides.sa'.format(datatmpd)
    logging.info(basename(guidessap))
    if not exists(guidessap) or cfg['force']:
        # cmd = "bwa aln -t %(cfg['cores'])s -o 0 -m %(bwaM)s -n %(cfg['mismatches_max'])d -k %(cfg['mismatches_max'])d -N -l %(guidel)d %(genomep)s %(guidesfap)s > %(guidessap)s" % locals()
        cmd = "bwa aln -t {} -o 0 -m {} -n {} -k {} -N -l {} {} {} > {}".format(cfg['cores'],bwaM,cfg['mismatches_max'],cfg['mismatches_max'],cfg['guidel'],genomep,guidesfap,guidessap)
        runbashcmd(cmd)

    guidessamp = '{}/guides.sam'.format(datatmpd)
    logging.info(basename(guidessamp))        
    if not exists(guidessamp) or cfg['force']:
        cmd = "bwa samse -n %(MAXOCC)d %(genomep)s %(guidessap)s %(guidesfap)s > %(guidessamp)s" % locals()
        runbashcmd(cmd)
    
    #----make tables-----------
    from beditor.lib.global_vars import bed_colns,gff_colns    

    alignmentbedp='{}/alignment.bed'.format(datatmpd)
    dalignbedp='{}/dalignbed.tsv'.format(datatmpd)
    logging.info(basename(dalignbedp))
    if not exists(alignmentbedp) or cfg['force']:
        samfile=pysam.AlignmentFile(guidessamp, "rb")
        dalignbed=pd.DataFrame(columns=bed_colns)
        for read in samfile.fetch():
            algnids=[]
            algnids.append('{}|{}{}|{}|{}'.format(read.reference_name,
                     ('-' if read.is_reverse else '+'),read.positions[0],read.cigarstring,read.get_tag('NM')))
            if read.has_tag('XA'):
                algnids+=['|'.join(s.split(',')) for s in read.get_tag('XA').split(';') if len(s.split(','))>1]
            chroms=[]
            starts=[]
            ends=[]
            ids=algnids
            NMs=[]
            strands=[]    
            for a in algnids:
                strand=a.split('|')[1][0]
                chroms.append(a.split('|')[0])
                if strand=='+':
                    offset=0
                elif strand=='-':
                    offset=0                    
                starts.append(int(a.split('|')[1][1:])+offset)
                ends.append(int(a.split('|')[1][1:])+str2num(a.split('|')[2])+offset)
                NMs.append(a.split('|')[3])
                strands.append(strand)
                del strand,offset
            col2dalignbed={'chromosome':chroms,
                           'start':starts,
                           'end':ends,
                           'id':ids,
                           'NM':NMs,
                           'strand':strands}
        #     col2dalignbed=dict(zip(cols,[a.split('|')[0],a.split('|')[1],a.split('|')[2],a,a.split('|')[3],a.split('|')[4] for a in algnids]))
            dalignbed_=pd.DataFrame(col2dalignbed)
            dalignbed_['guide: id']=read.qname.replace('_',' ')
            dalignbed = dalignbed.append(dalignbed_,ignore_index=True)
        #     break
        samfile.close()

        # filter bad asssembly junk genmomes
#         from beditor.lib.global_vars import host2contigs
#         dalignbed=dalignbed.loc[dalignbed['chromosome'].isin(host2contigs[cfg['host']]),:]
        dalignbed.to_csv(dalignbedp,sep='\t')
        from beditor.lib.io_nums import str2numorstr
        dalignbed['chromosome']=dalignbed.apply(lambda x : str2numorstr(x['chromosome']),axis=1)
        dalignbed=dalignbed.sort_values(['chromosome','start','end'], ascending=[True, True, True])
        dalignbed.loc[:,bed_colns].to_csv(alignmentbedp,sep='\t',
                        header=False,index=False,
                        chunksize=5000)
    else:
        dalignbed=pd.read_csv(dalignbedp,sep='\t')
        dalignbed=dalignbed.drop([c for c in dalignbed if 'Unnamed' in c],axis=1)
        
    alignmentbedsortedp=alignmentbedp+'.sorted.bed'
    logging.info(basename(alignmentbedsortedp))
    if not exists(alignmentbedsortedp) or cfg['force']:
        cmd='bedtools sort -i {} > {}'.format(alignmentbedp,alignmentbedsortedp)
        runbashcmd(cmd)
    
    genomegffsortedp=genomegffp+'.sorted.gff3.gz'
    logging.info(basename(genomegffsortedp))
    if not exists(genomegffsortedp) or cfg['force']:    
        cmd='bedtools sort -i {} > {}'.format(genomegffp,genomegffsortedp)
        runbashcmd(cmd)
        
    annotationsbedp='{}/annotations.bed'.format(datatmpd)
    logging.info(basename(annotationsbedp))
    if not exists(annotationsbedp) or cfg['force']:    
        cmd='bedtools intersect -wa -wb -loj -a {} -b {} > {}'.format(alignmentbedsortedp,genomegffsortedp,annotationsbedp)
        runbashcmd(cmd)
#     if the error in human, use: `cut -f 1 data/alignment.bed.sorted.bed | sort| uniq -c | grep -v CHR | grep -v GL | grep -v KI`
    dalignbedguidesp='{}/dalignbedguides.tsv'.format(datatmpd)
    logging.info(basename(dalignbedguidesp))
    if not exists(dalignbedguidesp) or cfg['force']:
        dalignbed=set_index(dalignbed,'guide: id').join(dguides)
        dalignbed=set_index(dalignbed,'guide: id')
        dalignbed.to_csv(dalignbedguidesp,'\t')
    else:
        dalignbed=pd.read_csv(dalignbedguidesp,'\t')
        dalignbed=del_Unnamed(dalignbed)
        
    # dalignid2seq=pd.DataFrame(columns=['sequence'])
    # dalignid2seq.index.name='id'
    dalignedfastap='{}/dalignedfasta.tsv'.format(datatmpd)
    logging.info(basename(dalignedfastap))
    if not exists(dalignedfastap) or cfg['force']:
        alignedfastap='{}/alignment.fa'.format(datatmpd)
        if not exists(alignedfastap) or cfg['force']:
            cmd='bedtools getfasta -s -name -fi {} -bed {} -fo {}'.format(genomep,alignmentbedp,alignedfastap)
            runbashcmd(cmd)

        dalignedfasta=fa2df(alignedfastap)
        dalignedfasta.columns=['aligned sequence']        
#         dalignedfasta=dalignedfasta.loc[[False if np.unique(list(s.upper()))=='N' else True for s in dalignedfasta['aligned sequence']],:]
        dalignedfasta.index=[i.split('(')[0] for i in dalignedfasta.index] # for bedtools 2.27, the fasta header now has hanging (+) or (-)
        dalignedfasta.index.name='id'
        dalignedfasta.to_csv(dalignedfastap,sep='\t')
    else:
        dalignedfasta=pd.read_csv(dalignedfastap,sep='\t')        
        dalignedfasta=dalignedfasta.drop([c for c in dalignbed if 'Unnamed' in c],axis=1)
        dalignedfasta=del_Unnamed(dalignedfasta)
        
    dalignbedguidesseqp='{}/dalignbedguidesseq.tsv'.format(datatmpd)
    logging.info(basename(dalignbedguidesseqp))
    if not exists(dalignbedguidesseqp) or cfg['force']:        
#         from beditor.lib.io_dfs import df2info
#         if cfg['test']:
#             df2info(dalignbed)
#             df2info(dalignedfasta)
            
#         dalignedfasta['guide: id']=dalignedfasta.apply(lambda x : x['guide: id'].replace('_',' '),axis=1)
        dalignbed=set_index(dalignbed,'id').join(set_index(dalignedfasta,'id'))
        dalignbed.index.name='id'
        dalignbed=dalignbed.drop_duplicates()
        dalignbed.to_csv(dalignbedguidesseqp,sep='\t')
    else:
        dalignbed=pd.read_csv(dalignbedguidesseqp,sep='\t',low_memory=False)
        dalignbed=del_Unnamed(dalignbed)

    dalignbedstatsp='{}/dalignbedstats.tsv'.format(datatmpd)  
    logging.info(basename(dalignbedstatsp))
    if not exists(dalignbedstatsp) or cfg['force']:
#         dalignbed['Hamming distance']=dalignbed.apply(lambda x : hamming_distance(x['guide sequence+PAM'], x['aligned sequence']),axis=1)
        df=dalignbed.apply(lambda x: align(x['guide sequence+PAM'],x['aligned sequence']),
                           axis=1).apply(pd.Series)
        df.columns=['alignment','alignment: score']
        dalignbed=dalignbed.join(df)
        del df
        dalignbed.to_csv(dalignbedstatsp,sep='\t')
    else:
        dalignbed=pd.read_csv(dalignbedstatsp,sep='\t',low_memory=False)
        dalignbed=del_Unnamed(dalignbed)
            
    daannotp='{}/dannot.tsv'.format(datatmpd)  
    dannotsaggp='{}/dannotsagg.tsv'.format(datatmpd)  
    logging.info(basename(daannotp))
    if (not exists(daannotp)) or (not exists(dannotsaggp)) or cfg['force']:
        dannots=pd.read_csv('{}/annotations.bed'.format(datatmpd),sep='\t',
                   names=bed_colns+[c+' annotation' if c in set(bed_colns).intersection(gff_colns) else c for c in gff_colns ],
                           low_memory=False)
        dannots=del_Unnamed(dannots)

        dannots=dannots.set_index('id')
        dannots['annotations count']=1
        # separate ids from attribute columns
        dannots=lambda2cols(dannots,lambdaf=gffatributes2ids,
                    in_coln='attributes',
                to_colns=['gene name','gene id','transcript id','protein id','exon id'])

        dannots['annotation coordinate']=dannots.apply(lambda x: '{}:{}-{}({})'.format(x['chromosome annotation'],x['start annotation'],x['end annotation'], x['strand annotation']),axis=1)
        dannots.to_csv(daannotp,sep='\t')
    else:
        dannots=pd.read_csv(daannotp,sep='\t',low_memory=False)
        dannots=del_Unnamed(dannots)

    logging.info(basename(dannotsaggp))
    if not exists(dannotsaggp) or cfg['force']:
        dannotsagg=pd.DataFrame(dannots.groupby('id')['annotations count'].agg('sum'))-1
        dannotsagg.loc[dannotsagg['annotations count']==0,'region']='intergenic'
        dannotsagg.loc[dannotsagg['annotations count']!=0,'region']='genic'

        dannots=dannots.reset_index()
        alignids=dannots['id'].unique()[:15]
        
        for alignidi in trange(len(alignids)):
            alignid=alignids[alignidi]
        #     if dannots.loc[i,'type'].value_counts().sum()==1:
            dannoti=dannots.loc[dannots['id']==alignid,:]
            if len(dannoti.shape)==1:
                dannoti=pd.DataFrame(dannoti).T
            dannoti=dannoti.loc[dannoti['type']!='chromosome',:].drop_duplicates(subset=['start annotation','end annotation'])
#             if cfg['test']:
# #                 print(alignid)
#                 print(dannoti.shape)
            for col in ['type','gene name','gene id','transcript id','protein id','exon id','annotation coordinate']:    
        #         print(";".join(dannoti[col].fillna('nan').tolist()))
                dannotsagg.loc[alignid,col+'s']=";".join(dannoti[col].fillna('nan').tolist())
        dannotsagg.to_csv(dannotsaggp,sep='\t')
    else:
        dannotsagg=pd.read_csv(dannotsaggp,sep='\t',low_memory=False)
        dannotsagg=del_Unnamed(dannotsagg)
        
    dalignbedannotp='{}/dalignbedannot.tsv'.format(datatmpd)  
    logging.info(basename(dalignbedannotp))
    if not exists(dalignbedannotp) or cfg['force']:
        dalignbedannot=set_index(dalignbed,'id').join(set_index(dannotsagg,'id'),
                                              rsuffix=' annotation')
        dalignbedannot.to_csv(dalignbedannotp,sep='\t')
    else:
        dalignbedannot=pd.read_csv(dalignbedannotp,sep='\t',low_memory=False)
        dalignbedannot=del_Unnamed(dalignbedannot)

    dofftargetsp='{}/dofftargets.tsv'.format(cfg['datad'])  
    logging.info(basename(dofftargetsp))
    if not exists(dofftargetsp) or cfg['force']:
        from beditor.lib.get_scores import get_beditorscore,get_CFDscore
        dalignbedannot['beditor score']=dalignbedannot.apply(lambda x : get_beditorscore(x['NM'], cfg['mismatches_max'], True if x['region']=='genic' else False, x['alignment']), axis=1) 
        dalignbedannot['CFD score']=dalignbedannot.apply(lambda x : get_CFDscore(x['guide sequence+PAM'].upper(), x['aligned sequence'].upper()), axis=1) 
        dalignbedannot.loc[:,['id',
 'guide: id',
#  'NM',
#  'chromosome',
#  'end',
#  'start',
#  'strand',
#  'Unnamed: 0',
#  'strategy',
 'guide sequence+PAM',
 'aligned sequence',
 'alignment',
#  'alignment: score',
 'annotations count',
 'region',
 'types',
#  'gene names',
#  'gene ids',
#  'transcript ids',
#  'protein ids',
#  'exon ids',
#  'annotation coordinates',
 'beditor score',
 'CFD score']].to_csv(dofftargetsp,sep='\t')
        
    # print('{}/dofftargets.tsv'.format(cfg['datad']))
    # dcombo.to_csv('{}/dofftargets.tsv'.format(cfg['datad']),sep='\t')


# if __name__ == '__main__':
#     main()
