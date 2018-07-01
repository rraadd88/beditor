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
    from Bio import pairwise2
    alignments = pairwise2.align.globalms(s1.upper(),s2.upper(),1,-1,-5,-5)
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
    maxMMs=0
    cores=8

    # maximum number of occurences in the genome to get flagged as repeats. 
    # This is used in bwa samse, when converting the sam file
    # and for warnings in the table output.
    MAXOCC = 60000

    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000/MAXOCC
    guidel=23
    PAMLEN=3
    pam='NGG'

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
        cmd = "bwa aln -t %(cores)s -o 0 -m %(bwaM)s -n %(maxMMs)d -k %(maxMMs)d -N -l %(guidel)d %(genomep)s %(guidesfap)s > %(guidessap)s" % locals()
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
                chroms.append(a.split('|')[0])
                starts.append(int(a.split('|')[1][1:])-1)
                ends.append(int(a.split('|')[1][1:])+str2num(a.split('|')[2])-1)
                NMs.append(a.split('|')[3])
                strands.append(a.split('|')[1][0])
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
        dalignbed.loc[:,bed_colns].to_csv(alignmentbedp,sep='\t',
                        header=False,index=False)
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
        dalignbed=dalignbed.set_index('guide: id').join(dguides)
        dalignbed=dalignbed.reset_index().set_index('guide: id')
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
        dalignbed=dalignbed.reset_index().set_index('id').join(dalignedfasta)
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
            
    dalignbedannotp='{}/dalignbedannot.tsv'.format(datatmpd)  
    logging.info(basename(dalignbedannotp))
    if not exists(dalignbedannotp) or cfg['force']:
        dannots=pd.read_csv('{}/annotations.bed'.format(datatmpd),sep='\t',
                   names=bed_colns+gff_colns,
                           low_memory=False)
        dannots=del_Unnamed(dannots)

        dcombo=set_index(dalignbed,'id').join(set_index(dannots,'id'),
                                              rsuffix='.2')
        dcombo=dcombo.drop_duplicates(subset=['type','guide: id','start.1','end.1'])
        #separate ids from attribute columns
#         dcombo=lambda2cols(dcombo,lambdaf=gffatributes2ids,
#                     in_coln='attributes',
#                 to_colns=['gene name','gene id','transcript id','protein id','exon id'])
        dcombo.to_csv(dalignbedannotp,sep='\t')
    else:
        dcombo=pd.read_csv(dalignbedannotp,sep='\t',low_memory=False)
        dcombo=del_Unnamed(dcombo)
    print('{}/dofftargets.tsv'.format(cfg['datad']))
    dcombo.to_csv('{}/dofftargets.tsv'.format(cfg['datad']),sep='\t')


# if __name__ == '__main__':
#     main()
