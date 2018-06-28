#!usr/bin/python

import logging
import subprocess
import logging
import re
import sys
import logging, operator

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

def pamIsCpf1(pam):
    " if you change this, also change bin/filterFaToBed! "
    return (pam in ["TTN", "TTTN", "TYCV", "TATV"])
        
def setupPamInfo(pam,addGenePlasmids,scoreNames):
    " modify a few globals based on the current pam "
    PAMLEN = len(pam)
    if pamIsCpf1(pam):
        logging.debug("switching on Cpf1 mode, guide length is 23bp")
        guidel = 23
        cpf1Mode = True
        scoreNames = cpf1ScoreNames
    elif pam=="NNGRRT" or pam=="NNNRRT":
        logging.debug("switching on S. aureus mode, guide length is 21bp")
        addGenePlasmids = addGenePlasmidsAureus
        guidel = 21
        cpf1Mode = False
        scoreNames = saCas9ScoreNames
    else:
        guidel = 20
        cpf1Mode = False
    return guidel,cpf1Mode,addGenePlasmids,PAMLEN,scoreNames

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
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))
def align(s1,s2,l,test=False):
    from Bio import pairwise2
    alignments = pairwise2.align.globalms(s1,s2,1,-1,-5,-5)
#     globalms(sequenceA, sequenceB, match, mismatch, open, extend)
#     if len(alignments)==0:
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

#--------------------------------------------
def dguides2offtargets(cfg):
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    
    #TODO get fa and gff and prepare them

    # cfg['datad']='../../../04_offtarget/data_cfg['test']_'
    # dguidesp='../../../04_offtarget/crisporWebsite/sampleFiles/mine/cfg['test'].csv'
#     dguidesp='../../../04_offtarget/data_cfg['test']_human_dguides.tsv'
#     cfg['datad']='../../../04_offtarget/data_cfg['test']_human_sabatani_cfg['test']'
#     cfg['host']='homo_sapiens'
#     genomefn='dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa'
#     genomegffp='pub/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh38.92.gff3.gz'

    stepn='04_offtargets'
    dguidesp='{}/dguideslin.csv'.format(cfg[cfg['step']-1])
#     host_="_".join(s for s in cfg['host'].split('_')).capitalize()

#     genomed='pub/release-{}/fasta/'.format(cfg['genomerelease'])
#     genomefn_='dna/{}.{}.dna_sm.*.fa'.format(host_,cfg['genomeassembly'])
#     genomep=glob('{}/{}/{}'.format(genomed,cfg['host'],genomefn))[0]

#     genomeannotd='pub/release-{}/gff3/'.format(cfg['genomerelease'])
#     genomegffp='{}/{}/{}.{}.{}.gff3.gz'.format(genomeannotd,cfg['host'],host_,cfg['genomeassembly'],cfg['genomerelease'])

    datatmpd='{}/tmp'.format(cfg['datad'],stepn)
    for dp in [datatmpd]:
        if not exists(dp):
            makedirs(dp)

    dguides=pd.read_csv(dguidesp)
#     dguides.to_csv('{}/{}'.format(cfg['datad'],basename(dguidesp)),sep='\t')
    dguides=dguides.set_index('guide: id')

    batchId='align'
    batchBase = join(datatmpd, batchId)
    otBedFname = batchBase+".bed"
    print(otBedFname)
    faFname = batchBase+".fa"
    if not exists(faFname) or cfg['force']:
        with open(faFname,'w') as f:
            for gi in dguides.index:
                f.write('>{}\n{}\n'.format(gi.replace(' ','_'),dguides.loc[gi,'guide sequence+PAM']))

    # BWA: allow up to X mismatches
    maxMMs=5

    # maximum number of occurences in the genome to get flagged as repeats. 
    # This is used in bwa samse, when converting the same file
    # and for warnings in the table output.
    MAXOCC = 60000

    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000/MAXOCC
    guidel=None
    PAMLEN=None
    cpf1Mode=None
    scoreNames = ["fusi", "crisprScan"]
    cpf1ScoreNames = ["seqDeepCpf1"]
    saCas9ScoreNames = ["najm"]
    cfg['test']=True
    pam='NGG'

    #FIXME prepare genome
    # bwa index pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa
    # samtools faidx pub/release-92/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa    

    guidel,cpf1Mode,addGenePlasmids,PAMLEN,scoreNames=setupPamInfo(pam,setupPamInfo,scoreNames)

    # get sequence
    # seqs = parseFasta(open(inSeqFname))

    matchesBedFname = batchBase+".matches.bed"
    saFname = batchBase+".sa"
    samp = batchBase+".sam"
    genomep=cfg['genomep']
    genomed = dirname(genomep) # make var local, see below
    genomegffp=cfg['genomegffp']
    
    open(matchesBedFname, "w") # truncate to 0 size

    # increase MAXOCC if there is only a single query, but only in CGI mode
    maxDiff = maxMMs
    seqLen = guidel

    bwaM = MFAC*MAXOCC # -m is queue size in bwa
    if not exists(saFname) or cfg['force']:
        cmd = "bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomep)s %(faFname)s > %(saFname)s" % locals()
        runbashcmd(cmd)
        
    if not exists(samp) or cfg['force']:
        cmd = "bwa samse -n %(MAXOCC)d %(genomep)s %(saFname)s %(faFname)s > %(samp)s" % locals()
        runbashcmd(cmd)
    
    #----make tables-----------
    gff_colns = ['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    bed_colns = ['chromosome','start','end','id','NM','strand']
    
    alignmentbedp='{}/alignment.bed'.format(datatmpd)
    dalignbedp='{}/dalignbed.tsv'.format(datatmpd)
    if not exists(alignmentbedp) or cfg['force']:
        samfile=pysam.AlignmentFile(samp, "rb")
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
        from beditor.lib.global_vars import host2contigs
        dalignbed=dalignbed.loc[dalignbed['chromosome'].isin(host2contigs[cfg['host']]),:]
        dalignbed.to_csv(dalignbedp,sep='\t')
        dalignbed.loc[:,bed_colns].to_csv(alignmentbedp,sep='\t',
                        header=False,index=False)
    else:
        dalignbed=pd.read_csv(dalignbedp,sep='\t')
        dalignbed=dalignbed.drop([c for c in dalignbed if 'Unnamed' in c],axis=1)
        
    alignmentbedsortedp=alignmentbedp+'.sorted.bed'
    if not exists(alignmentbedsortedp) or cfg['force']:
        cmd='bedtools sort -i {} > {}'.format(alignmentbedp,alignmentbedsortedp)
        runbashcmd(cmd)
    
    genomegffsortedp=genomegffp+'.sorted.gff3.gz'
    if not exists(genomegffsortedp) or cfg['force']:    
        cmd='bedtools sort -i {} > {}'.format(genomegffp,genomegffsortedp)
        runbashcmd(cmd)
        
    annotationsbedp='{}/annotations.bed'.format(datatmpd)
    if not exists(annotationsbedp) or cfg['force']:    
        cmd='bedtools intersect -wa -wb -loj -a {} -b {} > {}'.format(alignmentbedsortedp,genomegffsortedp,annotationsbedp)
        runbashcmd(cmd)
#     if the error in human, use: `cut -f 1 data/alignment.bed.sorted.bed | sort| uniq -c | grep -v CHR | grep -v GL | grep -v KI`
    dalignbedguidesp='{}/dalignbedguides.tsv'.format(datatmpd)
    if not exists(dalignbedguidesp) or cfg['force']:
        dalignbed=dalignbed.set_index('guide: id').join(dguides)
        dalignbed=dalignbed.reset_index().set_index('guide: id')
        dalignbed.to_csv(dalignbedguidesp,'\t')
    else:
        dalignbed=pd.read_csv(dalignbedguidesp,'\t')
        dalignbed=dalignbed.drop([c for c in dalignbed if 'Unnamed' in c],axis=1)
        
    # dalignid2seq=pd.DataFrame(columns=['sequence'])
    # dalignid2seq.index.name='id'
    dalignedfastap='{}/dalignedfasta.tsv'.format(datatmpd)
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
        
    dalignbedguidesseqp='{}/dalignbedguidesseq.tsv'.format(datatmpd)
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
        dalignbed=pd.read_csv(dalignbedguidesseqp,sep='\t')
        dalignbed=dalignbed.drop([c for c in dalignbed if 'Unnamed' in c],axis=1)
                
    dalignbed['Hamming distance']=dalignbed.apply(lambda x : hamming_distance(x['guide sequence+PAM'], x['aligned sequence']),axis=1)
    dalignbed['alignment','alignment: score']=dalignbed.apply(lambda x: align(x['guide sequence+PAM'],x['aligned sequence'],l=guidel),axis=1)

    # for i in dalignbed.index:
    #     dalignbed.loc[i,'alignment'],dalignbed.loc[i,'alignment: score']=align(dalignbed.loc[i,'guide sequence+PAM'],dalignbed.loc[i,'aligned sequence'],l=guidel)
    dannots=pd.read_csv('{}/annotations.bed'.format(datatmpd),sep='\t',
               names=bed_colns+gff_colns)

    dcombo=dalignbed.set_index('id').join(dannots.reset_index().set_index('id'),rsuffix='.2')
    dcombo.to_csv('{}/dofftargets.tsv'.format(cfg['datad']),sep='\t')


# if __name__ == '__main__':
#     main()
