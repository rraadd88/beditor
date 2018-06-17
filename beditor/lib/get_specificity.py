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

from Bio import SeqIO

def pamIsCpf1(pam):
    " if you change this, also change bin/filterFaToBed! "
    return (pam in ["TTN", "TTTN", "TYCV", "TATV"])
        
def pamIsSaCas9(pam):
    " only used for notes and efficiency scores, unlike its Cpf1 cousin function "
    return (pam in ["NNGRRT", "NNNRRT"])

def pamIsXCas9(pam):
    " "
    return (pam in ["NGK", "NGN"])
    
def pamIsSpCas9(pam):
    " only used for notes and efficiency scores, unlike its Cpf1 cousin function "
    return (pam in ["NGG", "NGA", "NGCG"])

def setupPamInfo(pam,addGenePlasmids,scoreNames):
    " modify a few globals based on the current pam "
    PAMLEN = len(pam)
    if pamIsCpf1(pam):
        logging.debug("switching on Cpf1 mode, guide length is 23bp")
        GUIDELEN = 23
        cpf1Mode = True
        scoreNames = cpf1ScoreNames
    elif pam=="NNGRRT" or pam=="NNNRRT":
        logging.debug("switching on S. aureus mode, guide length is 21bp")
        addGenePlasmids = addGenePlasmidsAureus
        GUIDELEN = 21
        cpf1Mode = False
        scoreNames = saCas9ScoreNames
    else:
        GUIDELEN = 20
        cpf1Mode = False
    return GUIDELEN,cpf1Mode,addGenePlasmids,PAMLEN,scoreNames

def runCmd(cmd):
    from beditor.lib.global_vars import dirs2ps 
    cmd = cmd.replace("$BIN", dirs2ps['binDir'])
    cmd = cmd.replace("$PYTHON", dirs2ps['pyp'])
    cmd = cmd.replace("$SCRIPT", dirs2ps['scriptDir'])
#     print(cmd)
    err=subprocess.call(cmd,shell=True)
    if err!=0:
        print('bash command error: {}\n{}\n'.format(err,cmd))
        sys.exit(1)
           

#-new code-----------
import pysam
import numpy as np

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
def align(s1,s2):
    from Bio import pairwise2
    alignments = pairwise2.align.localxx(s1, s2)
#     if len(alignments)==0:
#     print(len(alignments))
    alignsymb=np.nan
    score=np.nan
    for a in alignments:
        if len(a[0])==GUIDELEN and len(a[1])==GUIDELEN:
            alignstr=pairwise2.format_alignment(*a)
            alignsymb=alignstr.split('\n')[1]
            score=a[2]
            break
    return alignsymb,score

#--------------------------------------------
# def main(inSeqFname,genomeDir,org,outGuideFname,offtargetFname,genomefn):
def main():
    force=True
    dguidesp='data_test_human_dguides.tsv'

    # inSeqFname='../../../04_offtarget/data/04_specificity/in/sample.sacCer3.fa',
    # faFname='tmp/in/x50sPGMoTvUagv3zWjGg.fa'
    # genomeDir='tmp/in/genomes/'
    datad='data_test_human'
    dataind='{}/04_specificity/in'.format(datad)
    datatmpd='{}/04_specificity/tmp'.format(datad)
    dataoutd='{}/04_specificity/out'.format(datad)
    for dp in [dataind,datatmpd,dataoutd]: 
        if not exists(dp) or force:
            makedirs(dp,exist_ok=force)

    dguides=pd.read_csv(dguidesp,sep='\t')
    dguides.to_csv('{}/{}'.format(dataind,basename(dguidesp)),sep='\t')
    dguides=dguides.set_index('guideId')
    with open('{}/batchId.fa'.format(datatmpd),'w') as f:
        for gi in dguides.index:
            f.write('>{}\n{}\n'.format(gi,dguides.loc[gi,'targetSeq']))

    genomeDir='pub/release-92/fasta/'
    org='saccharomyces_cerevisiae'
    genomefn='dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.chromosome.I.fa'
    genomegffp='pub/release-92/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.92.gff3.gz'

    # def dguide2fltofftarget(cfg):
    # get dna and protein sequences 
    #--
    # junk
    # write debug output to stdout
    DEBUG = False
    #DEBUG = True

    # BWA: allow up to X mismatches
    maxMMs=4

    # maximum number of occurences in the genome to get flagged as repeats. 
    # This is used in bwa samse, when converting the same file
    # and for warnings in the table output.
    MAXOCC = 60000

    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000/MAXOCC

    # the length of the guide sequence, set by setupPamInfo
    GUIDELEN=None
    # length of the PAM sequence
    PAMLEN=None

    # are we doing a Cpf1 run?
    # this variable changes almost all processing and
    # has to be set on program start, as soon as we know 
    # the PAM we're running on
    cpf1Mode=None

    # global flag to indicate if we're run from command line or as a CGI

    # names/order of efficiency scores to show in UI
    scoreNames = ["fusi", "crisprScan"]
    # allScoreNames = ["fusi", "fusiOld", "chariRank", "ssc", "doench", "wang", "crisprScan", "aziInVitro"]

    cpf1ScoreNames = ["seqDeepCpf1"]

    saCas9ScoreNames = ["najm"]

    #
    DEBUG=True
    pam='NGG'
    genomep='{}/{}/{}'.format(genomeDir,org,genomefn)
    # genomep='crisporWebsite/genomes/sacCer3/sacCer3.fa'
    class options(object):
        def __init__(self, pam):
            self.skipAlign = False
            self.debug=DEBUG
    GUIDELEN,cpf1Mode,addGenePlasmids,PAMLEN,scoreNames=setupPamInfo(pam,setupPamInfo,scoreNames)
    skipAlign = False
    if options.skipAlign:
        skipAlign = True

    # get sequence
    # seqs = parseFasta(open(inSeqFname))
    batchId='batchId'
    batchBase = join(datatmpd, batchId)
    otBedFname = batchBase+".bed"
    print(otBedFname)
    faFname = batchBase+".fa"

    matchesBedFname = batchBase+".matches.bed"
    saFname = batchBase+".sa"
    samp = batchBase+".sam"
    genomeDir = dirname(genomep) # make var local, see below

    open(matchesBedFname, "w") # truncate to 0 size

    # increase MAXOCC if there is only a single query, but only in CGI mode
    maxDiff = maxMMs
    seqLen = GUIDELEN

    bwaM = MFAC*MAXOCC # -m is queue size in bwa
    cmd = "$BIN/bwa aln -o 0 -m %(bwaM)s -n %(maxDiff)d -k %(maxDiff)d -N -l %(seqLen)d %(genomep)s %(faFname)s > %(saFname)s" % locals()
    runCmd(cmd)

    cmd = "$BIN/bwa samse -n %(MAXOCC)d %(genomep)s %(saFname)s %(faFname)s > %(samp)s" % locals()
    runCmd(cmd)
    #----make tables-----------
    gff_colns = ['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    bed_colns = ['chromosome','start','end','id','NM','strand']

    samfile=pysam.AlignmentFile(samp, "rb")
    # with pysam.AlignmentFile(samp, "rb") as samfile:
    #     print(samfile.contigs)
    # cols=['chromosome','start','end','id','NM','strand']
    dalignbed=pd.DataFrame(columns=bed_colns)
    # rowi=1
    for read in samfile.fetch():
    #     print(read) 
    #     print(read.reference_name)
    #     algnid='{}|{}{}|{}|{}'.format(dalignbed.loc[rowi,'chromosome'],
    #              dalignbed.loc[rowi,'strand'],dalignbed.loc[rowi,'start'],read.cigarstring,dalignbed.loc[rowi,'NM'])
        algnids=[]
        algnids.append('{}|{}{}|{}|{}'.format(read.reference_name,
                 '-' if read.is_reverse else '+',read.positions[0],read.cigarstring,read.get_tag('NM')))
        algnids+=['|'.join(s.split(',')) for s in read.get_tag('XA').split(';') if len(s.split(','))>1]
        chroms=[]
        starts=[]
        ends=[]
        ids=algnids
        NMs=[]
        strands=[]    
        for a in algnids:
            chroms.append(a.split('|')[0])
            starts.append(int(a.split('|')[1][1:]))
            ends.append(int(a.split('|')[1][1:])+str2num(a.split('|')[2]))
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
        dalignbed_['grnaId']=read.qname
        dalignbed = dalignbed.append(dalignbed_,ignore_index=True)
    #     break
    samfile.close()

    alignmentbedp='{}/alignment.bed'.format(datatmpd)
    dalignbed.loc[:,bed_colns].to_csv(alignmentbedp,sep='\t',
                    header=False,index=False)
    # loc[:,cols].to
    # bedtools sort -i data_test_human/04_specificity/tmp/alignment.bed > data_test_human/04_specificity/tmp/alignment.sorted.bed
    # bedtools sort -i pub/release-92/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.92.gff3.gz > pub/release-92/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.92.sorted.gff3.gz
    # bedtools intersect -wa -wb -loj -a data_test_human/04_specificity/tmp/alignment.sorted.bed -b pub/release-92/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.92.sorted.gff3.gz > data_test_human/04_specificity/tmp/annotations.bed
    alignmentbedsortedp=alignmentbedp+'.sorted.bed'
    cmd='bedtools sort -i {} > {}'.format(alignmentbedp,alignmentbedsortedp)
    runCmd()
    genomegffsortedp=genomegffp+'.sorted.gff3.gz'
    cmd='bedtools sort -i {} > {}'.format(genomegffp,genomegffsortedp)
    runCmd()
    annotationsbedp='{}/annotations.bed'.format(datatmpd)
    cmd='bedtools intersect -wa -wb -loj -a {} -b {} > {}'.format(alignmentbedsortedp,genomegffsortedp,annotationsbedp)
    runCmd()

    # bedtools annotate -i pub/release-92/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.92.gff3.gz -files data_test_human/04_specificity/tmp/alignment.bed
    # [OPTIONS] -i <BED/GFF/VCF> -files FILE1 FILE2 FILE3 ... FILEn

    dannots=pd.read_csv('data_test_human/04_specificity/tmp/annotations.bed',sep='\t',
               names=bed_colns+gff_colns)

    dalignbed=dalignbed.set_index('grnaId').join(dguides)

    dalignbed=dalignbed.reset_index().set_index('id')

    # dalignid2seq=pd.DataFrame(columns=['sequence'])
    # dalignid2seq.index.name='id'
    alignedfastap='{}/alignment.fa'.format(datatmpd)
    for seq in SeqIO.parse(alignedfastap,"fasta"):
        dalignbed.loc[seq.id,'aligned sequence']=str(seq.seq)
    #     break

    dalignbed.loc[:,'Hamming distance']=[hamming_distance(dalignbed.loc[i,'targetSeq'], dalignbed.loc[i,'aligned sequence']) for i in dalignbed.index]

    for i in dalignbed.index:
        dalignbed.loc[i,'alignment'],dalignbed.loc[i,'alignment: score']=align(dalignbed.loc[i,'targetSeq'],dalignbed.loc[i,'aligned sequence'])

    dcombo=dalignbed.join(dannots.set_index('id'),rsuffix='.2').head()

    dcombo.to_csv('{}/dcombo.tsv'.fomar(dataoutd),sep='\t')
if __name__ == '__main__':
    main()
