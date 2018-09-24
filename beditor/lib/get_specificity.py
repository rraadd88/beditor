#!usr/bin/python

import logging
import subprocess
import re
import sys
import logging 

from collections import defaultdict
from os.path import join, basename, dirname, abspath, exists
from os import makedirs,stat

import pandas as pd
# import modin.pandas as pd

import pysam
import numpy as np
from glob import glob

from beditor.lib.global_vars import bed_colns,gff_colns    
from beditor.lib.io_sys import runbashcmd
from beditor.lib.io_seqs import fa2df 
from beditor.lib.io_dfs import set_index,del_Unnamed,df2info 

from beditor.lib.io_dfs import lambda2cols
from beditor.lib.io_nums import str2num
from beditor.lib.io_seqs import gffatributes2ids,hamming_distance,align

def dguides2guidessam(cfg,dguides):    
    """
    Aligns guides to genome and gets SAM file
    step#1

    :param cfg: configuration dict
    :param dguides: dataframe of guides
    """
    datatmpd=cfg['datatmpd']
    dguides=set_index(dguides,'guide: id')
    guidels=dguides.loc[:,'guide+PAM length'].unique()
    for guidel in guidels:
        logging.debug(f"now aligning guides of length {guidel}")
        guidesfap = f'{datatmpd}/01_guides_guidel{guidel:02}.fa'
        logging.info(basename(guidesfap))
        if not exists(guidesfap) or cfg['force']:
            with open(guidesfap,'w') as f:
                for gi in dguides.index:
                    f.write('>{}\n{}\n'.format(gi.replace(' ','_'),dguides.loc[gi,'guide+PAM sequence']))
        ## BWA alignment command is adapted from cripror 
        ## https://github.com/rraadd88/crisporWebsite/blob/master/crispor.py
        # BWA: allow up to X mismatches
        # maximum number of occurences in the genome to get flagged as repeats. 
        # This is used in bwa samse, when converting the sam file
        # and for warnings in the table output.
        MAXOCC = 60000

        # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
        MFAC = 2000000/MAXOCC

        genomep=cfg['genomep']
        genomed = dirname(genomep) # make var local, see below
        genomegffp=cfg['genomegffp']

        # increase MAXOCC if there is only a single query, but only in CGI mode
        bwaM = MFAC*MAXOCC # -m is queue size in bwa
        guidessap = f'{datatmpd}/01_guides_guidel{guidel:02}.sa'
        logging.info(basename(guidessap))
        if not exists(guidessap) or cfg['force']:
            cmd=f"{cfg['bwa']} aln -t 1 -o 0 -m {bwaM} -n {cfg['mismatches_max']} -k {cfg['mismatches_max']} -N -l {guidel} {genomep} {guidesfap} > {guidessap} 2> {guidessap}.log"
            runbashcmd(cmd)

        guidessamp = f'{datatmpd}/01_guides_guidel{guidel:02}.sam'
        logging.info(basename(guidessamp))        
        if not exists(guidessamp) or cfg['force']:
            cmd=f"{cfg['bwa']} samse -n {MAXOCC} {genomep} {guidessap} {guidesfap} > {guidessamp} 2> {guidessamp}.log"
            runbashcmd(cmd)
    return cfg

def guidessam2dalignbed(cfg):
    """
    Processes SAM file to get the genomic coordinates in BED format
    step#2

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    alignmentbedp=cfg['alignmentbedp']
    dalignbedp=cfg['dalignbedp']
    logging.info(basename(dalignbedp))
    if not exists(alignmentbedp) or cfg['force']:
        #input/s
        guidessamps=glob(f'{datatmpd}/01_guides_guidel*.sam')
        for guidessamp in guidessamps:
            if stat(guidessamp).st_size != 0:
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
                    dalignbed = dalignbed.append(dalignbed_,ignore_index=True,sort=True)
                #     break
                samfile.close()
            else:
                logging.warning(f"file is empty: {guidessamp}")
        dalignbed.to_csv(dalignbedp,sep='\t')
        from beditor.lib.io_nums import str2numorstr
        dalignbed['chromosome']=dalignbed.apply(lambda x : str2numorstr(x['chromosome']),axis=1)
        dalignbed=dalignbed.sort_values(['chromosome','start','end'], ascending=[True, True, True])
        dalignbed.loc[:,bed_colns].to_csv(alignmentbedp,sep='\t',
                        header=False,index=False,
                        chunksize=5000)
    return cfg

def dalignbed2annotationsbed(cfg):
    """
    Get annotations from the aligned BED file
    step#3

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    alignmentbedp=cfg['alignmentbedp']    
    alignmentbedsortedp=alignmentbedp+'.sorted.bed'
    logging.info(basename(alignmentbedsortedp))
    if not exists(alignmentbedsortedp) or cfg['force']:
        cmd='{} sort -i {} > {}'.format(cfg['bedtools'],alignmentbedp,alignmentbedsortedp)
        runbashcmd(cmd)
    
    genomegffsortedp=cfg['genomegffp']+'.sorted.gff3.gz'
    logging.info(basename(genomegffsortedp))
    if not exists(genomegffsortedp):    
        cmd=f"{cfg['bedtools']} sort -i {cfg['genomegffp']} > {genomegffsortedp}"
        runbashcmd(cmd)

    annotationsbedp='{}/03_annotations.bed'.format(datatmpd)
    cfg['annotationsbedp']=annotationsbedp
    logging.info(basename(annotationsbedp))
    if not exists(annotationsbedp) or cfg['force']:    
        cmd=f"{cfg['bedtools']} intersect -wa -wb -loj -a {alignmentbedsortedp} -b {genomegffsortedp} > {annotationsbedp}"
        runbashcmd(cmd)
    return cfg

def dalignbed2dalignbedguides(cfg):
    """
    Get guide seqeunces from the BED file
    step#4

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    dalignbed=del_Unnamed(pd.read_csv(cfg['dalignbedp'],sep='\t'))
    dguides=set_index(del_Unnamed(pd.read_csv(cfg['dguidesp'],sep='\t')),'guide: id')
    
#     if the error in human, use: `cut -f 1 data/alignment.bed.sorted.bed | sort| uniq -c | grep -v CHR | grep -v GL | grep -v KI`
    dalignbedguidesp=cfg['dalignbedguidesp']
    logging.info(basename(dalignbedguidesp))
    if not exists(dalignbedguidesp) or cfg['force']:
        dalignbed=pd.merge(dalignbed,dguides,on='guide: id',suffixes=('', '.1'))
        dalignbed.to_csv(dalignbedguidesp,'\t')
    return cfg

def alignmentbed2dalignedfasta(cfg):
    """
    Get sequences in FASTA format from BED file
    step#5

    :param cfg: configuration dict
    """    
    datatmpd=cfg['datatmpd']
    alignmentbedp=cfg['alignmentbedp']    
    dalignedfastap=cfg['dalignedfastap']
    logging.info(basename(dalignedfastap))
    if not exists(dalignedfastap) or cfg['force']:
        alignedfastap='{}/05_alignment.fa'.format(datatmpd)
        if not exists(alignedfastap) or cfg['force']:
            cmd=f"{cfg['bedtools']} getfasta -s -name -fi {cfg['genomep']} -bed {alignmentbedp} -fo {alignedfastap}"
            runbashcmd(cmd)

        dalignedfasta=fa2df(alignedfastap)
        dalignedfasta.columns=['aligned sequence']
        dalignedfasta=dalignedfasta.loc[(dalignedfasta.apply(lambda x: not 'N' in x['aligned sequence'],axis=1)),:] #FIXME bwa aligns to NNNNNs
        dalignedfasta.index=[i.split('(')[0] for i in dalignedfasta.index] # for bedtools 2.27, the fasta header now has hanging (+) or (-)
        dalignedfasta.index.name='id'
        dalignedfasta.to_csv(dalignedfastap,sep='\t')
    return cfg

def dalignbed2dalignbedguidesseq(cfg):
    """
    Get sequences from BED file
    step#6

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    dalignbedguides=del_Unnamed(pd.read_csv(cfg['dalignbedguidesp'],sep='\t'))
    dalignedfasta=del_Unnamed(pd.read_csv(cfg['dalignedfastap'],sep='\t'))
    dalignbedguidesseqp=cfg['dalignbedguidesseqp']
    logging.info(basename(dalignbedguidesseqp))
    if not exists(dalignbedguidesseqp) or cfg['force']:        
        dalignbedguidesseq=pd.merge(dalignbedguides,dalignedfasta,on='id',suffixes=('', '.2'))
        dalignbedguidesseq=dalignbedguidesseq.dropna(subset=['aligned sequence'],axis=0)

        # dalignbed.index.name='id'
        dalignbedguidesseq=dalignbedguidesseq.drop_duplicates()
        dalignbedguidesseq.to_csv(dalignbedguidesseqp,sep='\t')
    return cfg

def dalignbedguidesseq2dalignbedstats(cfg):
    """
    Gets scores for guides
    step#7

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    dalignbedguidesseq=del_Unnamed(pd.read_csv(cfg['dalignbedguidesseqp'],sep='\t'))
    
    dalignbedstatsp=cfg['dalignbedstatsp']  
    logging.info(basename(dalignbedstatsp))
    if not exists(dalignbedstatsp) or cfg['force']:
        df=dalignbedguidesseq.apply(lambda x: align(x['guide+PAM sequence'],x['aligned sequence']),
                           axis=1).apply(pd.Series)
        df.columns=['alignment','alignment: score']
        dalignbedstats=dalignbedguidesseq.join(df)
        del df
        dalignbedstats.to_csv(dalignbedstatsp,sep='\t')
    return cfg
def dannots2dalignbed2dannotsagg(cfg):
    """
    Aggregate annotations per guide
    step#8

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    
    daannotp=f'{datatmpd}/08_dannot.tsv'
    cfg['daannotp']=daannotp
    dannotsaggp=cfg['dannotsaggp']
    logging.info(basename(daannotp))
    if ((not exists(daannotp)) and (not exists(dannotsaggp))) or cfg['force']:
        dannots=pd.read_csv(cfg['annotationsbedp'],sep='\t',
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
        logging.debug('or this step takes more time?')
        dannots.to_csv(daannotp,sep='\t')
    else:
        dannots=pd.read_csv(daannotp,sep='\t',low_memory=False)
        dannots=del_Unnamed(dannots)
    logging.info(basename(dannotsaggp))
    if not exists(dannotsaggp) or cfg['force']:
        if not 'dannots' in locals():
            dannots=pd.read_table(daannotp,low_memory=False)
        dannots=del_Unnamed(dannots)
        dannots=dannots.reset_index()
        
        dannotsagg=pd.DataFrame(dannots.groupby('id')['annotations count'].agg('sum'))-1
        dannotsagg.loc[dannotsagg['annotations count']==0,'region']='intergenic'
        dannotsagg.loc[dannotsagg['annotations count']!=0,'region']='genic'

        alignids=dannots['id'].unique()#[:15]
        logging.debug('start of the slowest step')
        for alignidi in range(len(alignids)):
            alignid=alignids[alignidi]
            dannoti=dannots.loc[dannots['id']==alignid,:]
            if len(dannoti.shape)==1:
                dannoti=pd.DataFrame(dannoti).T
            dannoti=dannoti.loc[dannoti['type']!='chromosome',:].drop_duplicates(subset=['start annotation','end annotation'])
            for col in ['type','gene name','gene id','transcript id','protein id','exon id']:    
                dannotsagg.loc[alignid,col+'s']=";".join(np.unique(dannoti[col].fillna('nan').tolist()))
        logging.debug('end of the slowest step')
            
        del dannots    
        dannotsagg=dannotsagg.reset_index()
        dannotsagg.to_csv(dannotsaggp,sep='\t')
    return cfg

def dannotsagg2dannots2dalignbedannot(cfg):
    """
    Map aggregated annotations to guides
    step#9

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']
    
    dannotsagg=del_Unnamed(pd.read_csv(cfg['dannotsaggp'],sep='\t'))
    dalignbedstats=del_Unnamed(pd.read_csv(cfg['dalignbedstatsp'],sep='\t'))
    dalignbedannotp=cfg['dalignbedannotp']
    logging.info(basename(dalignbedannotp))
    if not exists(dalignbedannotp) or cfg['force']:
        # df2info(dalignbed)
        # df2info(dannotsagg)
        dalignbedannot=dalignbedstats.set_index('id').join(set_index(dannotsagg,'id'),
                                              rsuffix=' annotation')
        dalignbedannot['NM']=dalignbedannot['NM'].apply(int)
        from beditor.lib.get_scores import get_beditorscore_per_alignment,get_cfdscore
        dalignbedannot['beditor score']=dalignbedannot.apply(lambda x : get_beditorscore_per_alignment(NM=x['NM'],
                               genic=True if x['region']=='genic' else False,
                               alignment=x['alignment'],
                               pam_length=len(x['PAM']),
                               pam_position=x['original position'],
                               # test=cfg['test'],
                                ),axis=1) 
        dalignbedannot['CFD score']=dalignbedannot.apply(lambda x : get_cfdscore(x['guide+PAM sequence'].upper(), x['aligned sequence'].upper()), axis=1)            
        dalignbedannot.to_csv(dalignbedannotp,sep='\t')
    return cfg

def dalignbedannot2daggbyguide(cfg):
    """
    Aggregate annotations per alignment to annotations per guide.
    step#10

    :param cfg: configuration dict
    """
    datatmpd=cfg['datatmpd']

    dalignbedannot=del_Unnamed(pd.read_csv(cfg['dalignbedannotp'],sep='\t',low_memory=False))
    
    daggbyguidep='{}/10_daggbyguide.tsv'.format(datatmpd)      
    logging.info(basename(daggbyguidep))
    if not exists(daggbyguidep) or cfg['force']:
        daggbyguide=dalignbedannot.loc[(dalignbedannot['NM']==0),['guide: id','guide+PAM sequence','gene names', 'gene ids','transcript ids']].drop_duplicates(subset=['guide: id'])
        if len(daggbyguide)!=0:
            daggbyguide=set_index(daggbyguide,'guide: id')            
            guideids=daggbyguide.index.tolist()
            for gi in range(len(guideids)):
                gid=guideids[gi]
                dalignbedannoti=dalignbedannot.loc[dalignbedannot['guide: id']==gid,:]
                if len(dalignbedannoti.shape)==1:
                    dalignbedannoti=pd.DataFrame(dalignbedannoti).T
                for col in ['types','gene names','gene ids','transcript ids','protein ids','exon ids']:    
                    daggbyguide.loc[gid,col]=";".join(np.unique(dalignbedannoti[col].fillna('nan').tolist()))
            from beditor.lib.get_scores import get_beditorscore_per_guide
            for guideid in daggbyguide.index:
                dalignbedannotguide=dalignbedannot.loc[(dalignbedannot['guide: id']==guideid),:]
                daggbyguide.loc[guideid,'beditor score']=get_beditorscore_per_guide(guide_seq=dalignbedannotguide['guide+PAM sequence'].unique()[0], 
                                           strategy=dalignbedannotguide['strategy'].unique()[0],
                                           align_seqs_scores=dalignbedannotguide['beditor score'],
                                           BEs=cfg['BEs']
    #                                        test=cfg['test']
                                          )
                daggbyguide.loc[guideid,'CFD score']=dalignbedannotguide['CFD score'].mean() #FIXME if mean is not appropriate
            daggbyguide['beditor score (log10)']=daggbyguide['beditor score'].apply(np.log10)
            dalignbedannot['alternate alignments count']=1
            daggbyguide=daggbyguide.join(pd.DataFrame(dalignbedannot.groupby('guide: id')['alternate alignments count'].agg('sum')))
            daggbyguide.to_csv(daggbyguidep,sep='\t')
            daggbyguide.to_csv(cfg['dofftargetsp'],sep='\t')
    return cfg

def dguides2offtargets(cfg):
    """
    All the processes in offtarget detection are here.
    
    :param cfg: Configuration settings provided in .yml file
    """
    from beditor.lib.global_vars import saveemptytable
    
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    dofftargetsp='{}/dofftargets.tsv'.format(cfg['datad'])  
    
    stepn='04_offtargets'
    logging.info(stepn)
    dguidesp=f"{cfg[cfg['step']-1]}/dguides.tsv"
    if not exists(dguidesp):        
        logging.warning(f"not found {dguidesp}")
        return saveemptytable(cfg,dofftargetsp)
    dguides=pd.read_csv(dguidesp,sep='\t')
    if len(dguides)==0:
        logging.warning(f"dguides is empty.")
        return saveemptytable(cfg,dofftargetsp)       

    cfg['datatmpd']=f"{cfg['datad']}/tmp"
    for dp in [cfg['datatmpd']]:
        if not exists(dp):
            makedirs(dp)
            
    step2doutp={
    1:'01_guides_guidel*.fa',
    2:'02_dalignbed.tsv',
    3:'03_annotations.bed',
    4:'04_dalignbedguides.tsv',
    5:'05_dalignedfasta.tsv',
    6:'06_dalignbedguidesseq.tsv',
    7:'07_dalignbedstats.tsv',
    8:'08_dannotsagg.tsv',
    9:'09_dalignbedannot.tsv',
    10:'10_daggbyguide.tsv',
    }
    cfg['dguidesp']=dguidesp
    cfg['alignmentbedp']=f"{cfg['datatmpd']}/02_alignment.bed"
    cfg['dalignbedp']=f"{cfg['datatmpd']}/02_dalignbed.tsv"
    cfg['dalignbedguidesp']=f"{cfg['datatmpd']}/04_dalignbedguides.tsv"
    cfg['dalignedfastap']=f"{cfg['datatmpd']}/05_dalignedfasta.tsv"
    cfg['dalignbedguidesseqp']=f"{cfg['datatmpd']}/06_dalignbedguidesseq.tsv"
    cfg['dalignbedstatsp']=f"{cfg['datatmpd']}/07_dalignbedstats.tsv"
    cfg['dannotsaggp']=f"{cfg['datatmpd']}/08_dannotsagg.tsv"
    cfg['dalignbedannotp']=f"{cfg['datatmpd']}/09_dalignbedannot.tsv"
    cfg['daggbyguidep']=f"{cfg['datatmpd']}/10_daggbyguide.tsv"

    #check which step to process
    for step in range(2,10+1,1):
        if not exists(f"{cfg['datatmpd']}/{step2doutp[step]}"):
            if step==2:
                step='all'
            break
    logging.info(f'process from step:{step}')
    cfg['dofftargetsp']='{}/dofftargets.tsv'.format(cfg['datad'])
    if not exists(cfg['dofftargetsp']) or cfg['force']:
        if step==1 or step=='all':
            cfg=dguides2guidessam(cfg,dguides)
        if step==2 or step=='all' or (not cfg is None):
            cfg=guidessam2dalignbed(cfg)
        if step==3 or step=='all' or (not cfg is None):
            cfg=dalignbed2annotationsbed(cfg)
        if step==4 or step=='all' or (not cfg is None):
            cfg=dalignbed2dalignbedguides(cfg)
        if step==5 or step=='all' or (not cfg is None):
            cfg=alignmentbed2dalignedfasta(cfg)
        if step==6 or step=='all' or (not cfg is None):
            cfg=dalignbed2dalignbedguidesseq(cfg)
        if step==7 or step=='all' or (not cfg is None):
            cfg=dalignbedguidesseq2dalignbedstats(cfg)
        if step==8 or step=='all' or (not cfg is None):
            cfg=dannots2dalignbed2dannotsagg(cfg)
        if step==9 or step=='all' or (not cfg is None):
            cfg=dannotsagg2dannots2dalignbedannot(cfg)
        if step==10 or step=='all' or (not cfg is None):
            cfg=dalignbedannot2daggbyguide(cfg)

        if cfg is None:        
            logging.warning(f"no alignment found")
            cfg['step']=4
            return saveemptytable(cfg,cfg['dofftargetsp'])
        import gc
        gc.collect()