#!usr/bin/python

import logging
import subprocess
import re
import sys
import logging 

from collections import defaultdict
from os.path import join, basename, dirname, abspath, exists
from os import makedirs

import pandas as pd
# import modin.pandas as pd

import pysam
import numpy as np
from glob import glob

from beditor.lib.io_sys import runbashcmd
from beditor.lib.io_seqs import fa2df 
from beditor.lib.io_dfs import set_index,del_Unnamed,df2info 

from beditor.lib.io_dfs import lambda2cols
from beditor.lib.io_nums import str2num
from beditor.lib.io_seqs import gffatributes2ids,hamming_distance,align

def dguides2offtargets(cfg):
    """
    All the processes in offtarget detection are here.
    :param cfg: Configuration settings provided in .yml file
    """
    cfg['datad']=cfg[cfg['step']]
    cfg['plotd']=cfg['datad']
    
    stepn='04_offtargets'
    logging.info(stepn)
    dguidesp=f"{cfg[cfg['step']-1]}/dguides.tsv"
    datatmpd=f"{cfg['datad']}/tmp"
    for dp in [datatmpd]:
        if not exists(dp):
            makedirs(dp)

    dofftargetsp='{}/dofftargets.tsv'.format(cfg['datad'])  
    if not exists(dofftargetsp) or cfg['force']:
        step2dscr={1: 'align_guides',
                   2: 'align_guides',
                   3: 'align_guides',
                  }
        
        dguides=pd.read_csv(dguidesp,sep='\t')
    #     dguides.to_csv('{}/{}'.format(cfg['datad'],basename(dguidesp)),sep='\t')
        dguides=dguides.set_index('guide: id')
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
    #             cmd="{} aln -t {} -o 0 -m {} -n {} -k {} -N -l {} {} {} > {}".format(cfg['bwa'],1,bwaM,cfg['mismatches_max'],cfg['mismatches_max'],cfg['guidel'],genomep,guidesfap,guidessap)
    #             runbashcmd(cmd)
                cmd=f"{cfg['bwa']} aln -t 1 -o 0 -m {bwaM} -n {cfg['mismatches_max']} -k {cfg['mismatches_max']} -N -l {guidel} {genomep} {guidesfap} > {guidessap} 2> {guidessap}.log"
                runbashcmd(cmd)

            guidessamp = f'{datatmpd}/01_guides_guidel{guidel:02}.sam'
            logging.info(basename(guidessamp))        
            if not exists(guidessamp) or cfg['force']:
                cmd=f"{cfg['bwa']} samse -n {MAXOCC} {genomep} {guidessap} {guidesfap} > {guidessamp} 2> {guidessamp}.log"
                runbashcmd(cmd)

        #----make tables-----------
        from beditor.lib.global_vars import bed_colns,gff_colns    

        #output
        alignmentbedp=f'{datatmpd}/02_alignment.bed'
        dalignbedp=f'{datatmpd}/02_dalignbed.tsv'
        logging.info(basename(dalignbedp))
        if not exists(alignmentbedp) or cfg['force']:
            #input/s
            guidessamps=glob(f'{datatmpd}/01_guides_guidel*.sam')
            for guidessamp in guidessamps:
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
            cmd='{} sort -i {} > {}'.format(cfg['bedtools'],alignmentbedp,alignmentbedsortedp)
            runbashcmd(cmd)

        genomegffsortedp=genomegffp+'.sorted.gff3.gz'
        logging.info(basename(genomegffsortedp))
        if not exists(genomegffsortedp):    
            cmd='{} sort -i {} > {}'.format(cfg['bedtools'],genomegffp,genomegffsortedp)
            runbashcmd(cmd)

        annotationsbedp='{}/03_annotations.bed'.format(datatmpd)
        logging.info(basename(annotationsbedp))
        if not exists(annotationsbedp) or cfg['force']:    
            cmd='{} intersect -wa -wb -loj -a {} -b {} > {}'.format(cfg['bedtools'],alignmentbedsortedp,genomegffsortedp,annotationsbedp)
            runbashcmd(cmd)
    #     if the error in human, use: `cut -f 1 data/alignment.bed.sorted.bed | sort| uniq -c | grep -v CHR | grep -v GL | grep -v KI`
        dalignbedguidesp='{}/04_dalignbedguides.tsv'.format(datatmpd)
        logging.info(basename(dalignbedguidesp))
        if not exists(dalignbedguidesp) or cfg['force']:
            dalignbed=set_index(dalignbed,'guide: id').join(dguides,rsuffix='.1')
            dalignbed=set_index(dalignbed,'guide: id')
            dalignbed.to_csv(dalignbedguidesp,'\t')
        else:
            dalignbed=pd.read_csv(dalignbedguidesp,'\t')
            dalignbed=del_Unnamed(dalignbed)

        # dalignid2seq=pd.DataFrame(columns=['sequence'])
        # dalignid2seq.index.name='id'
        dalignedfastap='{}/05_dalignedfasta.tsv'.format(datatmpd)
        logging.info(basename(dalignedfastap))
        if not exists(dalignedfastap) or cfg['force']:
            alignedfastap='{}/05_alignment.fa'.format(datatmpd)
            if not exists(alignedfastap) or cfg['force']:
                cmd=f"{cfg['bedtools']} getfasta -s -name -fi {genomep} -bed {alignmentbedp} -fo {alignedfastap}"
                runbashcmd(cmd)

            dalignedfasta=fa2df(alignedfastap)
            dalignedfasta.columns=['aligned sequence']
            dalignedfasta=dalignedfasta.loc[(dalignedfasta.apply(lambda x: not 'N' in x['aligned sequence'],axis=1)),:] #FIXME bwa aligns to NNNNNs
    #         dalignedfasta=dalignedfasta.loc[[False if np.unique(list(s.upper()))=='N' else True for s in dalignedfasta['aligned sequence']],:]
            dalignedfasta.index=[i.split('(')[0] for i in dalignedfasta.index] # for bedtools 2.27, the fasta header now has hanging (+) or (-)
            dalignedfasta.index.name='id'
            dalignedfasta.to_csv(dalignedfastap,sep='\t')
        else:
            dalignedfasta=pd.read_csv(dalignedfastap,sep='\t')        
            dalignedfasta=dalignedfasta.drop([c for c in dalignbed if 'Unnamed' in c],axis=1)
            dalignedfasta=del_Unnamed(dalignedfasta)

        dalignbedguidesseqp='{}/06_dalignbedguidesseq.tsv'.format(datatmpd)
        logging.info(basename(dalignbedguidesseqp))
        if not exists(dalignbedguidesseqp) or cfg['force']:        
    #         from beditor.lib.io_dfs import df2info
    #         if cfg['test']:
    #             df2info(dalignbed)
    #             df2info(dalignedfasta)

    #         dalignedfasta['guide: id']=dalignedfasta.apply(lambda x : x['guide: id'].replace('_',' '),axis=1)
            dalignbed=set_index(dalignbed,'id').join(set_index(dalignedfasta,'id'))
            dalignbed=dalignbed.dropna(subset=['aligned sequence'],axis=0)

            dalignbed.index.name='id'
            dalignbed=dalignbed.drop_duplicates()
            dalignbed.to_csv(dalignbedguidesseqp,sep='\t')
        else:
            dalignbed=pd.read_csv(dalignbedguidesseqp,sep='\t',low_memory=False)
            dalignbed=del_Unnamed(dalignbed)

        dalignbedstatsp='{}/07_dalignbedstats.tsv'.format(datatmpd)  
        logging.info(basename(dalignbedstatsp))
        if not exists(dalignbedstatsp) or cfg['force']:
    #         dalignbed['Hamming distance']=dalignbed.apply(lambda x : hamming_distance(x['guide+PAM sequence'], x['aligned sequence']),axis=1)
            df=dalignbed.apply(lambda x: align(x['guide+PAM sequence'],x['aligned sequence']),
                               axis=1).apply(pd.Series)
            df.columns=['alignment','alignment: score']
            dalignbed=dalignbed.join(df)
            del df
            dalignbed.to_csv(dalignbedstatsp,sep='\t')
        else:
            dalignbed=pd.read_csv(dalignbedstatsp,sep='\t',low_memory=False)
            dalignbed=del_Unnamed(dalignbed)

        daannotp='{}/08_dannot.tsv'.format(datatmpd)  
        dannotsaggp='{}/08_dannotsagg.tsv'.format(datatmpd)  
        logging.info(basename(daannotp))
        if (not exists(daannotp)) or (not exists(dannotsaggp)) or cfg['force']:
            dannots=pd.read_csv(annotationsbedp,sep='\t',
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

        logging.info(basename(dannotsaggp))
        if not exists(dannotsaggp) or cfg['force']:
            dannotsagg=pd.DataFrame(dannots.groupby('id')['annotations count'].agg('sum'))-1
            dannotsagg.loc[dannotsagg['annotations count']==0,'region']='intergenic'
            dannotsagg.loc[dannotsagg['annotations count']!=0,'region']='genic'

            dannots=pd.read_csv(daannotp,sep='\t',low_memory=False)
            dannots=del_Unnamed(dannots)
            dannots=dannots.reset_index()
            alignids=dannots['id'].unique()#[:15]

            logging.debug('start of the slowest step')
            for alignidi in range(len(alignids)):
                alignid=alignids[alignidi]
            #     if dannots.loc[i,'type'].value_counts().sum()==1:
                dannoti=dannots.loc[dannots['id']==alignid,:]
                if len(dannoti.shape)==1:
                    dannoti=pd.DataFrame(dannoti).T
                dannoti=dannoti.loc[dannoti['type']!='chromosome',:].drop_duplicates(subset=['start annotation','end annotation'])
    #             if cfg['test']:
    # #                 print(alignid)
    #                 print(dannoti.shape)
                for col in ['type','gene name','gene id','transcript id','protein id','exon id']:    
            #         print(";".join(dannoti[col].fillna('nan').tolist()))
                    dannotsagg.loc[alignid,col+'s']=";".join(np.unique(dannoti[col].fillna('nan').tolist()))
            logging.debug('end of the slowest step')
                
            del dannots    
            dannotsagg.to_csv(dannotsaggp,sep='\t')
        else:
            dannotsagg=pd.read_csv(dannotsaggp,sep='\t',low_memory=False)
            dannotsagg=del_Unnamed(dannotsagg)

        dalignbedannotp='{}/09_dalignbedannot.tsv'.format(datatmpd)  
        logging.info(basename(dalignbedannotp))
        if not exists(dalignbedannotp) or cfg['force']:
            dalignbedannot=set_index(dalignbed,'id').join(set_index(dannotsagg,'id'),
                                                  rsuffix=' annotation')
            dalignbedannot['NM']=dalignbedannot['NM'].apply(int)
            
            from beditor.lib.get_scores import get_beditorscore_per_alignment,get_cfdscore
            dalignbedannot['beditor score']=dalignbedannot.apply(lambda x : get_beditorscore_per_alignment(x['NM'],cfg['mismatches_max'],
                            True if x['region']=='genic' else False,
                            x['alignment'],
                            #                                                                                                        test=cfg['test'],
                            ),axis=1) 
            dalignbedannot['CFD score']=dalignbedannot.apply(lambda x : get_cfdscore(x['guide+PAM sequence'].upper(), x['aligned sequence'].upper()), axis=1)            
            dalignbedannot.to_csv(dalignbedannotp,sep='\t')
        else:
            dalignbedannot=pd.read_csv(dalignbedannotp,sep='\t',low_memory=False)
            dalignbedannot=del_Unnamed(dalignbedannot)

        daggbyguidep='{}/10_daggbyguide.tsv'.format(datatmpd)      
        logging.info(basename(daggbyguidep))
        if not exists(daggbyguidep) or cfg['force']:
            daggbyguide=dalignbedannot.loc[(dalignbedannot['NM']==0),['guide: id','guide+PAM sequence','gene names', 'gene ids','transcript ids']].drop_duplicates(subset=['guide: id'])
#             if cfg['test']:
#                 df2info(daggbyguide)
            if len(daggbyguide)!=0:
                daggbyguide=set_index(daggbyguide,'guide: id')            
        #---
                guideids=daggbyguide.index.tolist()
        #         if cfg['test']:
        #             print(guideids)
        #             print(dalignbedannot['guide: id'].unique())
                for gi in range(len(guideids)):
                    gid=guideids[gi]

                    dalignbedannoti=dalignbedannot.loc[dalignbedannot['guide: id']==gid,:]
                    if len(dalignbedannoti.shape)==1:
                        dalignbedannoti=pd.DataFrame(dalignbedannoti).T
        #             if cfg['test']:
        #                 df2info(dalignbedannoti)
        #                 print(cdhb)
                    for col in ['types','gene names','gene ids','transcript ids','protein ids','exon ids']:    
                #         print(";".join(dannoti[col].fillna('nan').tolist()))
                        daggbyguide.loc[gid,col]=";".join(np.unique(dalignbedannoti[col].fillna('nan').tolist()))
        #---
#                 if cfg['test']:
#                     df2info(daggbyguide)
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
                daggbyguide.to_csv(dofftargetsp,sep='\t')
                # daggbyguide.loc[:,['guide+PAM sequence','beditor score','beditor score (log10)','alternate alignments count',
                #              'id',
                #              'gene names',
                #              'gene ids',
                #              'transcript ids']].to_csv(dofftargetsp,sep='\t')

        import gc
        gc.collect()
