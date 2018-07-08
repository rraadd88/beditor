#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  


import sys
from os.path import exists,splitext,abspath,dirname,basename    
from os import makedirs
from glob import glob
import pandas as pd
import subprocess
import logging
from beditor.lib.io_sys import runbashcmd    
    
# GET INPTS    
def get_deps(cfg):
    """
    """
    depsd="%s/deps" % abspath(dirname(__file__)) 
    if not exists(depsd):
        makedirs(depsd)
    deps=['samtools','bedtools','bwa']
    ddeps=pd.DataFrame(columns=['local path','download link'],index=deps)
#     ddeps=ddeps.set_index('name')
    ddeps.index.name='dep'
    ddeps.loc[:,'local path']=['{}/{}'.format(depsd,dep) for dep in ddeps.index]

    dep='samtools'
    ddeps.loc[dep,'download link']='https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2'
    ddeps.loc[dep,'executable']='{}/samtools-1.7/samtools'.format(ddeps.loc[dep,'local path'])
    ddeps.loc[dep,'install']='cd {};make;'.format(dirname(ddeps.loc[dep,'executable']))

    dep='bedtools'
    ddeps.loc[dep,'download link']='https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz'
    ddeps.loc[dep,'executable']='{}/bedtools2/bin/bedtools'.format(ddeps.loc[dep,'local path'])
    ddeps.loc[dep,'install']='cd {}/../;make;'.format(dirname(ddeps.loc[dep,'executable']))

    dep='bwa'
    ddeps.loc[dep,'download link']='https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2'
    ddeps.loc[dep,'executable']='{}/bwa-0.7.17/bwa'.format(ddeps.loc[dep,'local path'])
    ddeps.loc[dep,'install']='cd {};make;'.format(dirname(ddeps.loc[dep,'executable']))
    
    ddeps.loc[:,'ext']=['.tar.{}'.format(ddeps.loc[n,'download link'].split('.tar.')[-1]) for n in ddeps.index]
    ddeps.to_csv("{}/deps.tsv".format(depsd),sep='\t')

    logp="%s/deps.log" % (depsd)
    with open(logp,'a') as logf:
        for dep in ddeps.index:            
            if not exists(ddeps.loc[dep,'executable']):
                logging.info("configuring: {} to {}".format(dep,ddeps.loc[dep,'executable']))
                link=ddeps.loc[dep,'download link']
                path=ddeps.loc[dep,'local path']            
                tarp='{}/{}'.format(path,basename(link))
                if not exists(tarp):                    
                    runbashcmd("wget -q %s --directory-prefix=%s" % (link,path),logf=logf)
                if not exists(dirname(ddeps.loc[dep,'executable'])):                    
                    if ddeps.loc[dep,'ext']=='.tar.bz2':
                        tarcom='xvjf'
                    elif ddeps.loc[dep,'ext']=='.tar.gz':
                        tarcom='zxvf'    
                    runbashcmd("tar {} {} -C {}".format(tarcom,tarp,path),logf=logf)
                runbashcmd(ddeps.loc[dep,'install'],logf=logf)
    #                 ddeps.loc[dep,'executable']='{}/.{}'.format(srcd,dep)
    #             break
                        
    logging.info("dependencies are installed!")
    for dep in ddeps.index:
        cfg[dep]=ddeps.loc[dep,'executable']
    return cfg