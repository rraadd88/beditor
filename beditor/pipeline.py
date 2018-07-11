#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname,splitext,basename,realpath
from os import makedirs
import argparse
import pkg_resources

import pandas as pd
from multiprocessing import Pool

from beditor.lib.io_strs import get_logger
logging=get_logger()
from beditor.configure import get_deps,get_genomes
from beditor.lib.io_sys import runbashcmd

# GET INPTS 

def pipeline_chunks(cfgp):
    from beditor.lib.get_seq import din2dseq
    from beditor.lib.get_mutations import dseq2dmutagenesis 
    from beditor.lib.make_guides import dseq2dguides
    from beditor.lib.get_specificity import dguides2offtargets

    logging.info('processing: '+cfgp)
    import yaml
    cfg=yaml.load(open(cfgp, 'r'))

    #project dir
    cfg['prj']=splitext(basename(cfgp))[0]
    if dirname(cfgp)!='':
        cfg['prjd']=dirname(cfgp)+'/'+cfg['prj']
    else:
        cfg['prjd']=cfg['prj']

    #datads
    cfg[0]=cfg['prjd']+'/00_input/'
    cfg[1]=cfg['prjd']+'/01_sequences/'
    cfg[2]=cfg['prjd']+'/02_mutagenesis/'
    cfg[3]=cfg['prjd']+'/03_guides/'
    cfg[4]=cfg['prjd']+'/04_offtargets/'

    #make dirs
    for i in range(5):
        if not exists(cfg[i]):
            makedirs(cfg[i])
    
    #backup inputs
    cfgoutp='{}/cfg.yml'.format(cfg[0])    
    dinoutp='{}/din.tsv'.format(cfg[0])    
    if not exists(cfgoutp) or cfg['force']:
        with open(cfgoutp, 'w') as f:
            yaml.dump(cfg, f, default_flow_style=False) 
    if not exists(dinoutp) or cfg['force']:
        from shutil import copyfile
        copyfile(cfg['dinp'], dinoutp)
#         din.to_csv(dinoutp,sep='\t')
    
    if not exists(cfg['prjd']):
        makedirs(cfg['prjd'])
    for i in range(0,4+1,1):
        if not exists(cfg[i]):
            makedirs(cfg[i])            
    if cfg['step']==None:
        stepall=True
    else:
        stepall=False
    if cfg['step']==1 or stepall:
        cfg['step']=1
        cfg=get_deps(cfg)
        din2dseq(cfg)
    if cfg['step']==2 or stepall:
        cfg['step']=2
        dseq2dmutagenesis(cfg)
    if cfg['step']==3 or stepall:
        cfg['step']=3
        dseq2dguides(cfg)
    if cfg['step']==4 or stepall:
        cfg['step']=4
        from beditor.configure import get_deps
        cfg=get_deps(cfg)
        dguides2offtargets(cfg)
    if 'datad' in cfg.keys():
        logging.info("Location of output data: {}".format(cfg['datad']))
        logging.info("Location of output plot: {}".format(cfg['plotd']))


def pipeline(cfgp,step=None,test=False,force=False):        

    import yaml
    cfg=yaml.load(open(cfgp, 'r'))
    #project dir
    cfg['prj']=splitext(basename(cfgp))[0]
    if dirname(cfgp)!='':
        cfg['prjd']=dirname(cfgp)+'/'+cfg['prj']
    else:
        cfg['prjd']=cfg['prj']

    #step
    cfg['step']=step
    # basics
    cfg['test']=test
    cfg['force']=force
    cfg['cfgp']=cfgp
    
    cfg=get_deps(cfg)
    cfg=get_genomes(cfg)

    #datads
    cfg[0]=cfg['prjd']+'/00_input/'
    cfg[1]=cfg['prjd']+'/01_sequences/'
    cfg[2]=cfg['prjd']+'/02_mutagenesis/'
    cfg[3]=cfg['prjd']+'/03_guides/'
    cfg[4]=cfg['prjd']+'/04_offtargets/'

    #make dirs
    for i in range(5):
        if not exists(cfg[i]):
            makedirs(cfg[i])
    #backup combo inputs
    cfgoutp='{}/cfg.yml'.format(cfg[0])    
    dinoutp='{}/din.tsv'.format(cfg[0])    
    if not exists(cfgoutp) or cfg['force']:
        with open(cfgoutp, 'w') as f:
            yaml.dump(cfg, f
#                       default_flow_style=False
                     ) 
    if not exists(dinoutp) or cfg['force']:
        from shutil import copyfile
        copyfile(cfg['dinp'], dinoutp)
        
    from beditor.lib.io_dfs import df2chucks
    din=pd.read_csv(cfg['dinp'],sep='\t')
    din=din.loc[:,['aminoacid: position','transcript: id']].drop_duplicates()

    chunkps=df2chucks(din,chunksize=100,
                      outd='{}/chunks'.format(cfg['prjd']),
                      fn='din',return_fmt='\t',
                      force=cfg['force'])

    chunkcfgps=[]
    for ci,cp in enumerate(chunkps):
        cfg_=cfg.copy()
        cfg_['dinp']=cp
        cfgp='{}/chunk{:08d}.yml'.format(cfg['prjd'],ci+1)    
        cfg_['cfgp']=cfgp
        if not exists(cfgp):
            with open(cfgp, 'w') as f:
                yaml.dump(cfg_, f, default_flow_style=False) 
        chunkcfgps.append(cfgp)
        
    if len(chunkps)!=0:
        if cfg['test']:
            pipeline_chunks(chunkcfgps[0])
        else:
            pool=Pool(processes=cfg['cores']) # T : get it from xls
            pool.map(pipeline_chunks, chunkcfgps)
            pool.close(); pool.join()         

#     pipeline_chunks(cfgp)
    logging.shutdown()

def main():
    """
    This runs all analysis steps in tandem.

    From bash command line,

    .. code-block:: text

        python path/to/beditor/pipeline.py cfg.json
        
    :param cfg.json: path to configurations.

    """
    version_info='%(prog)s v{version}'.format(version=pkg_resources.require("beditor")[0].version)
    parser = argparse.ArgumentParser(description=version_info)
    parser.add_argument("cfg", help="path to project directory", 
                        action="store", default=False)    
    parser.add_argument("--step", help="1: get seqeucnces,\n2: get possible strategies,\n3: make guides,\n 4: identify offtargets", dest="step", 
                        type=float,action="store", choices=[1,2,3,4],default=None)  
    parser.add_argument("--test", help="Debug mode on", dest="test", 
                        action='store_true', default=False)    
    parser.add_argument("--force", help="Debug mode on", dest="force", 
                        action='store_true', default=False)    
    parser.add_argument('-v','--version', action='version',version=version_info)
#    parser.add_argument('-h', '--help', action='help', #default=argparse.SUPPRESS,
#                    help='Show this help message and exit. \n Version info: %s' % version_info)
    args = parser.parse_args()
    logging.info("start")
    pipeline(args.cfg,step=args.step,
        test=args.test,force=args.force)

if __name__ == '__main__':
    main()
