#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname,splitext,basename,realpath
from os import makedirs
import argparse
import pkg_resources

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

from multiprocessing import Pool

import logging
from beditor.configure import get_deps,get_genomes
from beditor.lib.io_sys import runbashcmd

# GET INPTS 

def pipeline_chunks(cfgp):
    # from beditor.configure import get_deps

    import yaml
    cfg=yaml.load(open(cfgp, 'r'))
    logging.info('processing: '+cfgp)

#     print(cfg)    
#     deps and genome are only needed if running step =1 or 4
    cfg['step2ignoredl']=[2,3]
    if not cfg['step'] in cfg['step2ignoredl']:
        print('installing dependencies')
        cfg=get_deps(cfg)
        print('installing genomes')
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
    
    #backup inputs
    cfgoutp='{}/cfg.yml'.format(cfg[0])    
    dinoutp='{}/dinput.tsv'.format(cfg[0])    
    if not exists(cfgoutp) or cfg['force']:
        with open(cfgoutp, 'w') as f:
            yaml.dump(cfg, f, default_flow_style=False)
    logging.info(cfg)
    if cfg['step']==None or cfg['step']==1:
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
    # print(cfg['step'],stepall)
    print('processing: '+cfgp)
    if cfg['step']==1 or stepall:
        from beditor.lib.get_seq import din2dseq
        cfg['step']=1
        din2dseq(cfg)
    if cfg['step']==2 or stepall:
        from beditor.lib.get_mutations import dseq2dmutagenesis 
        cfg['step']=2
        dseq2dmutagenesis(cfg)
    if cfg['step']==3 or stepall:
        from beditor.lib.make_guides import dseq2dguides
        cfg['step']=3
        dseq2dguides(cfg)
    if cfg['step']==4 or stepall:
        from beditor.lib.get_specificity import dguides2offtargets
        cfg['step']=4
        dguides2offtargets(cfg)
    if 'datad' in cfg.keys():
        print("Location of output data: {}".format(cfg['datad']))
        print("Location of output plot: {}".format(cfg['plotd']))
    # if not exists(f"{cfg['prjd']}/04_offtargets/dofftargets.tsv"):
    # else:
    #     print(f"skipped: {cfg['cfgp']}")

from beditor.lib.io_dfs import fhs2data_combo_appended
from beditor.lib.global_vars import stepi2name
from beditor.lib.io_nums import str2num
from os.path import basename
def collect_chunks(cfg,chunkcfgps):
    """
    #collects chunks
    """    
    for step in stepi2name.keys():
        doutp=f"{cfg['prjd']}/{step:02d}_{stepi2name[step]}/d{stepi2name[step]}.tsv"
        if not exists(doutp) or cfg['force']:
            dps=[] 
            outd=cfg[step]
            for chunkcfgp in chunkcfgps:
                chunkprjd=chunkcfgp.replace('.yml','')
                dps.append(f"{chunkprjd}/{step:02d}_{stepi2name[step]}/d{stepi2name[step]}.tsv")
            if len(dps)!=0:
                dout=fhs2data_combo_appended(dps,sep='\t',
                                             labels=[str2num(basename(p)) for p in chunkcfgps],
                                             labels_coln='chunk#')
                dout.to_csv(doutp,sep='\t')
                del dout

from glob import glob
from beditor.lib.plot_res import plot_vizbysteps
def make_outputs(cfg,plot=True):
    stepi2cols={0: ['transcript: id','aminoacid: position','amino acid mutation'],
            1: ['transcript: id','aminoacid: position','aminoacid: wild-type'],
            3: ['transcript: id','aminoacid: position','aminoacid: wild-type','amino acid mutation','guide: id','guide+PAM sequence'],
            4: ['guide: id','beditor score','alternate alignments count','CFD score'],
           }
    prjd=cfg['prjd']
    #make one output table and stepwise plots
    datad=f"{prjd}/05_output"
    #table
    doutputp=f"{datad}/doutput.tsv" #FIXME if steps are added
    if not exists(doutputp) or cfg['force']:
        from beditor.lib.io_dfs import del_Unnamed
        if 'doutput' in locals():
            del doutput
        for stepi in range(5):
            if stepi!=2:
                dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
                if exists(dstepp):
                    logging.info(f'combining {stepi}')
                    dstep=del_Unnamed(pd.read_table(dstepp)).loc[:,stepi2cols[stepi]].drop_duplicates()
                    if not 'doutput' in locals():
                        doutput=dstep.copy()
                        del dstep
                    else:
                        cols_on=list(set(doutput.columns.tolist()).intersection(dstep.columns.tolist()))
                        if len(cols_on)!=0:                    
                            doutput=pd.merge(doutput,dstep,on=cols_on,how='left')
                        else:
                            logging.error(f'output of step {stepi-1} is missing.')
        #         if stepi==4:
        #             break
        makedirs(dirname(doutputp),exist_ok=True)
        doutput.to_csv(doutputp,sep='\t')
    else:
        doutput=pd.read_table(doutputp)
    # plot
    if plot:
        plot_vizbysteps(cfg)
    return doutput 

def pipeline(cfgp,step=None,test=False,force=False):        

    import yaml
    from glob import glob

    cfg=yaml.load(open(cfgp, 'r'))
    # check inputs
    if not exists(cfg['dinp']):
        logging.error(f"input file {cfg['dinp']} is not found.")
        sys.exit(1)
    if (cfg['mutations']=='substitutions'):    
        if not exists(cfg['dsubmap_preferred_path']):
            logging.critical(f"dsubmap_preferred_path is {cfg['dsubmap_preferred_path']}")
            logging.critical(cfg)
            sys.exit(1)

    # get names right
    import pyensembl
    cfg['host']=pyensembl.species.normalize_species_name(cfg['host'])        
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
    
#     deps and genome are only needed if running step =1 or 4
    cfg['step2ignoredl']=[2,3]
    if not cfg['step'] in cfg['step2ignoredl']:
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
    dinoutp='{}/dinput.tsv'.format(cfg[0])    
    if not exists(cfgoutp) or cfg['force']:
        with open(cfgoutp, 'w') as f:
            yaml.dump(cfg, f
#                       default_flow_style=False
                     ) 
    if (not '/chunk' in cfgp) and (step==1 or (step is None)):
        from beditor.lib.io_dfs import df2chucks
        din=pd.read_csv(dinoutp,sep='\t')
        din=din.drop_duplicates()
        if not exists(dinoutp) or cfg['force']:
            if not 'amino acid mutation' in din:
                din.to_csv(dinoutp,sep='\t')
        chunkps=df2chucks(din,chunksize=cfg['chunksize'],
                          outd='{}/chunks'.format(cfg['prjd']),
                          fn='din',return_fmt='\t',
                          force=cfg['force'])
        chunkcfgps=[]
        for ci,cp in enumerate(chunkps):
            cfg_=cfg.copy()
            cfg_['dinp']=cp
            chunkcfgp=f"{cfg['prjd']}/chunks/chunk{(ci+1):08}.yml"    
            cfg_['cfgp']=chunkcfgp
            
            #project dir
            cfg_['prj']=splitext(basename(cfg_['cfgp']))[0]
#             if dirname(chunkcfgp)!='':
            cfg_['prjd']=f"{dirname(chunkcfgp)}/{cfg_['prj']}"
#             else:
#                 cfg_['prjd']=f"./chunks/{cfg_['prj']}"
            cfg_['step']=None    
            if cfg['test'] and cfg['force']:
                cfg_['test']=True 
                cfg_['force']=True    
            elif cfg['test'] and not cfg['force']:
                cfg_['test']=True 
                cfg_['force']=False
            elif not cfg['test'] and cfg['force']:
                cfg_['test']=False 
                cfg_['force']=True
            else:
                cfg_['test']=False    
                cfg_['force']=False    
            if (not exists(chunkcfgp)) or cfg['force'] or cfg['test']:
                with open(chunkcfgp, 'w') as f:
                    print(f"created {chunkcfgp}")
                    yaml.dump(cfg_, f, default_flow_style=False) 
            chunkcfgps.append(chunkcfgp)
        # sys.exit(1)
    else:
        chunkcfgps=glob('{}/chunk*.yml'.format(cfg['prjd']))
    
    chunkcfgps=np.sort(chunkcfgps)
    if len(chunkcfgps)!=0 and (not '/chunk' in cfgp):
        if cfg['test']:
            pipeline_chunks(chunkcfgps[0])
        else:
            pool=Pool(processes=cfg['cores']) # T : get it from xls
            pool.map(pipeline_chunks, chunkcfgps)
            pool.close(); pool.join()         
            collect_chunks(cfg,chunkcfgps)
    else:
        pipeline_chunks(cfgoutp)

    # get_outputs
    _=make_outputs(cfg)        

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
    parser.add_argument("--cfg", help="path to project directory", 
                        action="store", default=False)    
    parser.add_argument("--step", help="1: get seqeucnces,\n2: get possible strategies,\n3: make guides,\n 4: identify offtargets \n else all the steps in tandem.", dest="step", 
                        type=float,action="store", choices=[1,2,3,4],default=None)  
    parser.add_argument("--list", help="list something", dest="lister", 
                        action="store", default=None)
    parser.add_argument("--test", help="Debug mode on", dest="test", 
                        action='store_true', default=False)    
    parser.add_argument("--force", help="Overwrite existing outputs.", dest="force", 
                        action='store_true', default=False)    
    parser.add_argument('-v','--version', action='version',version=version_info)
#    parser.add_argument('-h', '--help', action='help', #default=argparse.SUPPRESS,
#                    help='Show this help message and exit. \n Version info: %s' % version_info)
    args = parser.parse_args()

    lists=['pams','editors']
    if args.lister is not None:
        if args.lister in lists:
            from .lib.io_dfs import del_Unnamed
            if args.lister.lower()=='pams':
                d=pd.read_table(f"{dirname(realpath(__file__))}/data/dpam.tsv")
                print(del_Unnamed(d.loc[:,['PAM','Description']]).set_index('PAM'))                
            elif args.lister.lower()=='editors':
                d=pd.read_table(f"{dirname(realpath(__file__))}/data/dBEs.tsv")
                print(del_Unnamed(d.loc[(d['strand']=='+'),['method', 'nucleotide', 'nucleotide mutation']]).set_index('method'))
        else:
            logging.error("args.lister should be one these {','.join(lists)}")  
    else:

        from beditor.lib.io_strs import get_logger
        if args.test:
            level=logging.INFO
        else: 
            level=logging.ERROR
        logp=get_logger(program='beditor',
                   argv=list(vars(args).values()),
                   level=level,
                   dp=None)
        
        logging.info(f"start\nlog file: {logp}")
        print(f"start\nlog file: {logp}")
        pipeline(args.cfg,step=args.step,
            test=args.test,force=args.force)

if __name__ == '__main__':
    main()
