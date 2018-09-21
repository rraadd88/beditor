#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
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
from beditor.lib.io_strs import get_datetime
import yaml 

def pipeline_chunks(cfgp=None,cfg=None):
    if cfg is None: 
        cfg=yaml.load(open(cfgp, 'r'))
    #datads
    cfg[0]=cfg['prjd']+'/00_input/'
    cfg[1]=cfg['prjd']+'/01_sequences/'
    cfg[2]=cfg['prjd']+'/02_mutagenesis/'
    cfg[3]=cfg['prjd']+'/03_guides/'
    cfg[4]=cfg['prjd']+'/04_offtargets/'

    if not ('step2ignore' in cfg):
        cfg['step2ignore']=None
        
    if not cfg['step2ignore'] is None:
        step_last=cfg['step2ignore']-1
    else:
        step_last=4

    dstep_last_outputp=f"{cfg[step_last]}/d{cfg[step_last].replace('/','').split('_')[-1]}.tsv"
    if not exists(dstep_last_outputp):
        print(f"{get_datetime()}: processing: {basename(cfg['prjd'])}")
        logging.info(f"processing: {cfg['prjd']}")
    #     deps and genome are only needed if running step =1 or 4
        cfg['step2ignoredl']=[2,3]
        if not cfg['step'] in cfg['step2ignoredl']:
            if cfg['test']:
                logging.info('installing dependencies')
            cfg=get_deps(cfg)
            if cfg['test']:
                logging.info('installing genomes')
            cfg=get_genomes(cfg)

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
        logging.info(cfgp)
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
        if (cfg['step']==1 or stepall)  and cfg['step2ignore']!=1:
            from beditor.lib.get_seq import din2dseq
            cfg['step']=1
            din2dseq(cfg)
        if (cfg['step']==2 or stepall)  and cfg['step2ignore']!=2:
            from beditor.lib.get_mutations import dseq2dmutagenesis 
            cfg['step']=2
            dseq2dmutagenesis(cfg)
        if (cfg['step']==3 or stepall)  and cfg['step2ignore']!=3:
            from beditor.lib.make_guides import dseq2dguides
            cfg['step']=3
            dseq2dguides(cfg)
        if (cfg['step']==4 or stepall)  and cfg['step2ignore']!=4:
            from beditor.lib.get_specificity import dguides2offtargets
            cfg['step']=4
            dguides2offtargets(cfg)

from beditor.lib.io_dfs import fhs2data_combo_appended
from beditor.lib.global_vars import stepi2name
from beditor.lib.io_nums import str2num
from os.path import basename
def collect_chunks(cfg,chunkcfgps):
    """
    #collects chunks
    """    
    print(f"{get_datetime()}: collecting chunks")    
    for step in stepi2name.keys():
        doutp=f"{cfg['prjd']}/{step:02d}_{stepi2name[step]}/d{stepi2name[step]}.tsv"
        if not exists(doutp) or cfg['force']:
            dps=[] 
            outd=cfg[step]
            for chunkcfgp in chunkcfgps:
                chunkprjd=chunkcfgp.replace('.yml','')
                dps.append(f"{chunkprjd}/{step:02d}_{stepi2name[step]}/d{stepi2name[step]}.tsv")
#             dps=[p if exists(p) else np.nan for p in dps]
            tdpcp=[(dp,cp) for dp,cp in zip(dps,chunkcfgps) if exists(dp)]            
            if len(tdpcp)!=0:
                dps_,chunkcfgps_=zip(*tdpcp)
                if len(dps_)!=0:                
                    dout=fhs2data_combo_appended(dps_,sep='\t',
                                                 labels=[str2num(basename(p)) for p in chunkcfgps_],
                                                 labels_coln='chunk#')
                    dout.to_csv(doutp,sep='\t')
                    del dout
                else:
                    logging.warning(f"no chunks found for step {step}: {stepi2name[step]}")                    
            else:
                logging.warning(f"no chunks found for step {step}: {stepi2name[step]}")

def collectchuckfiles(cfg,fpinchunk,force=False):
    from beditor.lib.io_dfs import fhs2data_combo_appended
    from beditor.lib.global_vars import stepi2name
    from beditor.lib.io_nums import str2num
    from os.path import basename
    from glob import glob
    doutp=f"{cfg['prjd']}/{fpinchunk}"
    if not exists(doutp) or force:
        dps_=glob(f"{cfg['prjd']}/chunks/chunk0*/{fpinchunk}")
#         print(dps_)
        dout=fhs2data_combo_appended(dps_,sep='\t',
                                     labels_coln='chunk#')
        makedirs(dirname(doutp),exist_ok=True)
        dout.to_csv(doutp,sep='\t')
    else:
        dout=pd.read_table(doutp)        
    return dout 

from glob import glob
from beditor.lib.plot_res import plot_vizbysteps
def make_outputs(cfg,plot=True):
    print(f"{get_datetime()}: generating outputs")        
    from beditor.lib.global_vars import stepi2colsoutput
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
            if stepi!=2 and cfg['step2ignore']!=stepi:
                dstepp=f"{cfg[stepi]}/d{cfg[stepi].replace('/','').split('_')[-1]}.tsv"
                if exists(dstepp):
                    logging.info(f'combining {stepi}')
                    colsoutput=stepi2colsoutput[stepi]
                    dstep=del_Unnamed(pd.read_table(dstepp))
                    if 'reverse_mutations' in cfg:
                        if cfg['reverse_mutations']:
                            if stepi==0:
                                continue
                    colsoutput=[col for col in colsoutput if col in dstep] 
                    dstep=dstep.loc[:,colsoutput].drop_duplicates()
                    if not 'doutput' in locals():
                        doutput=dstep.copy()
                        del dstep
                    else:
                        cols_on=list(set(doutput.columns.tolist()).intersection(dstep.columns.tolist()))
#                         print('left',doutput.columns.tolist())
#                         print('right',dstep.columns.tolist())
#                         print('common',cols_on)
                        if len(cols_on)!=0:         
                            doutput=pd.merge(doutput,dstep,on=cols_on,how='left')
                        else:
                            logging.error(f'output of step {stepi-1} or {stepi} are missing.')
                        del dstep
        if cfg['mutation_format']=='nucleotide':
            doutput=doutput.drop([c for c in doutput if (('codon' in c) or ('amino' in c) or ('transcript' in c))],axis=1)
        makedirs(dirname(doutputp),exist_ok=True)
        doutput.to_csv(doutputp,sep='\t')
    else:
        doutput=pd.read_table(doutputp)
    # plot
    if plot:
        plot_vizbysteps(cfg)
    logging.info(f"Outputs are located at {datad}")
    return doutput 

def validcfg(cfg): 
    from beditor.lib.global_vars import cfgoption2allowed,cfgoption2reguired
    opt_validity=[]
    for opt in ['mutations','mutation_format']:
        opts=cfgoption2allowed[opt]
        if cfg[opt] in opts:
            opt_validity.append(True)
            if opt in cfgoption2reguired:
                for opt_opt in cfgoption2reguired[opt]:
                    opt_opt_opt=cfgoption2reguired[opt][opt_opt]
                    if cfg[opt]==opt_opt: 
                        if (opt_opt_opt in cfg):
                            opt_validity.append(True)
                            if 'path' in opt_opt_opt:
                                if not cfg[opt_opt_opt] is None:
                                    if exists(cfg[opt_opt_opt]):
                                        opt_validity.append(True)
                                    else:
                                        opt_validity.append(False)
                                        logging.error(f"{opt_opt_opt}:{cfg[opt_opt_opt]} is not found.")
                                else:
                                    opt_validity.append(False)
                                    logging.error(f"path {opt_opt_opt} is {cfg[opt_opt_opt]}.")
                        else:
                            opt_validity.append(False)
                            logging.error(f"{opt_opt} is not an option variable")
        else:
            opt_validity.append(False)
            logging.error(f"invalid option: {cfg[opt]} is not in [{','.join([s if not s is None else 'None' for s in opts])}]")
    return all(opt_validity)

def validinput(cfg,din): 
    from beditor.lib.global_vars import mutation_format2cols
    opt_validity=[]
    for col in mutation_format2cols[cfg['mutation_format']]:
        if col in din:
            opt_validity.append(True)
        else:
            opt_validity.append(False)
            logging.error(f"invalid column name: {col} is not in [{','.join(mutation_format2cols[cfg['mutation_format']])}]")
    return all(opt_validity)

def pipeline(cfgp,step=None,test=False,force=False):
    import yaml
    from glob import glob
    cfgp=abspath(cfgp)
    cfg=yaml.load(open(cfgp, 'r'))
    
    # check inputs
    if not exists(cfg['dinp']):
        logging.error(f"input file {cfg['dinp']} is not found.")
        sys.exit(1)

    if not validcfg(cfg):
        logging.error(f"configuration file {cfgp} is not valid.")
        print(cfg)
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
    cfg['prjd']=dirname(cfgp)+'/'+cfg['prj']
    cfg['beditor']=abspath(__file__)

    #step
    cfg['step']=step
    # basics
    cfg['test']=test
    cfg['force']=force
    cfg['cfgp']=cfgp
    # step 04 offtargets
    if not 'mismatches_max' in cfg:
        cfg['mismatches_max']=2
    elif cfg['mismatches_max'] is None:
        cfg['mismatches_max']=2
        # logging.info(f"setting mismatches_max to {cfg['mismatches_max']}") 
    if not 'reverse_mutations' in cfg:
        cfg['reverse_mutations']=False
    elif cfg['reverse_mutations'] is None:
        cfg['reverse_mutations']=False

    if not 'step2ignore' in cfg:
            cfg['step2ignore']=None

    # deps and genome are only needed if running step =1 or 4
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
    cfgoutp=f'{cfg[0]}/cfg.yml'
    dinoutp=f'{cfg[0]}/dinput.tsv'    
    if not exists(cfgoutp) or cfg['force']:
        with open(cfgoutp, 'w') as f:
            yaml.dump(cfg, f
#                       default_flow_style=False
                     ) 
    if (not '/chunk' in cfgp) and (step==1 or (step is None)):
        from beditor.lib.io_dfs import df2chucks
        din=pd.read_csv(cfg['dinp'],sep='\t',low_memory=False)
        if not validinput(cfg,din):
            logging.error(f"configuration file {cfgp} is not valid.")
            sys.exit(1)
        
        din=din.drop_duplicates()
        if not exists(dinoutp) or cfg['force']:
            if not 'amino acid mutation' in din:
                din.to_csv(dinoutp,sep='\t')
        cfg_=cfg.copy()
        cfg_['step']=1 #gotta run step 1 isolated because of memory buid up otherwise
        pipeline_chunks(cfg=cfg_)

        dseq=pd.read_csv(f"{cfg[1]}/dsequences.tsv",sep='\t',low_memory=False)
        chunkps=df2chucks(dseq,chunksize=cfg['chunksize'],
                          outd=f"{cfg['prjd']}/chunks",
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
                    yaml.dump(cfg_, f, default_flow_style=False) 
            chunkcfgps.append(chunkcfgp)
        # sys.exit(1)
    else:
        chunkcfgps=glob('{}/chunk*.yml'.format(cfg['prjd']))
    
    chunkcfgps=np.sort(chunkcfgps)
    if len(chunkcfgps)!=0 and (not '/chunk' in cfgp):
        print(f"{get_datetime()}: processing: {len(chunkcfgps)} chunks.")
        if cfg['test']:
            for chunkcfgp in chunkcfgps:
                pipeline_chunks(cfgp=chunkcfgp)
        else:
            pool=Pool(processes=cfg['cores']) # T : get it from xls
            pool.map(pipeline_chunks, chunkcfgps)
            pool.close(); pool.join()         
            collect_chunks(cfg,chunkcfgps)
            # get_outputs
            _=make_outputs(cfg)        
    else:
        pipeline_chunks(cfg=cfg)
        if not '/chunk' in cfgp:
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
