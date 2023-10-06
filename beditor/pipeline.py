#!usr/bin/python
"""Workflow management"""

import logging
import sys
from os.path import exists,splitext,dirname,splitext,basename,abspath
from os import makedirs
import argparse
import yaml 
import pkg_resources
import numpy as np
import pandas as pd
from multiprocessing import Pool

from beditor.configure import get_genomes,validcfg,validinput
from beditor.lib.io_dfs import to_table, fhs2data_combo_appended
from beditor.lib.io_strs import get_datetime
from beditor.lib.io_nums import str2num
from beditor.lib.global_vars import stepi2name

pd.options.mode.chained_assignment = None

def pipeline_chunks(
    cfgp=None,
    cfg=None,
    ):
    """
    Runs indivudual chunk.

    :param cfgp: path to configuration file. 
    :param cfg: configuration dict
    :returns:
    """
    if cfg is None: 
        cfg=yaml.load(open(cfgp, 'r'))
    #datads
    cfg[0]=cfg['prjd']+'/00_input/'
    cfg[1]=cfg['prjd']+'/01_sequences/'
    cfg[2]=cfg['prjd']+'/02_mutagenesis/'
    cfg[3]=cfg['prjd']+'/03_guides/'
    cfg[4]=cfg['prjd']+'/04_offtargets/'

    # back compatibility
    if not 'make_control_neg' in cfg:
        cfg['make_control_neg']=False
    if not 'make_control_pos' in cfg:
        cfg['make_control_pos']=False        
        
    if not ('step2ignore' in cfg):
        cfg['step2ignore']=None            
    if not cfg['step2ignore'] is None:
        step_last=cfg['step2ignore']-1
    else:
        step_last=4

    dstep_last_outputp=f"{cfg[step_last]}/d{cfg[step_last].replace('/','').split('_')[-1]}.tsv"
    if not exists(dstep_last_outputp):
        logging.info(f"{get_datetime()}: processing: {basename(cfg['prjd'])}")
        logging.info(f"processing: {cfg['prjd']}")
    #     deps and genome are only needed if running step =1 or 4
        cfg['step2ignoredl']=[2,3]
        if not cfg['step'] in cfg['step2ignoredl']:
            if cfg['test']:
                logging.info('installing dependencies')
            if cfg['test']:
                logging.info('installing genomes')
            cfg=get_genomes(cfg)

        #make dirs
        for i in range(5):
            if not exists(cfg[i]):
                makedirs(cfg[i])

        #backup inputs
        cfgoutp=f'{cfg[0]}/cfg.yml'    
        dinoutp=f'{cfg[0]}/dinput.tsv'    
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
        stepi2time={}
        stepi2time[0]=get_datetime(ret_obj=True)
        if (cfg['step']==1 or stepall)  and cfg['step2ignore']!=1:
            from beditor.lib.get_seq import din2dseq
            cfg['step']=1
            din2dseq(cfg)
            stepi2time[1]=get_datetime(ret_obj=True)            
        if (cfg['step']==2 or stepall)  and cfg['step2ignore']!=2:
            from beditor.lib.get_mutations import dseq2dmutagenesis 
            cfg['step']=2
            dseq2dmutagenesis(cfg)
            stepi2time[2]=get_datetime(ret_obj=True)            
        if (cfg['step']==3 or stepall)  and cfg['step2ignore']!=3:
            from beditor.lib.make_guides import dseq2dguides
            cfg['step']=3
            dseq2dguides(cfg)
            stepi2time[3]=get_datetime(ret_obj=True)
        if (cfg['step']==4 or stepall)  and cfg['step2ignore']!=4:
            from beditor.lib.get_specificity import dguides2offtargets
            cfg['step']=4
            dguides2offtargets(cfg)
            stepi2time[4]=get_datetime(ret_obj=True)
        
        for stepi in range(1,5,1):
            if (stepi in stepi2time) and (stepi-1 in stepi2time):
                logging.info(f"time taken by step#{stepi:02d} = {str(stepi2time[stepi]-stepi2time[stepi-1])}")

def collect_chunks(cfg,chunkcfgps):
    """
    Collects analysed chunks

    :param cfg: main configuration dict.
    :param chunkcfgps: paths to all configuration files of chunks
    """    
    logging.info(f"{get_datetime()}: collecting chunks")    
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
                    if step==3:
                        # save the control guides
                        for control_type in ['pos','neg']:
                            dguides_controlp=f"{cfg['prjd']}/05_output/dguides.tsv.{control_type}_control.tsv"
                            if cfg[f'make_control_{control_type}']:
                                if not exists(dguides_controlp):
                                    from beditor.output import collectchuckfiles
                                    from beditor.lib.global_vars import stepi2colsoutput
                                    dguides_control=collectchuckfiles(cfg, fpinchunk=f'03_guides/dguides.tsv.{control_type}_control.tsv')
                                    if not dguides_control is None:
                                        dguides_control=dguides_control.loc[:,list(set(dguides_control.columns).intersection(stepi2colsoutput[3]))]    
                                        to_table(dguides_control,dguides_controlp)
                else:
                    logging.warning(f"no chunks found for step {step}: {stepi2name[step]}")                    
            else:
                logging.warning(f"no chunks found for step {step}: {stepi2name[step]}")


from glob import glob
from beditor.lib.plot_res import plot_vizbysteps
def make_outputs(cfg,plot=True):
    """
    Cobines stepwise analysis files into a pretty table.

    :param cfg: main configuration dict
    :param plot: if True creates visualizations
    """
    logging.info(f"{get_datetime()}: generating outputs")        
    from beditor.lib.global_vars import stepi2colsoutput
    prjd=cfg['prjd']
    #make one output table and stepwise plots
    datad=f"{prjd}/05_output"
    makedirs(datad, exist_ok=True)
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
                    dstep=del_Unnamed(pd.read_table(dstepp,keep_default_na=False))
                    if 'reverse_mutations' in cfg:
                        if cfg['reverse_mutations']:
                            if stepi==0:
                                continue
                    colsoutput=[col for col in colsoutput if col in dstep] 
                    dstep=dstep.loc[:,colsoutput]
                    if len(dstep)!=0:
                        dstep=dstep.drop_duplicates()
                    if not 'doutput' in locals():
                        doutput=dstep.copy()
                        del dstep
                    else:
                        cols_on=list(set(doutput.columns.tolist()).intersection(dstep.columns.tolist()))
                        if len(cols_on)!=0:         
                            doutput=pd.merge(doutput,dstep,on=cols_on,how='left')
                        else:
                            logging.error(f'output of step {stepi-1} or {stepi} are missing.')
                            return None
                        del dstep
        if cfg['mutation_format']=='nucleotide':
            doutput=doutput.drop([c for c in doutput if (('codon' in c) or ('amino' in c) or ('transcript' in c))],axis=1)
        if len(doutput)!=0 and 'guide+PAM sequence' in doutput:
            from beditor.lib.io_seqs import get_polyt_length
            doutput['length of polyT stretch']=doutput['guide+PAM sequence'].apply(lambda x : get_polyt_length(x))
        makedirs(dirname(doutputp),exist_ok=True)
        doutput.to_csv(doutputp,sep='\t')
    else:
        doutput=pd.read_table(doutputp,keep_default_na=False)
    # plot
    if plot:
        plot_vizbysteps(cfg)
    logging.info(f"Outputs are located at {datad}")
    return doutput 

def pipeline(
    cfgp,
    step=None,
    test=False,
    force=False,
    ):
    """
    Runs steps of the analysis workflow in tandem.
    
    :param cfgp: path to configuration file
    :param step: step number
    :param test: if True uses only one core, linear processing with verbose allowed
    :param force: if True overwrites outputs
    """
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
        logging.info(cfg)
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
    
    #back compatible
    deps=['samtools','bedtools','bwa']
    if not 'deps' in cfg:
        cfg['deps']=deps
    for dep in deps:
        if not dep in cfg:
            cfg[dep]=dep

    # defaults
    if not 'gui' in cfg:
        cfg['gui']=False
    if not 'chunksize' in cfg:
        cfg['chunksize']=200        
    if not 'max_subs_per_codon' in cfg:
        cfg['max_subs_per_codon']=1      
        
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

    if not 'reverse_mutations' in cfg:
        cfg['reverse_mutations']=False
    elif cfg['reverse_mutations'] is None:
        cfg['reverse_mutations']=False

    if not 'make_control_neg' in cfg:
        cfg['make_control_neg']=False
    if not 'make_control_pos' in cfg:
        cfg['make_control_pos']=False
        
    if cfg['make_control_pos']:
        cfg['keep_mutation_nonsense']=True

    if not 'step2ignore' in cfg:
        cfg['step2ignore']=None

    # deps and genome are only needed if running step =1 or 4
    cfg['step2ignoredl']=[2,3]
    if not cfg['step'] in cfg['step2ignoredl']:
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
    # prep be and pam
    cfg=validcfg(cfg,outcfg=True)
            
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
        din=pd.read_csv(cfg['dinp'],sep='\t',low_memory=False,keep_default_na=False)
        if not validinput(cfg,din):
            logging.error(f"input mutation file {cfg['dinp']} is not valid.")
            sys.exit(1)
        
        din=din.drop_duplicates()
        if not exists(dinoutp) or cfg['force']:
            if not 'amino acid mutation' in din:
                din.to_csv(dinoutp,sep='\t')
        cfg_=cfg.copy()
        cfg_['step']=1 #gotta run step 1 isolated because of memory buid up otherwise
        pipeline_chunks(cfg=cfg_)

        dseq=pd.read_csv(f"{cfg[1]}/dsequences.tsv",sep='\t',low_memory=False,keep_default_na=False)
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
        logging.info(f"{get_datetime()}: processing: {len(chunkcfgps)} chunks.")
        if cfg['test']:
            for chunkcfgp in chunkcfgps:
                pipeline_chunks(cfgp=chunkcfgp)
        else:
            pool=Pool(processes=cfg['cores']) # T : get it from xls
            pool.map(pipeline_chunks, chunkcfgps)
            pool.close(); pool.join()         
            collect_chunks(cfg,chunkcfgps)
            # get_outputs
#             if cfg['step2ignore'] is None:
            _=make_outputs(cfg)        
    else:
        pipeline_chunks(cfg=cfg)
        if not '/chunk' in cfgp:
#             if cfg['step2ignore'] is None:
            _=make_outputs(cfg)        

#     pipeline_chunks(cfgp)
    logging.shutdown()

def main():
    """
    Provides command-line inputs to the pipeline.

    For checking the command-lineinputs,

    .. code-block:: text

        beditor --help

    """
    version_info='%(prog)s v{version}'.format(version=pkg_resources.require("beditor")[0].version)
    parser = argparse.ArgumentParser(description=version_info)
    parser.add_argument("--cfg", help="path to configuration file in YAML format.", 
                        action="store", default=None, dest="cfg",type=str,
                       )    
    parser.add_argument("--step", help="1: Get genomic loci flanking the target site,\n2: Get possible mutagenesis strategies,\n3: Design guides,\n 4: Check offtarget-effects \n else all the steps are run in tandem.", dest="step", 
                        type=float,action="store", choices=[1,2,3,4],default=None)  
    parser.add_argument("--test", help="Debug mode on", dest="test", 
                        action='store_true', default=False)    
    parser.add_argument("--force", help="Overwrite existing outputs.", dest="force", 
                        action='store_true', default=False)    
    parser.add_argument('-v','--version', action='version',version=version_info)
    args = parser.parse_args()

    logging.info(args)
    
    if not args.cfg is None:
        from beditor.lib.io_strs import get_logger
        if args.test:
            level=logging.INFO
        else: 
            level=logging.ERROR
        logp=get_logger(program='beditor',
                   argv=list(vars(args).values()),
                   level=level,
                   dp=None)
        time_ini=get_datetime(ret_obj=True)
        logging.info(f"start. log file: {logp}")
        logging.info(f"start. log file: {logp}")
        pipeline(args.cfg,step=args.step,
            test=args.test,force=args.force)
        logging.info(f'end. time taken={str(get_datetime(ret_obj=True)-time_ini)}')
    else:
        from .gui import gui
        gui(test=args.test)
#     else:
#         parser.print_help()
#         sys.exit(1)        

if __name__ == '__main__':
    main()
