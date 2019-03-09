#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
from os import makedirs
from glob import glob
import pandas as pd
import subprocess
import logging
from beditor.lib.io_sys import runbashcmd    
from beditor.lib.io_dfs import * 

def get_genomeurls(name,release,test=False):
    name=name.lower().replace(' ','_')
    release=int(release)
    import ftplib
    host='ftp.ensembl.org'
    ftp=ftplib.FTP(host,user='anonymous', passwd='anonymous')
    contig_mito=['MTDNA','MITO','MT']    
    dtype=['dna','gff3']
    subds=[ f'fasta/{name}/dna/',f'gff3/{name}/']
    dtype2subd=dict(zip(dtype,subds))
    dtype2url={}
    for dtype in dtype2subd:
        subd=dtype2subd[dtype]
        ext=f'/pub/release-{int(release)}/{subd}'
#         if test:
#             print(ext)
        ftp.cwd(ext)
        if dtype=='dna':
            file_fmt='.fa.gz'        
            if test:
                print(len(ftp.nlst()))
            fns_=[p for p in ftp.nlst() if p.endswith(file_fmt) and ('dna_sm' in p) and (not 'nonchromosomal' in p) and (not '.primary_assembly.' in p) and (not '.chr.' in p)]
            if len(fns_)>1:
                fns_=[fn for fn in fns_ if (not 'dna_sm.toplevel' in fn) and (not '.alt.' in fn)]
                if test:
                    print(fns_)
                #contigs are present
                if sum(['dna_sm.chromosome.' in fn for fn in fns_])!=0:
                    #remove mito and contigs
                    fns=[]
                    for p in fns_:
                        contig=p.split('dna_sm.chromosome.')[1].split(file_fmt)[0]
                        if (not '.' in contig) and (not len(contig)>5) and (not contig in contig_mito):
                            fns.append(p)
                else:
                    brk
            else:
                fns=[fn for fn in fns_ if 'dna_sm.toplevel' in fn]
            if len(fns)==0:
                logging.warning(f'no genome files found: {dtype}')                
                print(ftp.nlst())
            urls=[f"ftp://{host}{ext}{fn}" for fn in fns]
            dtype2url[dtype]=urls
        elif dtype=='gff3':
            file_fmt=f'{release}.gff3.gz'                    
            print(file_fmt)
#             fns=[p for p in ftp.nlst()]
            fns=[p for p in ftp.nlst() if p.endswith(file_fmt) and (not 'chromosome' in p) and (not 'abinitio' in p) and (not 'patch' in p) and (not 'scaff' in p)  and (not '.chr.' in p)]
            print(fns)
            if len(fns)>1:
                logging.warning(f'two files instead of one found for {dtype}')
#                 if test:
#                     print(fns)                
            url=f"ftp://{host}{ext}{fns[0]}"
            dtype2url[dtype]=url
    if test:
        print(dtype2url)
    return dtype2url

def get_genomes(cfg):
    """
    Installs genomes
    
    :param cfg: configuration dict
    """
#     print('checking if genome is installed/ downloading if necessary.')
    logging.info('pyensembl: checking if genome is installed/ downloading if necessary.')
    runbashcmd(f"pyensembl install --reference-name {cfg['genomeassembly']} --release {cfg['genomerelease']} --species {cfg['host']}")

    # if 'step2ignore' in cfg:
    #     if cfg['step2ignore']==4:z`
    #         return cfg
        
    #download genome for step 5 specificity
    host_="_".join(s for s in cfg['host'].split('_')).capitalize()
    ensembl_fastad='pub/release-{}/fasta/{}/dna/'.format(cfg['genomerelease'],cfg['host'])
    genome_fastad='{}/{}'.format(dirname(realpath(__file__)),ensembl_fastad)
    cfg['genomep']='{}/genome.fa'.format(genome_fastad)
    if not exists(cfg['genomep']):
        logging.error(f"not found: {cfg['genomep']}")
        logging.info(f"downloading file: {cfg['genomep']}")
        #back compatible
        if not 'gui' in cfg:
            cfg['gui']=False        
        if (not '/test_beditor/' in cfg['cfgp']) or (not cfg['gui']):
            ifdlref = input("Download genome at {}?[Y/n]: ".format(genome_fastad))
        else:
            ifdlref='Y'
        if ifdlref=='Y':
        # #FIXME download contigs and cat and get index, sizes
            try:
                contigurls=get_genomeurls(cfg['host'],cfg['genomerelease'],
                                      test=False)['dna']
            except:
                contigurls=get_genomeurls(cfg['host'],cfg['genomerelease'],
                                      test=False)['dna']            
            logging.info(f"{len(contigurls)} contigs/chromosomes in the genome")
            logging.info(contigurls)
            for contigurl in contigurls:
                fn=basename(contigurl)
                fp=f'{ensembl_fastad}/{fn}'
                logging.info(f"downloading: {contigurl}")
                if not exists(fp):
                    cmd=f"wget -q -x -nH {contigurl} -P {dirname(realpath(__file__))}"
                    runbashcmd(cmd,test=cfg['test'])
            # make the fa ready
            if not exists(cfg['genomep']):
                cmd='gunzip {}*.fa.gz;cat {}/*.fa > {}/genome.fa;'.format(genome_fastad,genome_fastad,genome_fastad)
                runbashcmd(cmd,test=cfg['test'])
        else:
            logging.error('abort')
            sys.exit(1)
    if not exists(cfg['genomep']+'.bwt'):
        cmd='{} index {}'.format(cfg['bwa'],cfg['genomep'])
        runbashcmd(cmd,test=cfg['test'])
    else:        
        logging.info('bwa index is present')
    if not exists(cfg['genomep']+'.fai'):
        cmd='{} faidx {}'.format(cfg['samtools'],cfg['genomep'])
        runbashcmd(cmd,test=cfg['test'])
    else:
        logging.info('samtools index is present')
    if not exists(cfg['genomep']+'.sizes'):
        cmd='cut -f1,2 {}.fai > {}.sizes'.format(cfg['genomep'],cfg['genomep'])            
        runbashcmd(cmd,test=cfg['test'])
    else:
        logging.info('sizes of contigs are present')
        
    #download gff3
    ensembl_gff3d='pub/release-{}/gff3/{}/'.format(cfg['genomerelease'],cfg['host'])    
    genome_gff3d=f'{dirname(realpath(__file__))}/{ensembl_gff3d}'
    cfg['genomegffp']='{}/genome.gff3'.format(genome_gff3d)
    if not exists(cfg['genomegffp']):
        logging.error(f"not found: {cfg['genomegffp']}")
        logging.info(f"downloading file: {cfg['genomegffp']}")
        
        if (not '/test_beditor/' in cfg['cfgp']) or (not cfg['gui']):
            ifdlref = input(f"Download genome annotations at {genome_gff3d}?[Y/n]: ")
        else:
            ifdlref = 'Y'
        if ifdlref=='Y':
        # #FIXME download contigs and cat and get index, sizes
            fn='{}.{}.{}.gff3.gz'.format(cfg['host'].capitalize(),cfg['genomeassembly'],cfg['genomerelease'])
            fp='{}/{}'.format(ensembl_gff3d,fn)
            try: 
                ensembl_gff3p=get_genomeurls(cfg['host'],cfg['genomerelease'],
                                                             test=False)['gff3']
            except:
                ensembl_gff3p=get_genomeurls(cfg['host'],cfg['genomerelease'],
                                                             test=False)['gff3']                
            if not exists(fp):
                cmd=f'wget -x -nH {ensembl_gff3p} -P {dirname(realpath(__file__))}'
                runbashcmd(cmd,test=cfg['test'])
            # move to genome.gff3
                cmd='cp {}/{} {}'.format(genome_gff3d,fn,cfg['genomegffp'])
                runbashcmd(cmd,test=cfg['test'])
        else:
            logging.error('abort')
            sys.exit(1)
    logging.info('genomes are installed!')
    return cfg

def get_distance_of_mutation_from_pam(pam_position='down',
                                        window_ini=4,
                                        window_end=8,
                                     guide_length=20):
    """
    down,20
    8:13
    4:17 
    down,21
    12:10
    3:19    
    """    
    if pam_position=='down':
        mn=guide_length-window_end+1
        mx=mn+window_end-window_ini
    elif pam_position=='up':
        mn=window_ini
        mx=window_end        
    return mn,mx

def get_distance_of_codon_from_PAM(d_pammin,d_pammax,pam_position):
    if pam_position=='down':
        mn,mx=d_pammin,d_pammax+2
    elif pam_position=='up':
        mn,mx=d_pammin+2,d_pammax
    return mn,mx

# validity of inputs
def validcfg(cfg,outcfg=False): 
    """
    Checks if configuration dict is valid i.e. contains all the required fields

    :param cfg: configuration dict
    """
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
    #rename
    if 'BEs' in cfg:
        cfg['BE names']=cfg['BEs']
        del cfg['BEs']
    if 'pams' in cfg:
        cfg['PAMs']=cfg['pams']
        del cfg['pams']
    #debug
    for option in ['BE name and PAM','PAMs','BE names']:
        if not (option in cfg): cfg[option]=None
    if (cfg['PAMs'] is None) and (cfg['BE names'] is None):
        if not cfg['BE name and PAM'] is None:
            if isinstance(cfg['BE name and PAM'],str):
                cfg['BE name and PAM']=cfg['BE name and PAM']
        else:
            opt_validity.append(False)       
            logging.error(f"invalid option: BE PAM not specified")
    elif (cfg['PAMs'] is None) or (cfg['BE names'] is None):
        opt_validity.append(False)       
        logging.error(f"invalid option: BE PAM not specified")     
    else:
        from itertools import product
        cfg['BE name and PAM']=[list(t) for t in product(cfg['BE names'],cfg['PAMs'])]
    BE_names=list(np.unique([t[0] for t in cfg['BE name and PAM']]))
    PAMs=list(np.unique([t[1] for t in cfg['BE name and PAM']]))
    cfg['BE names']=[str(s) for s in BE_names]
    cfg['PAMs']=[str(s) for s in PAMs]
    if 'guide length' in cfg:
        cfg['guide length']=[int(i) for i in  cfg['guide length']]
    if outcfg:
        #save run specific debepams
        cfg['dbepamsp']=f"{cfg[0]}/dbepams.tsv"
        dbepams=pd.read_table(f"{dirname(realpath(__file__))}/data/dbepams.tsv",keep_default_na=False)
        # check if gui added custom be pam. add that to the prj table 
        if cfg['gui']:
            # there is only one pair of be and pam added by gui at a time
            if (not cfg['BE names'][0] in dbepams['method'].tolist()) or (not cfg['PAMs'][0] in dbepams['PAM'].tolist()):
                # if custom BE and PAM
                # got a new be and/or pam
                # f"method:{vals2['BE name']} editing window:{int(vals2['editing window min'])}-{int(vals2['editing window max'])}bp"
                if ('BE editing window' in cfg) and ('BE type' in cfg):
                    row={'method':cfg['BE names'][0],
                    'nucleotide':cfg['BE type'][0][0],
                    'nucleotide mutation':cfg['BE type'][0][1],
                    'window start':cfg['BE editing window'][0][0],
                    'window end':cfg['BE editing window'][0][1],
                    'PAM':cfg['PAMs'][0],
                    'PAM position':cfg['PAM position'][0],
                    'guide length':cfg['guide length'][0],
                        }

                    row['distance of mutation from PAM: minimum'],row['distance of mutation from PAM: maximum']=get_distance_of_mutation_from_pam(pam_position=row['PAM position'],window_ini=row['window start'],window_end=row['window end'],guide_length=row['guide length'])
                    row['distance of codon start from PAM: minimum'],row['distance of codon start from PAM: maximum']=get_distance_of_codon_from_PAM(row['distance of mutation from PAM: minimum'],row['distance of mutation from PAM: maximum'],pam_position=row['PAM position'])
                    dbepams=dbepams.append(row,ignore_index=True)
                else:
                    logging.error('BE editing window and BE type is needed for the new BE')                
                
        dbepams['strand']='+'
        dbepams=dbes2dbes_strands(dbepams)
        dbepams=dbepams.loc[(dbepams['method'].isin(cfg['BE names']) & dbepams['PAM'].isin(cfg['PAMs'])),:]
        to_table(dbepams,cfg['dbepamsp'])            
        return cfg
    else:
        return all(opt_validity)

def validinput(cfg,din): 
    """
    Checks if input file is valid i.e. contains all the required columns.

    :param cfg: configuration dict
    :param din: dataframe containing input data 
    """
    from beditor.lib.global_vars import mutation_format2cols
    opt_validity=[]
    for col in mutation_format2cols[cfg['mutation_format']]:
        if col in din:
            opt_validity.append(True)
        else:
            
            opt_validity.append(False)
            logging.error(f"invalid column name: {col} is not in [{','.join(din.columns.tolist())}]")
    return all(opt_validity)

# related to be and pam

from beditor.lib.global_vars import cols_dpam
from beditor.lib.io_seqs import reverse_complement_multintseq,reverse_complement_multintseqreg,str2seq
from beditor.lib.io_dfs import *
from beditor.lib.io_strs import s2re
from beditor.lib.global_vars import multint2reg,multint2regcomplement,nt2complement

def dbes2dbes_strands(dBEs):
    """
    add reverse strand editing methods in the dBEs dataframe
    
    :param dBEs: pandas dataframe with BE methods
    """
    from beditor.lib.global_vars import nt2complement
    dBEs_=dBEs.copy()
    dBEs_['strand']='-'
    for col_nt in ['nucleotide','nucleotide mutation']:
        dBEs_[col_nt]=dBEs_[col_nt].apply(lambda x : nt2complement[x])
    dBEs=dBEs.append(dBEs_,sort=True)
    return dBEs

def dpam2dpam_strands(dpam,pams):
    """
    Duplicates dpam dataframe to be compatible for searching PAMs on - strand

    :param dpam: dataframe with pam information
    :param pams: pams to be used for actual designing of guides.
    """
    
    dpam=del_Unnamed(dpam)
    dpam['rPAM']=dpam.apply(lambda x : s2re(x['PAM'],multint2reg) ,axis=1)
    dpam=set_index(dpam,'PAM')
    dpam['strand']='+'
    dpamr=pd.DataFrame(columns=dpam.columns)
    dpam.loc[:,'reverse complement']=np.nan
    dpam.loc[:,'original']=np.nan
    for pam in dpam.index:
        pamr=reverse_complement_multintseq(pam,nt2complement)
        dpam.loc[pam,'reverse complement']=pamr
        dpam.loc[pam,'original']=pam
        dpamr.loc[pamr,'original']=pam
        dpam.loc[pam,'original position']=dpam.loc[pam,'PAM position']
        dpamr.loc[pamr,'original position']=dpam.loc[pam,'PAM position']
        dpamr.loc[pamr,['PAM position','guide length']]=dpam.loc[pam,['PAM position','guide length']]
        dpamr.loc[pamr,['rPAM']]=reverse_complement_multintseqreg(pam,multint2regcomplement,nt2complement)    
    dpamr['PAM position']= dpamr.apply(lambda x: 'up' if x['PAM position']=='down' else 'down',axis=1)
    dpamr['strand']='-'
    dpam_strands=dpam.append(dpamr,sort=True)
    dpam_strands.index.name='PAM'
    dpam_strands.loc[:,'is a reverse complement']=pd.isnull(dpam_strands.loc[:,'reverse complement'])
    pams_strands=pams+dpam_strands.loc[pams,'reverse complement'].dropna().tolist()
    dpam_strands=dpam_strands.loc[pams_strands,:]
    return dpam_strands

def get_be2dpam(din,test=False):
    """
    make BE to dpam mapping i.e. dict
    
    :param din: df with BE and PAM info all cols_dpam needed
    """
    be2dpam={}
    be2pam=din.loc[:,['method','PAM']].drop_duplicates().set_index('method').to_dict()['PAM']
#     if test:
#         print(be2pam)
    for be in be2pam:
        pam=be2pam[be]
        dpam=din.loc[((din['PAM']==pam) & (din['method']==be) & (din['strand']=='+')),cols_dpam]
        dpam_strands=dpam2dpam_strands(dpam,pams=[pam])
        be2dpam[be]=set_index(dpam_strands,'PAM')                
    return be2dpam

