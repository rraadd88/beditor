#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
from os import makedirs
from glob import glob
import pandas as pd
import subprocess
import logging
from beditor.lib.io_sys import runbashcmd    
    
def get_genomes(cfg):
    """
    Installs genomes
    
    :param cfg: configuration dict
    """
    
    runbashcmd(f"pyensembl install --reference-name {cfg['genomeassembly']} --release {cfg['genomerelease']} --species {cfg['host']}")

    import pyensembl
    ensembl = pyensembl.EnsemblRelease(species=pyensembl.species.Species.register(
                                        latin_name=cfg['host'],
                                        synonyms=[cfg['host']],
                                        reference_assemblies={
                                            cfg['genomeassembly']: (cfg['genomerelease'], cfg['genomerelease']),
                                        }),release=cfg['genomerelease'])
    contig_mito=['MTDNA','MITO','MT']
    contigs=[c for c in ensembl.contigs() if ((not '.' in c) and (len(c)<5) and (c not in contig_mito))]    
    if len(contigs)==0:
        contigs=[c for c in ensembl.contigs()]
        # logging.error([c for c in ensembl.contigs()])
        # logging.error('no contigs identified by pyensembl; aborting')
        # sys.exit(0)
    logging.info(f"{len(contigs)} contigs/chromosomes in the genome")
    logging.info(contigs)
    # raw genome next
    if 'human' in cfg['host'].lower():
        cfg['host']='homo_sapiens'
    if 'yeast' in cfg['host'].lower():
        cfg['host']='saccharomyces_cerevisiae'
    host_="_".join(s for s in cfg['host'].split('_')).capitalize()
    ensembl_fastad='pub/release-{}/fasta/{}/dna/'.format(cfg['genomerelease'],cfg['host'])
    genome_fastad='{}/{}'.format(dirname(realpath(__file__)),ensembl_fastad)
    cfg['genomep']='{}/genome.fa'.format(genome_fastad)
    if not exists(cfg['genomep']):
        logging.error('not found: {}'.format(cfg['genomep']))
        if not '/test_beditor/' in cfg['cfgp']:
            ifdlref = input("Download genome at {}?[Y/n]: ".format(genome_fastad))
        else:
            ifdlref='Y'
        if ifdlref=='Y':
        # #FIXME download contigs and cat and get index, sizes
            for contig in contigs:
                if 'GRCh37' in cfg['genomeassembly']:
                    #Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz
                    fn=f"{cfg['host'].capitalize()}.{cfg['genomeassembly']}.{cfg['genomerelease']}.dna_sm.chromosome.{contig}.fa.gz"
                else:
                    fn=f"{cfg['host'].capitalize()}.{cfg['genomeassembly']}.dna_sm.chromosome.{contig}.fa.gz"
                fp='{}/{}'.format(ensembl_fastad,fn)
                if not exists(fp):
                    cmd=f"wget -q -x -nH ftp://ftp.ensembl.org/{fp} -P {dirname(realpath(__file__))}"
                    try:
                        runbashcmd(cmd,test=cfg['test'])
                    except:
                        fn=f"{cfg['host'].capitalize()}.{cfg['genomeassembly']}.dna_sm.toplevel.fa.gz"
                        fp='{}/{}'.format(ensembl_fastad,fn)                        
                        if not exists(fp):
                            cmd=f"wget -q -x -nH ftp://ftp.ensembl.org/{fp} -P {dirname(realpath(__file__))}"
                            # print(cmd)
                            runbashcmd(cmd,test=cfg['test'])
                            break
            #                 break
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

    ensembl_gff3d='pub/release-{}/gff3/{}/'.format(cfg['genomerelease'],cfg['host'])    
    genome_gff3d='{}/{}'.format(dirname(realpath(__file__)),ensembl_gff3d)
    cfg['genomegffp']='{}/genome.gff3'.format(genome_gff3d)
    if not exists(cfg['genomegffp']):
        logging.error('not found: {}'.format(cfg['genomegffp']))
        if not '/test_beditor/' in cfg['cfgp']:
            ifdlref = input("Download genome annotations at {}?[Y/n]: ".format(genome_gff3d))
        else:
            ifdlref = 'Y'
        if ifdlref=='Y':
        # #FIXME download contigs and cat and get index, sizes
            fn='{}.{}.{}.gff3.gz'.format(cfg['host'].capitalize(),cfg['genomeassembly'],cfg['genomerelease'])
            fp='{}/{}'.format(ensembl_gff3d,fn)
            if not exists(fp):
                cmd='wget -x -nH ftp://ftp.ensembl.org/{} -P {}'.format(fp,dirname(realpath(__file__)))
                runbashcmd(cmd,test=cfg['test'])
            # move to genome.gff3
                cmd='cp {}/{} {}'.format(genome_gff3d,fn,cfg['genomegffp'])
                runbashcmd(cmd,test=cfg['test'])

        else:
            logging.error('abort')
            sys.exit(1)
    logging.info('genomes are installed!')
    return cfg

# related to be and pam

from beditor.lib.global_vars import cols_dpam
from beditor.lib.io_seqs import reverse_complement_multintseq,reverse_complement_multintseqreg,str2seq
from beditor.lib.io_dfs import *
from beditor.lib.io_strs import s2re
from beditor.lib.global_vars import multint2reg,multint2regcomplement,nt2complement

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
    if test:
        print(be2pam)
    for be in be2pam:
        pam=be2pam[be]
        dpam=din.loc[((din['PAM']==pam) & (din['method']==be)),cols_dpam]
        dpam_strands=dpam2dpam_strands(dpam,pams=[pam])
        be2dpam[be]=set_index(dpam_strands,'PAM')                
    return be2dpam