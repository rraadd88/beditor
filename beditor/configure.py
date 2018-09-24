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
    
# GET INPTS    
def get_deps(cfg):
    """
    Installs dependencies of `beditor`

    :param cfg: configuration dict
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
    ddeps.loc[dep,'install']='cd {};./configure --disable-lzma;make;'.format(dirname(ddeps.loc[dep,'executable']))

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
    contigs=[c for c in ensembl.contigs() if ((not '.' in c) and (c not in contig_mito))]    

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
                fn='{}.{}.dna_sm.chromosome.{}.fa.gz'.format(cfg['host'].capitalize(),cfg['genomeassembly'],contig)
                fp='{}/{}'.format(ensembl_fastad,fn)
                if not exists(fp):
                    cmd='wget -q -x -nH ftp://ftp.ensembl.org/{} -P {}'.format(fp,dirname(realpath(__file__)))
                    runbashcmd(cmd,test=cfg['test'])
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