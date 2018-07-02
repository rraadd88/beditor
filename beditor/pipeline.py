#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname,splitext,basename,realpath
from os import makedirs
import argparse
import pkg_resources
from beditor.lib.io_strs import get_logger
logging=get_logger()

# GET INPTS    
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
    parser.add_argument("--step", help="1: get seqeucnces,\n2: get possible strategies,\n3: make guides,\n 4: check offtargets", dest="step", 
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

def pipeline(cfgp,step=None,test=False,force=False):        
    from beditor.lib.get_seq import din2dseq
    from beditor.lib.get_mutations import dseq2dmutagenesis 
    from beditor.lib.make_guides import dseq2dguides
    from beditor.lib.get_specificity import dguides2offtargets
    from glob import glob

    import yaml
    cfg=yaml.load(open(cfgp, 'r'))
    # basics
    cfg['prj']=splitext(basename(cfgp))[0]
    if dirname(cfgp)!='':
        cfg['prjd']=dirname(cfgp)+'/'+cfg['prj']
    else:
        cfg['prjd']=cfg['prj']
    cfg['test']=test
    cfg['force']=force
    cfg['cfgp']=cfgp

# refs
    if 'human' in cfg['host'].lower():
        cfg['host']='homo_sapiens'
    if 'yeast' in cfg['host'].lower():
        cfg['host']='saccharomyces_cerevisiae'
    host_="_".join(s for s in cfg['host'].split('_')).capitalize()
    genomed_ensembl_fasta='pub/release-{}/fasta/{}/dna/'.format(cfg['genomerelease'],cfg['host'])
    genomed='{}/{}'.format(dirname(realpath(__file__)),genomed_ensembl_fasta)
    cfg['genomep']='{}/genome.fa'.format(genomed)
    if not exists(cfg['genomep']):
        logging.error('not found: {}'.format(cfg['genomep']))
        ifdlref = input("\nDownload genome at {}?[Y/n]: ".format(genomed))
        if ifdlref=='Y':
        # #FIXME download contigs and cat and get index, sizes
            from beditor.lib.global_vars import host2contigs
            from beditor.lib.io_sys import runbashcmd
            for contig in host2contigs[cfg['host']]:
                fn='{}.{}.dna_sm.chromosome.{}.fa.gz'.format(cfg['host'].capitalize(),cfg['genomeassembly'],contig)
                fp='{}/{}'.format(genomed_ensembl_fasta,fn)
                if not exists(fp) or cfg['force']:
                    cmd='wget -x -nH ftp://ftp.ensembl.org/{} -P {}'.format(fp,dirname(realpath(__file__)))
                    runbashcmd(cmd,test=cfg['test'])
#                 break
            # make the fa ready
            if not exists(cfg['genomep']) or cfg['force']:
                cmd='gunzip {}*.fa.gz;cat {}/*.fa > {}/genome.fa;'.format(genomed,genomed,genomed)
                runbashcmd(cmd,test=cfg['test'])
            if not exists(cfg['genomep']+'.bwt') or cfg['force']:
                cmd='bwa index {}'.format(cfg['genomep'])
                runbashcmd(cmd,test=cfg['test'])
            if not exists(cfg['genomep']+'.fai') or cfg['force']:
                cmd='samtools faidx {}'.format(cfg['genomep'])
                runbashcmd(cmd,test=cfg['test'])
            if not exists(cfg['genomep']+'.sizes') or cfg['force']:
                cmd='cut -f1,2 {}.fai > {}.sizes'.format(cfg['genomep'],cfg['genomep'])            
                runbashcmd(cmd,test=cfg['test'])
        else:
            sys.exit(1)
    
    genomeannotd='{}/lib/pub/release-{}/gff3/'.format(dirname(realpath(__file__)),cfg['genomerelease'])
    cfg['genomegffp']='{}/{}/{}.{}.{}.gff3.gz'.format(genomeannotd,cfg['host'],host_,cfg['genomeassembly'],cfg['genomerelease'])

   
    #datads
    cfg[0]=cfg['prjd']+'/00_input/'
    cfg[1]=cfg['prjd']+'/01_sequences/'
    cfg[2]=cfg['prjd']+'/02_mutagenesis/'
    cfg[3]=cfg['prjd']+'/03_guides/'
    cfg[4]=cfg['prjd']+'/04_offtargets/'

    if not exists(cfg['prjd']):
        makedirs(cfg['prjd'])
    for i in range(0,4+1,1):
        if not exists(cfg[i]):
            makedirs(cfg[i])
    if step==1 or step==None:
        cfg['step']=1
        din2dseq(cfg)
    if step==2 or step==None:
        cfg['step']=2
        dseq2dmutagenesis(cfg)
    if step==3 or step==None:
        cfg['step']=3
        dseq2dguides(cfg)
    if step==4 or step==None:
        cfg['step']=4
        dguides2offtargets(cfg)

    logging.info("Location of output data: {}".format(cfg['datad']))
    logging.info("Location of output plot: {}".format(cfg['plotd']))

    logging.shutdown()

if __name__ == '__main__':
    main()
