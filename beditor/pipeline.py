#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

import sys
from os.path import exists,splitext,dirname
from os import makedirs
import argparse
import pkg_resources
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_dh+'.log'

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
    parser.add_argument("--step", help="1: get seqeucnces,\n2: get possible strategies,\n3: make guides", dest="step", 
                        type=float,action="store", choices=[1,2,3],default=None)  
    parser.add_argument("--test", help="Debug mode on", dest="test", 
                        action='store_true', default=False)    
    parser.add_argument("--force", help="Debug mode on", dest="force", 
                        action='store_true', default=False)    
    parser.add_argument('-v','--version', action='version',version=version_info)
#    parser.add_argument('-h', '--help', action='help', #default=argparse.SUPPRESS,
#                    help='Show this help message and exit. \n Version info: %s' % version_info)
    args = parser.parse_args()
    logging.info("start")
    pipeline(args.cfg,step==args.step,
        test=args.test,force=args.force)

def pipeline(cfgp,test=False,force=False):        
    from beditor.lib.get_seq import din2dseq
    from beditor.lib.get_mutations import dseq2dmutagenesis 
    from beditor.lib.make_guides import dseq2dguides

    with open(cfgp, 'r') as f:
        import json
        cfg = json.load(f)   
    cfg['prj']=basename(cfgp)
    cfg['prjd']=dirname(cfgp)+cfg['prj']
    cfg['datad']=cfg['prjd']+'/data'
    cfg['plotd']=cfg['prjd']+'/plot'
    if not exists(prjd) or force:
        makedirs(cfg['prjd'])
        makedirs(cfg['datad'])
        makedirs(cfg['plotd'])
    cfg['test']=test
    cfg['force']=force
    cfg['cfgp']=cfgp

    if exists(cfg):
        if step==1 or step==None:
            din2dseq(cfg)
        if step==2 or step==None:
            dseq2dmutagenesis(cfg)
        if step==3 or step==None:
            dseq2dguides(cfg)
        if step==None:
            logging.info("Location of output data: {}".format(cfg['datad']))
            logging.info("Location of output plot: {}".format(cfg['datad']))
    else:
        logging.error('configuration file (cfg): {} not found'.format(cfg))                  
    logging.shutdown()

if __name__ == '__main__':
    main()
