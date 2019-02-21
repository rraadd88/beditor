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
import yaml

def collectchuckfiles(cfg,fpinchunk,force=False):
    """
    Collects minor chunk files

    :param cfg: configuration dict
    :param fpinchunk: path inside chuck's project directory
    :param force: if True overwrites the outputs 
    """
    from beditor.lib.io_dfs import fhs2data_combo_appended
    from beditor.lib.global_vars import stepi2name
    from beditor.lib.io_nums import str2num
    from os.path import basename
    from glob import glob
    dps_p=f"{cfg['prjd']}/chunks/chunk0*/{fpinchunk}"
    dps_=glob(dps_p)
#         print(dps_)
    if len(dps_)>0:
        dout=fhs2data_combo_appended(dps_,sep='\t',
                                     labels_coln='chunk#')
        return dout 
    else:
        logging.info(f'no files found: {dps_p}')
