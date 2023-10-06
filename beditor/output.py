#!usr/bin/python
"""Output management"""

import logging
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
    from glob import glob
    dps_p=f"{cfg['prjd']}/chunks/chunk0*/{fpinchunk}"
    dps_=glob(dps_p)
#         logging.info(dps_)
    if len(dps_)>0:
        dout=fhs2data_combo_appended(dps_,sep='\t',
                                     labels_coln='chunk#')
        return dout 
    else:
        logging.info(f'no files found: {dps_p}')
