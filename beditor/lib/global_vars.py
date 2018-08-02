import pandas as pd
import numpy as np
from os.path import exists, dirname
import subprocess

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

# steps
# 0:'input',
stepi2name= {
 1: 'sequences',
 2: 'mutagenesis',
 3: 'guides',
 4: 'offtargets'}

# # regex
multint2reg={'R':'[AG]',
'Y':'[CT]',
'S':'[GC]',
'W':'[AT]',
'K':'[GT]',
'M':'[AC]',
'B':'[CGT]',
'D':'[AGT]',
'H':'[ACT]',
'V':'[ACG]',
'N':'[ATGC]',}
multint2regcomplement={'R':'[TC]',
'Y':'[GA]',
'S':'[GC]',
'W':'[AT]',
'K':'[CA]',
'M':'[TG]',
'B':'[^A]',
'D':'[^C]',
'H':'[^G]',
'V':'[^T]',
'N':'[ATGC]',}

def get_nt2complement(): 
    nt2complement={'A':'T',
                  'G':'C',
                  'N':'N',
                  'R':'Y',
                  'S':'W',
                  'K':'M',
                   'B':'b',
                   'D':'d',
                   'H':'h',
                   'V':'v',
                   'N':'N',
                  }
    nt2complement.update(dict(zip(nt2complement.values(),nt2complement.keys())))
    return nt2complement
nt2complement=get_nt2complement()

# EXT 
import beditor 
dirs2ps={'pyp':str(subprocess.check_output('which python3'.split(' '))).replace("b'",'').replace("\\n'",''),
#'binDir':dirname(str(subprocess.check_output('which bwa'.split(' '))).replace("b'",'').replace("\\n'",'')),
'binDir':dirname(beditor.__file__)+'/bin', 
'scriptDir':dirname(beditor.__file__)+'/bin',
}

#use good old bash programs for speed
# bed_colns = ['chromosome','start','end','id','NM','strand']
gff_colns = ['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
bed_colns = ['chromosome','start','end','id','NM','strand']

## ref coordinate system
# UCSC : 0 based 
# ENSEMBL : 1 based
