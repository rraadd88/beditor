import pandas as pd
import numpy as np
from os.path import exists, dirname
import subprocess

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

hosts={"coli":11, # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
"saccharomyces_cerevisiae":12,
"homo_sapiens":1}

BEs={'Target-AID on + strand':['C',['T','G']],
'ABE on + strand':['A',['G']],
'Target-AID on - strand':['G',['A','C']],
'ABE on - strand':['T',['C']],
}


pos_muts=pd.DataFrame(index=['ABE','Target-AID'],columns=['Position of mutation from PAM: minimum','Position of mutation from PAM: maximum'])
pos_muts.index.name='method'
pos_muts.loc['ABE','Position of mutation from PAM: minimum']=-17
pos_muts.loc['ABE','Position of mutation from PAM: maximum']=-13
pos_muts.loc['Target-AID','Position of mutation from PAM: minimum']=-19
pos_muts.loc['Target-AID','Position of mutation from PAM: maximum']=-17

pos_muts['Position of codon start from PAM: minimum']=pos_muts['Position of mutation from PAM: minimum']-2
pos_muts['Position of codon start from PAM: maximum']=pos_muts['Position of mutation from PAM: maximum']

# regex

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
    return nt2complement.update(dict(zip(nt2complement.values(),nt2complement.keys())))
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
host2contigs={'homo_sapiens':['1','2','3','4','5','6','7','8','9','10','11','12','13','14',
                              '15','16','17','18','19','20','21','22','X','Y','MT'],
'saccharomyces_cerevisiae':['I',
'II',
'III',
'IV',
'IX',
'Mito',
'V',
'VI',
'VII',
'VIII',
'X',
'XI',
'XII',
'XIII',
'XIV',
'XV',
'XVI']}
