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

stepi2cols={
3: ['PAM','Description','guide length','original','original position','position','rPAM','reverse complement','strand','is a reverse complement','guide+PAM sequence','guide sequence','PAM sequence','position of PAM ini','position of PAM end','distance of codon from PAM','codon','guide sequence length','index','transcript: id','aminoacid: position','aminoacid: wild-type','codon: wild-type','id','amino acid mutation','Unnamed: 0','method','amino acid','position of mutation in codon','codon mutation','nucleotide','nucleotide mutation','mutation on strand','codon mutation usage Fraction','codon mutation usage Frequency','nucleotide mutation: count','Position of mutation from PAM: minimum','Position of mutation from PAM: maximum','Position of codon start from PAM: minimum','Position of codon start from PAM: maximum','distance of mutation from PAM: minimum','distance of mutation from PAM: maximum','distance of codon start from PAM: minimum','distance of codon start from PAM: maximum','activity sequence','distance of mutation in codon from PAM','strategy','guide: id','guide+PAM length'],
4: ['guide: id','id','guide+PAM sequence','gene names','gene ids','transcript ids','types','protein ids','exon ids','beditor score','CFD score','beditor score (log10)','alternate alignments count'],    
}

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
