import pandas as pd
import numpy as np
from os.path import exists, dirname
import subprocess

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

import logging

# input 
# mutation_format and must-have columns
mutation_format2cols={
    'aminoacid': ['transcript: id','aminoacid: position','amino acid mutation'],
    'nucleotide':['genome coordinate','nucleotide mutation'],
}

# configuration 
# allowed options
cfgoption2allowed={
'mutations':['mutations','substitutions','mimetic',None],
'mutation_format':['aminoacid','nucleotide'],
}
cfgoption2reguired={
'mutations':{'substitutions':'dsubmap_preferred_path'},
}

#output 
# cols
stepi2colsoutput={
    0: mutation_format2cols['nucleotide']+mutation_format2cols['aminoacid'],
    1: mutation_format2cols['nucleotide']+['nucleotide wild-type','codon: wild-type']+['transcript: id','aminoacid: position','aminoacid: wild-type'],
    3: ['nucleotide mutation', 'nucleotide wild-type']+['transcript: id','aminoacid: position','aminoacid: wild-type',
        'amino acid mutation','guide: id','guide+PAM sequence'],
    4: ['guide: id','beditor score','alternate alignments count','CFD score'],
   }


flankaac=7
flankntc=22
flankaal=15
flankntl=45

# steps
# 0:'input',
stepi2name= {
 1: 'sequences',
 2: 'mutagenesis',
 3: 'guides',
 4: 'offtargets'}

stepi2cols_nucleotide={
1:['genome coordinate','nucleotide mutation','id','transcript: sequence','nucleotide wild-type','codon: wild-type','codon: mutation','transcript: id']}

stepi2cols={
1: ['aminoacid: position', 'gene: id', 'gene: name', 'protein: id', 'transcript: id', 'transcript: sequence', 'aminoacid: wild-type', 'codon: wild-type', 'contig', 'strand', 'start', 'end', 'codon start', 'codon end'],
2:  ['amino acid',
 'amino acid mutation',
 'codon',
 'codon mutation',
 'codon mutation usage Fraction',
 'codon mutation usage Frequency',
 'method',
 'mutation on strand',
 'nucleotide',
 'nucleotide mutation',
 'nucleotide mutation: count',
 'position of mutation in codon'],    
3: ['PAM','Description','guide length','original','original position','position','rPAM','reverse complement','strand','is a reverse complement','guide+PAM sequence','guide sequence','PAM sequence','position of PAM ini','position of PAM end','distance of codon from PAM','codon','guide sequence length','index','transcript: id','aminoacid: position','aminoacid: wild-type','codon: wild-type','id','amino acid mutation','method','amino acid','position of mutation in codon','codon mutation','nucleotide','nucleotide mutation','mutation on strand','codon mutation usage Fraction','codon mutation usage Frequency','nucleotide mutation: count','distance of mutation from PAM: minimum','distance of mutation from PAM: maximum','distance of codon start from PAM: minimum','distance of codon start from PAM: maximum','activity sequence','distance of mutation from PAM','strategy','guide: id','guide+PAM length'],
4: ['guide: id','id','guide+PAM sequence','gene names','gene ids','transcript ids','types','protein ids','exon ids','beditor score','CFD score','beditor score (log10)','alternate alignments count'],    
5:['transcript: id',
 'aminoacid: position',
 'amino acid mutation',
 'aminoacid: wild-type',
 'guide: id',
 'guide+PAM sequence',
 'beditor score',
 'alternate alignments count',
 'CFD score'],
}

def saveemptytable(cfg,doutp):
    from beditor.lib.global_vars import stepi2cols
    dout=pd.DataFrame(columns=stepi2cols[cfg['step']])
    logging.warning(f"saved enmpty table {doutp}")
    dout.to_csv(doutp,sep='\t')


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


# common 
aminoacids=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","*"] #for indexing
aminoacids_3letter=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
codons=["TTT",    "TTC",    "TTA",  "TTG",  "TCT",  "TCC",  "TCA",  "TCG",  "TAT",  "TAC",  "TAA",  "TAG",  "TGT",  "TGC",  "TGA",  "TGG",  "CTT",  "CTC",  "CTA",  "CTG",  "CCT",  "CCC",  "CCA",  "CCG",  "CAT",  "CAC",  "CAA",  "CAG",  "CGT",  "CGC",  "CGA",  "CGG",  "ATT",  "ATC",  "ATA",  "ATG",  "ACT",  "ACC",  "ACA",  "ACG",  "AAT",  "AAC",  "AAA",  "AAG",  "AGT",  "AGC",  "AGA",  "AGG",  "GTT",  "GTC",  "GTA",  "GTG",  "GCT",  "GCC",  "GCA",  "GCG",  "GAT",  "GAC",  "GAA",  "GAG",  "GGT",  "GGC",  "GGA",  "GGG"]
