import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from Bio import SeqIO, Alphabet, Data, Seq, SeqUtils
from Bio import motifs,Seq,AlignIO

hosts={"coli":11, # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
"yeast":12,
"human":1}

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