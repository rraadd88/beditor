import numpy as np
import pandas as pd

from beditor.lib.io_nums import rescale
def get_ppamdist(guidelength,pamlength,pam_position,ppamdist_min,pmutatpam):
    """
    Get penalties set based on distances of the mismatch/es from PAM 

    :param guidelength: length of guide sequence
    :param pamlength: length of PAM sequence
    :param pam_position: PAM location 3' or 5'
    :param ppamdist_min: minimum penalty
    :param pmutatpam: penalty for mismatch at PAM
    """
    x=np.arange(pamlength)
#     ppamdist_pam=np.zeros(pamlength)
    ppamdist_pam=np.full([1, pamlength], pmutatpam)[0]
    
    guide_positions=np.arange(guidelength)
    guide_positions=rescale(guide_positions)    
    coeffs=np.array([ 0.86947513, -1.15561009, -0.0826878 ,  0.77644198])
    ppamdist=np.polyval(coeffs,guide_positions)
    ppamdist=rescale(ppamdist,mn=ppamdist_min)
    ppamdist=np.concatenate((ppamdist,ppamdist_pam))
#     ppamdist=1-ppamdist
    if pam_position=='up':
        ppamdist=ppamdist[::-1]
    elif pam_position=='down':
        ppamdist=ppamdist
    else:
        raise(ValueError('ppamdist is invalid: {ppamdist}'))
    return ppamdist

def get_beditorscore_per_alignment(NM,genic,alignment,pam_length,pam_position,
                    pentalty_genic=0.5,
                    pentalty_intergenic=0.9,
                    pentalty_dist_from_pam=0.1,
                    test=False,debug=False):
    """
    Calculates beditor score per alignment between guide and genomic DNA.

    :param NM: Hamming distance
    :param mismatches_max: Maximum mismatches allowed in alignment
    :param genic: True if guide aligns to genic regions, else (intergenic) False.
    :param alignment: Symbol '|' means a match, '.' means mismatch and ' ' means gap. e.g. |||||.||||||||||.||||.| 
    :param pentalty_genic: penalty for genic alignment
    :param pentalty_intergenic: penalty for intergenic alignment
    :param pentalty_dist_from_pam: maximum pentalty for a mismatch at PAM () 
    :returns: beditor score per alignment.
    """
    pmutatpam=1e-300
    if not pd.isnull(alignment):
        if NM!=0:            
            pentalty_region_of_alignment=pentalty_genic if genic else pentalty_intergenic 

            mutations_penalties=np.array([0 if s=='|' else 1 for s in alignment])
            dist_from_pam_penalties=get_ppamdist(guidelength=len(alignment)-pam_length,
                                                 pamlength=pam_length,
                                                 pam_position=pam_position,
                                                 ppamdist_min=pentalty_dist_from_pam,
                                                pmutatpam=pmutatpam)
            mutations_penalties_multi=mutations_penalties*dist_from_pam_penalties
            mutations_penalties_multi=mutations_penalties_multi[mutations_penalties_multi != 0]
            if len(mutations_penalties_multi)==0:
                # all the penalties are 0
                penality_cum_dist_from_pam=0
            else:
                if pmutatpam in mutations_penalties_multi:
                    # mutation at pam
                    penality_cum_dist_from_pam=0
                else:
                    penality_cum_dist_from_pam=np.prod(mutations_penalties_multi)
            if test:
                print(dist_from_pam_penalties)
                print(mutations_penalties)
                print(mutations_penalties_multi)            
                print(penality_cum_dist_from_pam)            
            if debug:
                return pentalty_region_of_alignment,penality_cum_dist_from_pam                
            return pentalty_region_of_alignment*penality_cum_dist_from_pam
        else:
            if debug:
                return np.nan,np.nan                
            return 1
    else:
        if debug:
            return np.nan,np.nan                
        return np.nan

def get_beditorscore_per_guide(guide_seq, strategy,
                               align_seqs_scores,
                              BEs,
                              penalty_activity_window=0.5,
                               test=False,
                              ):
    """
    Calculates beditor score per guide.
    
    :param guide_seq: guide seqeunce 23nts
    :param strategy: strategy string eg. ABE;+;@-14;ACT:GCT;T:A;
    :param align_seqs_scores: list of beditor scores per alignments for all the alignments between guide and genomic DNA
    :param penalty_activity_window: if editable base is not in activity window, penalty_activity_window=0.5
    :returns: beditor score per guide.
    """
    
    #create BEs and pos_muts for back-compatibility
    from os.path import dirname,realpath
    dBEs=pd.read_table(f"{dirname(realpath(__file__))}/../data/dBEs.tsv")
    dBEs=dBEs.loc[dBEs['method'].isin(BEs),:]
    pos_muts=dBEs.loc[:,['method','distance of mutation from PAM: minimum',
     'distance of mutation from PAM: maximum',
     'distance of codon start from PAM: minimum',
     'distance of codon start from PAM: maximum']].drop_duplicates().set_index('method')

    pos_mut=int(strategy.split(';')[2].replace('@',''))
    method=strategy.split(';')[0]
    penalty_activity_window=1 if (pos_muts.loc[method,'distance of mutation from PAM: minimum']<=pos_mut<=pos_muts.loc[method,'distance of mutation from PAM: maximum']) else penalty_activity_window
    penalty_align_seqs_scores=np.prod(align_seqs_scores)
    if test:
        print(list(align_seqs_scores))
        print(penalty_align_seqs_scores)
    return penalty_activity_window*penalty_align_seqs_scores

        
#Calculates the Cutting Frequency Determination score
#Requirements: 1. Pickle file with mismatch scores in working directory
#              2. Pickle file containing PAM scores in working directory 
#Input: 1. 23mer WT sgRNA sequence
#       2. 23mer Off-target sgRNA sequence
#Output: CFD score
import re

#Reverse complements a given string
def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)

#Calculates CFD score
def calc_cfd(wt,sg,pam):
    mm_scores = {'rA:dA,1': 1.0,
 'rA:dA,10': 0.882352941,
 'rA:dA,11': 0.307692308,
 'rA:dA,12': 0.333333333,
 'rA:dA,13': 0.3,
 'rA:dA,14': 0.533333333,
 'rA:dA,15': 0.2,
 'rA:dA,16': 0.0,
 'rA:dA,17': 0.133333333,
 'rA:dA,18': 0.5,
 'rA:dA,19': 0.538461538,
 'rA:dA,2': 0.727272727,
 'rA:dA,20': 0.6,
 'rA:dA,3': 0.705882353,
 'rA:dA,4': 0.636363636,
 'rA:dA,5': 0.363636364,
 'rA:dA,6': 0.7142857140000001,
 'rA:dA,7': 0.4375,
 'rA:dA,8': 0.428571429,
 'rA:dA,9': 0.6,
 'rA:dC,1': 1.0,
 'rA:dC,10': 0.5555555560000001,
 'rA:dC,11': 0.65,
 'rA:dC,12': 0.7222222220000001,
 'rA:dC,13': 0.6521739129999999,
 'rA:dC,14': 0.46666666700000003,
 'rA:dC,15': 0.65,
 'rA:dC,16': 0.192307692,
 'rA:dC,17': 0.176470588,
 'rA:dC,18': 0.4,
 'rA:dC,19': 0.375,
 'rA:dC,2': 0.8,
 'rA:dC,20': 0.764705882,
 'rA:dC,3': 0.611111111,
 'rA:dC,4': 0.625,
 'rA:dC,5': 0.72,
 'rA:dC,6': 0.7142857140000001,
 'rA:dC,7': 0.705882353,
 'rA:dC,8': 0.7333333329999999,
 'rA:dC,9': 0.666666667,
 'rA:dG,1': 0.857142857,
 'rA:dG,10': 0.333333333,
 'rA:dG,11': 0.4,
 'rA:dG,12': 0.263157895,
 'rA:dG,13': 0.21052631600000002,
 'rA:dG,14': 0.214285714,
 'rA:dG,15': 0.272727273,
 'rA:dG,16': 0.0,
 'rA:dG,17': 0.176470588,
 'rA:dG,18': 0.19047619,
 'rA:dG,19': 0.20689655199999998,
 'rA:dG,2': 0.7857142859999999,
 'rA:dG,20': 0.22727272699999998,
 'rA:dG,3': 0.428571429,
 'rA:dG,4': 0.352941176,
 'rA:dG,5': 0.5,
 'rA:dG,6': 0.454545455,
 'rA:dG,7': 0.4375,
 'rA:dG,8': 0.428571429,
 'rA:dG,9': 0.571428571,
 'rC:dA,1': 1.0,
 'rC:dA,10': 0.9411764709999999,
 'rC:dA,11': 0.307692308,
 'rC:dA,12': 0.538461538,
 'rC:dA,13': 0.7,
 'rC:dA,14': 0.7333333329999999,
 'rC:dA,15': 0.066666667,
 'rC:dA,16': 0.307692308,
 'rC:dA,17': 0.46666666700000003,
 'rC:dA,18': 0.642857143,
 'rC:dA,19': 0.46153846200000004,
 'rC:dA,2': 0.9090909090000001,
 'rC:dA,20': 0.3,
 'rC:dA,3': 0.6875,
 'rC:dA,4': 0.8,
 'rC:dA,5': 0.636363636,
 'rC:dA,6': 0.9285714290000001,
 'rC:dA,7': 0.8125,
 'rC:dA,8': 0.875,
 'rC:dA,9': 0.875,
 'rC:dC,1': 0.913043478,
 'rC:dC,10': 0.38888888899999996,
 'rC:dC,11': 0.25,
 'rC:dC,12': 0.444444444,
 'rC:dC,13': 0.13636363599999998,
 'rC:dC,14': 0.0,
 'rC:dC,15': 0.05,
 'rC:dC,16': 0.153846154,
 'rC:dC,17': 0.058823529000000006,
 'rC:dC,18': 0.133333333,
 'rC:dC,19': 0.125,
 'rC:dC,2': 0.695652174,
 'rC:dC,20': 0.058823529000000006,
 'rC:dC,3': 0.5,
 'rC:dC,4': 0.5,
 'rC:dC,5': 0.6,
 'rC:dC,6': 0.5,
 'rC:dC,7': 0.470588235,
 'rC:dC,8': 0.642857143,
 'rC:dC,9': 0.6190476189999999,
 'rC:dT,1': 1.0,
 'rC:dT,10': 0.8666666670000001,
 'rC:dT,11': 0.75,
 'rC:dT,12': 0.7142857140000001,
 'rC:dT,13': 0.384615385,
 'rC:dT,14': 0.35,
 'rC:dT,15': 0.222222222,
 'rC:dT,16': 1.0,
 'rC:dT,17': 0.46666666700000003,
 'rC:dT,18': 0.538461538,
 'rC:dT,19': 0.428571429,
 'rC:dT,2': 0.727272727,
 'rC:dT,20': 0.5,
 'rC:dT,3': 0.8666666670000001,
 'rC:dT,4': 0.842105263,
 'rC:dT,5': 0.571428571,
 'rC:dT,6': 0.9285714290000001,
 'rC:dT,7': 0.75,
 'rC:dT,8': 0.65,
 'rC:dT,9': 0.857142857,
 'rG:dA,1': 1.0,
 'rG:dA,10': 0.8125,
 'rG:dA,11': 0.384615385,
 'rG:dA,12': 0.384615385,
 'rG:dA,13': 0.3,
 'rG:dA,14': 0.26666666699999997,
 'rG:dA,15': 0.14285714300000002,
 'rG:dA,16': 0.0,
 'rG:dA,17': 0.25,
 'rG:dA,18': 0.666666667,
 'rG:dA,19': 0.666666667,
 'rG:dA,2': 0.636363636,
 'rG:dA,20': 0.7,
 'rG:dA,3': 0.5,
 'rG:dA,4': 0.363636364,
 'rG:dA,5': 0.3,
 'rG:dA,6': 0.666666667,
 'rG:dA,7': 0.571428571,
 'rG:dA,8': 0.625,
 'rG:dA,9': 0.533333333,
 'rG:dG,1': 0.7142857140000001,
 'rG:dG,10': 0.4,
 'rG:dG,11': 0.428571429,
 'rG:dG,12': 0.529411765,
 'rG:dG,13': 0.42105263200000004,
 'rG:dG,14': 0.428571429,
 'rG:dG,15': 0.272727273,
 'rG:dG,16': 0.0,
 'rG:dG,17': 0.235294118,
 'rG:dG,18': 0.47619047600000003,
 'rG:dG,19': 0.448275862,
 'rG:dG,2': 0.692307692,
 'rG:dG,20': 0.428571429,
 'rG:dG,3': 0.384615385,
 'rG:dG,4': 0.529411765,
 'rG:dG,5': 0.7857142859999999,
 'rG:dG,6': 0.681818182,
 'rG:dG,7': 0.6875,
 'rG:dG,8': 0.615384615,
 'rG:dG,9': 0.538461538,
 'rG:dT,1': 0.9,
 'rG:dT,10': 0.933333333,
 'rG:dT,11': 1.0,
 'rG:dT,12': 0.933333333,
 'rG:dT,13': 0.923076923,
 'rG:dT,14': 0.75,
 'rG:dT,15': 0.9411764709999999,
 'rG:dT,16': 1.0,
 'rG:dT,17': 0.933333333,
 'rG:dT,18': 0.692307692,
 'rG:dT,19': 0.7142857140000001,
 'rG:dT,2': 0.846153846,
 'rG:dT,20': 0.9375,
 'rG:dT,3': 0.75,
 'rG:dT,4': 0.9,
 'rG:dT,5': 0.8666666670000001,
 'rG:dT,6': 1.0,
 'rG:dT,7': 1.0,
 'rG:dT,8': 1.0,
 'rG:dT,9': 0.642857143,
 'rU:dC,1': 0.956521739,
 'rU:dC,10': 0.5,
 'rU:dC,11': 0.4,
 'rU:dC,12': 0.5,
 'rU:dC,13': 0.260869565,
 'rU:dC,14': 0.0,
 'rU:dC,15': 0.05,
 'rU:dC,16': 0.346153846,
 'rU:dC,17': 0.117647059,
 'rU:dC,18': 0.333333333,
 'rU:dC,19': 0.25,
 'rU:dC,2': 0.84,
 'rU:dC,20': 0.176470588,
 'rU:dC,3': 0.5,
 'rU:dC,4': 0.625,
 'rU:dC,5': 0.64,
 'rU:dC,6': 0.571428571,
 'rU:dC,7': 0.588235294,
 'rU:dC,8': 0.7333333329999999,
 'rU:dC,9': 0.6190476189999999,
 'rU:dG,1': 0.857142857,
 'rU:dG,10': 0.533333333,
 'rU:dG,11': 0.666666667,
 'rU:dG,12': 0.947368421,
 'rU:dG,13': 0.7894736840000001,
 'rU:dG,14': 0.28571428600000004,
 'rU:dG,15': 0.272727273,
 'rU:dG,16': 0.666666667,
 'rU:dG,17': 0.705882353,
 'rU:dG,18': 0.428571429,
 'rU:dG,19': 0.275862069,
 'rU:dG,2': 0.857142857,
 'rU:dG,20': 0.090909091,
 'rU:dG,3': 0.428571429,
 'rU:dG,4': 0.647058824,
 'rU:dG,5': 1.0,
 'rU:dG,6': 0.9090909090000001,
 'rU:dG,7': 0.6875,
 'rU:dG,8': 1.0,
 'rU:dG,9': 0.923076923,
 'rU:dT,1': 1.0,
 'rU:dT,10': 0.857142857,
 'rU:dT,11': 0.75,
 'rU:dT,12': 0.8,
 'rU:dT,13': 0.692307692,
 'rU:dT,14': 0.6190476189999999,
 'rU:dT,15': 0.578947368,
 'rU:dT,16': 0.9090909090000001,
 'rU:dT,17': 0.533333333,
 'rU:dT,18': 0.666666667,
 'rU:dT,19': 0.28571428600000004,
 'rU:dT,2': 0.846153846,
 'rU:dT,20': 0.5625,
 'rU:dT,3': 0.7142857140000001,
 'rU:dT,4': 0.47619047600000003,
 'rU:dT,5': 0.5,
 'rU:dT,6': 0.8666666670000001,
 'rU:dT,7': 0.875,
 'rU:dT,8': 0.8,
 'rU:dT,9': 0.9285714290000001}
    pam_scores = {'AA': 0.0,
 'AC': 0.0,
 'AG': 0.25925925899999996,
 'AT': 0.0,
 'CA': 0.0,
 'CC': 0.0,
 'CG': 0.107142857,
 'CT': 0.0,
 'GA': 0.06944444400000001,
 'GC': 0.022222222000000003,
 'GG': 1.0,
 'GT': 0.016129031999999998,
 'TA': 0.0,
 'TC': 0.0,
 'TG': 0.038961038999999996,
 'TT': 0.0}
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

# wt = 'ATTCAATCCTATGCTGACTTGGG'
# off = 'ATTCAATCCTATGATGACTTGGG'
def get_cfdscore(wt,off):
    m_wt = re.search('[^ATCG]',wt)
    m_off = re.search('[^ATCG]',off)    
    if (m_wt is None) and (m_off is None):
        if len(wt)== 23 and len(off)== 23:
            pam = off[-2:]
            sg = off[:-3]
            cfd_score = calc_cfd(wt,sg,pam)
    #         print("CFD score: "+str(cfd_score))
            return cfd_score
    else:
        np.nan