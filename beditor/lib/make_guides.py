def make_guides(dintseqflt01,
               test=False,
               dbug=False):
    dintseqguides=dintseqflt01.copy()
    for strand in dmutagenesis.loc[:,'mutation on strand'].unique():
        for subi,sub in zip(dintseqguides.index,dintseqguides['Substrate'].tolist()):
            seq=dintseqguides.loc[subi,'30[P-codon]30']
            pos_codon=int(dintseqguides.loc[subi,'30[P-codon]30: P-codon position'])
            codon=dintseqguides.loc[subi,'P-codon']
            if strand=='- strand':
                seq=str(Seq.Seq(seq,Alphabet.generic_dna).reverse_complement())
                codon=str(Seq.Seq(codon,Alphabet.generic_dna).reverse_complement())
                pos_codon=len(seq)-(pos_codon)-3
            posGGs=[i for i in range(len(seq)) if seq.startswith('GG',i)]
            if len(posGGs)!=0:
                for posGGi,posGG in enumerate(posGGs):
                    seq_target=seq[posGG-21:posGG-1]
                    pos_codon_from_PAM=pos_codon-(posGG)+1
                    for muti in dmutagenesis[(dmutagenesis['codon']==dintseqguides.loc[subi,'P-codon']) & (dmutagenesis['mutation on strand']==strand)].index:
                        method=dmutagenesis.loc[muti,'method']
                        seq_activity=seq_target[20+dmutagenesis.loc[muti,'Position of mutation from PAM: minimum']:20+1+dmutagenesis.loc[muti,'Position of mutation from PAM: maximum']]
                        if seq_activity.count(dmutagenesis.loc[muti,'nucleotide'])==1:
                            if (pos_codon_from_PAM>=dmutagenesis.loc[muti,'Position of codon start from PAM: minimum']) and (pos_codon_from_PAM<=dmutagenesis.loc[muti,'Position of codon start from PAM: maximum']):
                                if strand=='+ strand':
                                    pos_mut_from_PAM=pos_codon_from_PAM-1+dmutagenesis.loc[muti,'position of mutation in codon']                    
                                elif strand=='- strand':
                                    pos_mut_from_PAM=pos_codon_from_PAM-1+4-dmutagenesis.loc[muti,'position of mutation in codon']                    
                                if (pos_mut_from_PAM>=dmutagenesis.loc[muti,'Position of mutation from PAM: minimum']) and (pos_mut_from_PAM<=dmutagenesis.loc[muti,'Position of mutation from PAM: maximum']):
                                    strategy='{}; {}: {} to {}; {} to {}, codon position={}; mutation position={};'.format(dmutagenesis.loc[muti,'mutation on strand'],
                                                                                 method,
                                                                                dmutagenesis.loc[muti,'codon'],
                                                                                dmutagenesis.loc[muti,'codon mutation'],
                                                                                dmutagenesis.loc[muti,'amino acid'],
                                                                                dmutagenesis.loc[muti,'amino acid mutation'],
                                                                                pos_codon_from_PAM,
                                                                                pos_mut_from_PAM
                                                                                )        
                                    codon_mut=dmutagenesis.loc[muti,'codon mutation']
                                    if strand=='- strand':
                                        codon_mut=str(Seq.Seq(codon_mut,Alphabet.generic_dna).reverse_complement())
    #                                 seq_mutated='{}{}{}'.format(seq_target[:(20+pos_codon_from_PAM)],codon_mut,seq_target[(20+pos_codon_from_PAM)+3:])
    #                                 dintseqguides.loc[subi,'PAM sequence ({0}; PAM position #{1:02d})'.format(strand,posGGi+1)]=seq[posGG-1:posGG+2]
    #                                 dintseqguides.loc[subi,'target codon ({0})'.format(strand)]=codon
    #                                 dintseqguides.loc[subi,'mutated codon ({0}; PAM position #{1:02d})'.format(strategy,posGGi+1)]=codon_mut
    #                                 dintseqguides.loc[subi,'target codon position from PAM ({0}; PAM position #{1:02d})'.format(strand,posGGi+1)]=pos_codon_from_PAM
    #                                 dintseqguides.loc[subi,'mutation position from PAM ({0}; PAM position #{1:02d})'.format(strategy,posGGi+1)]=pos_mut_from_PAM
    #                                 dintseqguides.loc[subi,'target sequence ({0}; PAM position #{1:02d})'.format(strand,posGGi+1)]=seq_target
                                    dintseqguides.loc[subi,'guide sequence ({0})'.format(strategy)]=seq_target
                                    dintseqguides.loc[subi,'guide sequence+PAM({0})'.format(strategy)]=seq_target+seq[posGG-1:posGG+2]
                                    if test:
                                        print('{}:pos_mut_from_PAM={};pos_codon_from_PAM={};seq_activity={};{}'.format(sub,pos_mut_from_PAM,pos_codon_from_PAM,seq_activity,strategy))
    #                                 if sub=="YMR072W":
    #                                 if pos_mut_from_PAM==-13:
    #                                 if strand=='- strand':
    #                                     dbug=True
    #                     `
                                    if dbug:
                                        print(posGG)
                                        print(strand)
                                        print(codon)
                                        print(seq)
                                        print(seq[30:33])
                                        print(seq_target)
                                        print(seq_activity)
                                        print(pos_codon_from_PAM)
                                        print(pos_mut_from_PAM)
        #                                     print(seq_guide)
                                        print(sdf)
                                        break
    #                     if posGG==45:
    #                         break
    #                 break
    #             break
    #         break
    #     break
    return dintseqguides