#!usr/bin/python
"""Designing the sgRNAs"""

import logging

import pandas as pd
import regex as re

from roux.lib.set import get_alt
import roux.lib.df as rd #noqa

from beditor.lib.utils import get_pos, parse_locus, to_locus

def get_guide_pam(
    match: str,
    pam_stream: str,
    guidel: int,
    seq: str,
    pos_codon: int = None,
):
    if pam_stream == "down":
        # IIIIIIINGG
        # 0123456789
        seq_guidepam = seq[match.span()[0] - guidel : match.span()[1]]
        seq_guide = seq[match.span()[0] - guidel : match.span()[0]]
        seq_pam = seq[match.span()[0] : match.span()[1]]
        if pos_codon is not None:
            dist_codon = pos_codon - match.span()[0]
    #             if test:
    #                 print(match.span()[0]-pos_codon)
    elif pam_stream == "up":
        # TTNIIIIIII
        # 0123456789
        seq_guidepam = seq[match.span()[0] : match.span()[1] + guidel]
        seq_guide = seq[match.span()[1] : match.span()[1] + guidel]
        seq_pam = seq[match.span()[0] : match.span()[1]]
        if pos_codon is not None:
            dist_codon = pos_codon - match.span()[1]
    if seq_pam != match.group():
        logging.error(
            f"indexing is wrong:seq_guidepam: {seq_guidepam}, seq_guide: {seq_guide}, seq_pam: {seq_pam},match.group(): {match.group()}"
        )
    return (
        seq_guidepam,
        seq_guide,
        seq_pam,
        abs(dist_codon) if pos_codon is not None else None,
    )


def get_pam_searches(
    dpam: pd.DataFrame,
    seq: str,
    pos_codon: int,
    # test=False
) -> pd.DataFrame:
    """
    Search PAM occurance

    :param dpam: dataframe with PAM sequences
    :param seq: target sequence
    :param pos_codon: reading frame
    :param test: debug mode on
    :returns dpam_searches: dataframe with positions of pams
    """
    res = []
    for _, x in dpam.iterrows():
        matchesiter = re.finditer(x["rPAM"], seq, overlapped=True)
        for match in matchesiter:
            d1 = {}
            d1["position of PAM ini"], d1["position of PAM end"] = match.span()
            d1["position of PAM end"] = d1["position of PAM end"] - 1
            (
                d1["guide+PAM sequence"],
                d1["guide sequence"],
                d1["PAM sequence"],
                d1["distance of codon from PAM"],
            ) = get_guide_pam(
                match,
                x["PAM position"],
                x["guide length"],
                seq=seq,
                pos_codon=pos_codon,
            )
            d1["strand"] = x["strand"]
            res.append(d1)
    if len(res) == 0:
        return
    df1 = pd.DataFrame(res)
    if pos_codon is not None:
        df1["codon: from pam search"] = seq[pos_codon : pos_codon + 3]
    return df1.assign(
        **{
            "guide sequence": lambda df: df["guide sequence"].fillna(""),
            "guide sequence length": lambda df: df.apply(
                lambda x: len(x["guide sequence"]), axis=1
            ),
        }
    )


def get_guides(
    data: pd.DataFrame,
    dpam: pd.DataFrame,
    guide_len: int,
    base_fraction_max: float = 0.8,
) -> pd.DataFrame:
    """Get guides

    Args:
        data (pd.DataFrame): input table
        dpam (pd.DataFrame): table with PAM info
        guide_len (int): guide length
        base_fraction_max (float, optional): base fraction max. Defaults to 0.8.

    Returns:
        pd.DataFrame: output table
    """
    from beditor.lib.utils import get_orep, str2seq

    return (
        data.loc[:, ["sequence flanking"]]
        .drop_duplicates()
        .groupby("sequence flanking")
        .apply(
            lambda df: get_pam_searches(
                dpam=dpam,
                seq=df.name,
                pos_codon=None,
            )
        )
        .reset_index(0)
        # filter by length
        .log.query(expr=f"`guide sequence length`=={guide_len}")
        # filter by overrepresentation
        .assign(
            **{
                "base fraction max": lambda df: df.apply(
                    lambda x: get_orep(x["guide sequence"]) / guide_len, axis=1
                ),
            },
        )
        .log.query(expr=f"`base fraction max` < {base_fraction_max}")
        .drop(["guide sequence length"], axis=1)
        # .drop(['original','original position'],axis=1)
        .dropna(axis=1, how="all")
        .rename(
            columns={
                "position of PAM ini": "start PAM in flanking",
                "position of PAM end": "end PAM in flanking",
                "strand": "guide strand",
                "guide sequence": "guide sequence raw",
                "guide+PAM sequence": "guide+PAM sequence raw",  # not rev com
                "PAM sequence": "PAM sequence raw",
            }
        )
        # .drop(['guide sequence'],axis=1) #because needed for the alignment
        .assign(
            **{
                "guide sequence": lambda df: df.apply(
                    lambda x: str(str2seq(x["guide sequence raw"]).reverse_complement())
                    if x["guide strand"] == "-"
                    else x["guide sequence raw"],
                    axis=1,
                ),
                "PAM sequence": lambda df: df.apply(
                    lambda x: str(str2seq(x["PAM sequence raw"]).reverse_complement())
                    if x["guide strand"] == "-"
                    else x["PAM sequence raw"],
                    axis=1,
                ),
                "guide+PAM sequence": lambda df: df.apply(
                    lambda x: str(
                        str2seq(x["guide+PAM sequence raw"]).reverse_complement()
                    )
                    if x["guide strand"] == "-"
                    else x["guide+PAM sequence raw"],
                    axis=1,
                ),
            }
        )
    )


def to_locusby_pam(
    chrom: str,
    pam_start: int,
    pam_end: int,
    pam_position: str,
    strand: str,
    length: int,  ## for farthest from pam
    start_off: int = 0,  ## for nearest from pam
) -> str:
    """To locus by PAM from PAM coords.

    Args:
        chrom (str): chrom
        pam_start (int): PAM start
        pam_end (int): PAM end
        pam_position (str): PAM position
        strand (str): strand
        length (int): length

    Returns:
        str: locus
    """
    # print(pam_position,strand)
    if pam_position == "up":  #
        # PAM---
        # if strand=='+':
        # PAM---
        start, end = (
            pam_end + start_off,
            pam_end + length + (1 if start_off != 0 else 0),
        )
        # elif strand=='-':
        #     # ---PAM
        #     start,end=(pam_start-length)-1,pam_start-1
    elif pam_position == "down":  #
        # ---PAM
        # if strand=='+':
        # ---PAM
        start, end = (
            (pam_start - length) - 1 - (1 if start_off != 0 else 0),
            pam_start - 1 - start_off,
        )
        # elif strand=='-':
        # PAM---
        # start,end=pam_end,pam_end+length
    # return f"{chrom}:{int(start)}-{int(end)}"
    return to_locus(
        chrom,
        int(start),
        int(end),
        strand,
    )


def to_pam_coord(
    startf: int,
    endf: int,
    startp: int,
    endp: int,
    strand: str,
) -> tuple:
    """To PAM coordinates

    Args:
        startf (int): start flank start
        endf (int): start flank end
        startp (int): start PAM start
        endp (int): start PAM end
        strand (str): strand

    Returns:
        tuple: start,end
    """
    # print(startf,endf,startp,endp,strand)
    # if strand=='+':
    #     start=startf+startp
    #     end=startf+endp
    # elif strand=='-':
    r = range(startf, endf)[startp:endp]
    start = endf - endp
    end = endf - startp
    start, end = r[0], r[-1] + (1 if strand == "+" else 0)
    return start, end


def get_distances(
    df2: pd.DataFrame,
    df3: pd.DataFrame,
    cfg_method: dict,
) -> pd.DataFrame:
    """Get distances

    Args:
        df2 (pd.DataFrame): input table #1
        df3 (pd.DataFrame): input table #2
        cfg_method (dict): config for the method

    Returns:
        pd.DataFrame: output table
    """
    from roux.lib.set import get_alt

    ## note: start < end, even for the negative strand
    return (
        df2.log.merge(
            right=df3,
            on="sequence flanking",
            how="inner",
        )
        .assign(
            **{
                ## convert pam positions to genome coords
                "start,end PAM": lambda df: df.apply(
                    lambda x: to_pam_coord(
                        startf=x["start flanking"] - 1,
                        endf=x["end flanking"],
                        startp=x["start PAM in flanking"] + 1,
                        endp=x["end PAM in flanking"] + 1,
                        strand=x[
                            "guide strand"
                        ],  # get_alt(['+','-'],x['strand']) if x['is a reverse complement'] else x['strand'],
                    ),
                    axis=1,
                ),
                "start PAM": lambda df: df["start,end PAM"].apply(lambda x: x[0]),
                "end PAM": lambda df: df["start,end PAM"].apply(lambda x: x[1]),
                "guide locus": lambda df: df.apply(
                    lambda x: to_locusby_pam(
                        chrom=x["chrom"],
                        pam_start=x["start PAM"],
                        pam_end=x["end PAM"],
                        pam_position=get_alt(["up", "down"], cfg_method["PAM position"])
                        if x["guide strand"] == "-"
                        else cfg_method["PAM position"],
                        strand=x["guide strand"],
                        length=cfg_method["guide length"],
                    ),
                    axis=1,
                ),
            },
        )
        .drop(
            [
                "start,end PAM",
                "sequence target",  # todo translate earlier
                "sequence flanking",
                "guide+PAM sequence raw",
                "guide sequence raw",
                # 'guide strand',
                # 'PAM',"rPAM","is a reverse complement",'reverse complement'
            ],
            axis=1,
        )
    )


## base editing
def get_windows_seq(
    s: str,
    l: str,
    wl: str,
    verbose: bool = False,
) -> str:
    """Sequence by guide strand

    Args:
        s (str): sequence
        l (str): locus
        wl (str): window locus
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        str: window sequence
    """
    wl = parse_locus(wl)
    pos = get_pos(
        s=s,
        l=l,
    )
    if verbose:
        print(pos.to_dict())
    # return pos.loc[wl[1]+(0 if wl[3]=='+' else 0):wl[2]+(-1 if wl[3]=='+' else 0)].sum() # ordered by positions
    return pos.loc[wl[1] + 1 : wl[2]].sum()  # ordered by positions


def filter_guides(
    df1: pd.DataFrame,
    cfg_method: dict,
    verbose: bool = False,
) -> pd.DataFrame:
    """Filter sgRNAs

    Args:
        df1 (pd.DataFrame): input table
        cfg_method (dict): config of the method
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        pd.DataFrame: output table
    """
    df1 = df1.assign(
        **{
            "window locus": lambda df: df.apply(
                lambda x: to_locusby_pam(
                    chrom=x["chrom"],
                    pam_start=x["start PAM"],
                    pam_end=x["end PAM"],
                    pam_position=get_alt(["up", "down"], cfg_method["PAM position"])
                    if x["guide strand"] == "-"
                    else cfg_method["PAM position"],
                    strand=x["guide strand"],
                    length=cfg_method["distance of mutation from PAM: maximum"]-1,
                    start_off=cfg_method["distance of mutation from PAM: minimum"]-1,
                ),
                axis=1,
            ),
            "window sequence": lambda df: df.apply(
                lambda x: get_windows_seq(
                    s=x["guide sequence"],
                    l=x["guide locus"],
                    wl=x["window locus"],
                    verbose=verbose,
                ),
                axis=1,
            ),
            "window editable": lambda df: df["window sequence"].apply(
                lambda x: cfg_method["nucleotide"] in x
            ),
        },
    )
    return df1.log.query(expr="`window editable`==True"), df1.log.query(
        expr="`window editable`==False"
    )


def get_window_target_overlap(
    tstart: int,
    tend: int,
    wl: str,
    ws: str,
    nt: str,
    verbose: bool = False,
) -> tuple:
    """Get window target overlap

    Args:
        tstart (int): target start
        tend (int): target end
        wl (str): window locus
        ws (str): window sequence
        nt (str): nucleotide
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        tuple: window_overlaps_the_target,wts,nt_in_overlap,wtl
    """
    wpos = get_pos(
        ws,
        wl,
        reverse=False,
    )
    if verbose:
        print(wpos.to_dict())

    ## overlapping pos
    wtpos = wpos.loc[tstart:tend]  # because target is 1-based

    if len(wtpos) == 0:
        return False, None, None, None
    window_overlaps_the_target = True

    ## overlapping sequence
    wts = wtpos.sum()

    ## nt in overlapping seq
    if nt not in wts:
        return window_overlaps_the_target, wts, False, None
    nt_in_overlap = True

    from roux.lib.str import get_bracket

    ## overlapping locus
    wtl = to_locus(
        wl.split(":")[0],
        wtpos.index.min() - 1,  ## 0 based locus
        wtpos.index.max(),
        get_bracket(wl),
    )
    return window_overlaps_the_target, wts, nt_in_overlap, wtl


def get_mutated_codon(
    ts: str,
    tl: str,
    tes: str,
    tel: str,
    strand: str,
    verbose: bool = False,
) -> str:
    """Get mutated codon

    Args:
        ts (str): target sequence
        tl (str): target locus
        tes (str): target edited sequence
        tel (str): target edited locus
        strand (str): strand
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        str: mutated codon
    """

    ## get location of the edited sequence within the codon
    tp = get_pos(
        ts,
        tl,
        zero_based=False,
    )
    ##
    ep = get_pos(
        tes,
        tel,
        zero_based=True,
    )

    if verbose:
        print(tp.to_dict(), ep.to_dict())
    return (
        pd.concat([tp, ep])
        .reset_index()
        .drop_duplicates(subset=["index"], keep="last")
        .set_index("index")[0]
        .sort_values(ascending=strand != "+")
    ).sum()


def get_coedits_base(
    ws: str,
    wl: str,
    wts: str,
    wtl: str,
    nt: str,
    verbose: bool = False,
) -> str:
    """Get co-edited bases

    Args:
        ws (str): window sequence
        wl (str): window locus
        wts (str): window target overlap sequence
        wtl (str): window target overlap locus
        nt (str): nucleotide
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        str: coedits
    """
    chrom_, _, _, strand_ = parse_locus(wl)
    chrom, start, end, strand = parse_locus(wtl)

    assert chrom == chrom_, (chrom, chrom_)
    assert strand == strand_, (strand, strand_)

    wpos = get_pos(
        ws,
        wl,
        reverse=False,
    )
    ## overlap
    wtpos = get_pos(
        wts,
        wtl,
        reverse=False,
    )
    ## validate the editable nt position
    if verbose:
        print(wl, wpos.to_dict(), wtl, wtpos.to_dict())
    assert any([wtpos[i] == nt for i in range(start + 1, end + 1)]), (
        wl,
        wpos.to_dict(),
        wtl,
        wtpos.to_dict(),
    )
    wnott_pos = wpos.loc[[i for i in wpos.index if i not in wtpos.index]]

    coedits = wnott_pos.loc[lambda x: x == nt]
    return ";".join(
        [to_locus(chrom, k - 1, k, strand) for k, v in coedits.to_dict().items()]
    )
