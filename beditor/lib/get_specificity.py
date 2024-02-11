#!usr/bin/python
"""Specificities"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

import roux.lib.df as rd #noqa

from beditor.lib.utils import get_sequences, to_locus

def run_alignment(
    src_path: str,
    genomep: str,
    guidesfap: str,
    guidessamp: str,  # output
    guidel: int,  # length
    mismatches_max: int = 2,
    threads: int = 1,
    force: bool = False,
    verbose: bool = False,
) -> str:
    """Run alignment

    Args:
        src_path (str): source path
        genomep (str): genome path
        guidesfap (str): guide fasta path
        guidessamp (str): guide sam path
        threads (int, optional): threads. Defaults to 1.
        force (bool, optional): force. Defaults to False.
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        str: alignment file.
    """
    ## validate inputs
    assert isinstance(threads, int)

    from beditor.lib.utils import runbashcmd

    guidessap = str(Path(guidessamp).with_suffix("")) + ".sa"
    if not Path(genomep + ".bwt").exists():
        genome_is_large = Path(genomep).stat().st_size > 200000000
        if genome_is_large:
            logging.info("large genome detected using '-a bwtsw' for indexing")
        cmd = f"{src_path} index {'-a bwtsw' if genome_is_large else ''} {genomep} 2> {genomep}.log"
        logging.debug(cmd)
        if verbose:
            print("Running `bwa index` to index database sequences in the FASTA format")
            print(
                "Note: Depending on the size of the genome, this step can take a while."
            )
            if genome_is_large:
                print("Note: For human genome, this time could be in hours.")
                print(
                    "Note: This step is necessary only for the first-time use of a genome."
                )
                print(
                    "Note: At later use of the genome, beditor will use the already indexed file."
                )
                print(
                    f"Note: This step can be skipped if pre-indexed BWA indices are available, placed next to the genome file ({genomep})."
                )
        runbashcmd(cmd)

    ## BWA alignment command is adapted from cripror
    ## crispor.py
    # BWA: allow up to X mismatches
    # maximum number of occurences in the genome to get flagged as repeats.
    # This is used in bwa samse, when converting the sam file
    # and for warnings in the table output.
    # upto MAXOCC matching alignments (e.g 20M) are stored in the XA tag
    MAXOCC = 60000
    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000 / MAXOCC
    # increase MAXOCC if there is only a single query, but only in CGI mode
    bwaM = MFAC * MAXOCC  # -m is queue size in bwa
    # Find the SA coordinates of the input reads. Maximum maxSeedDiff differences are allowed in the first seedLen subsequence and maximum maxDiff differences are allowed in the whole sequence.
    cmd = f"{src_path} aln -t 1 -o 0 -m {bwaM} -n {mismatches_max} -k {mismatches_max} -N -l {guidel} {genomep} {guidesfap} -t {threads} > {guidessap} 2> {guidessap}.log"
    logging.debug(cmd)
    if verbose:
        print("Running alignment of the sgRNAs")
    runbashcmd(cmd)
    # Generate alignments in the SAM format given single-end reads.
    # if verbose:
    #     print("Running alignment of the sgRNAs")
    cmd = f"{src_path} samse -n {MAXOCC} {genomep} {guidessap} {guidesfap} > {guidessamp} 2> {guidessamp}.log"
    logging.debug(cmd)
    runbashcmd(cmd)
    return guidessamp


def read_sam(
    align_path: str,
) -> pd.DataFrame:
    """read alignment file

    Args:
        align_path (str): path to the alignment file

    Returns:
        pd.DataFrame: output table

    Notes:
        Tag	Meaning
        NM	Edit distance
        MD	Mismatching positions/bases
        AS	Alignment score
        BC	Barcode sequence
        X0	Number of best hits
        X1	Number of suboptimal hits found by BWA
        XN	Number of ambiguous bases in the referenece
        XM	Number of mismatches in the alignment
        XO	Number of gap opens
        XG	Number of gap extentions
        XT	Type: Unique/Repeat/N/Mate-sw
        XA	Alternative hits; format: (chr,pos,CIGAR,NM;)*
        XS	Suboptimal alignment score
        XF	Support from forward/reverse alignment
        XE	Number of supporting seeds

    Reference:
        https://bio-bwa.sourceforge.net/bwa.shtml
    """
    ## alignment values
    from Bio.Align.sam import AlignmentIterator

    aligns = {}
    for a in AlignmentIterator(align_path):
        if a.target is None:
            logging.warning(f"a.target is None, a.query={a.query}")
            continue
        aligns[a.query.id] = pd.Series(
            {
                **{
                    "mapq": a.mapq,
                    "chrom": a.target.id,
                    "start": a.target.seq.defined_ranges[0][0],
                    "end": a.target.seq.defined_ranges[0][1],
                    "sequence": str(
                        a.target.seq[
                            a.target.seq.defined_ranges[0][
                                0
                            ] : a.target.seq.defined_ranges[0][1]
                        ]
                    ),
                },
                **a.annotations,
            }
        )

    return (
        pd.concat(aligns, names=["guide sequence"])
        .unstack()
        .add_prefix("aligned ")
        .reset_index()
        .drop_duplicates()
    )


## alignements in sam format
def parse_XA(
    XA: str,
) -> pd.DataFrame:
    """Parse XA tags

    Args:
        XA (str): XA tag

    Notes:
        format: (chr,pos,CIGAR,NM;)

    Example:
        XA='4,+908051,23M,0;4,+302823,23M,0;4,-183556,23M,0;4,+1274932,23M,0;4,+207765,23M,0;4,+456906,23M,0;4,-1260135,23M,0;4,+454215,23M,0;4,-1177442,23M,0;4,+955254,23M,1;4,+1167921,23M,1;4,-613257,23M,1;4,+857893,23M,1;4,-932678,23M,2;4,-53825,23M,2;4,+306783,23M,2;'
    """
    import re

    # Split the string into individual elements based on the semicolon delimiter
    elements = XA.split(";")

    parsed_data = []
    for element in elements:
        # Split each element based on the comma delimiter
        parts = element.split(",")
        # Extract the values for chr, pos, CIGAR, and NM
        if len(parts) == 4:
            chr_val = parts[0]
            pos_val = parts[1]
            cigar_val = parts[2]
            nm_val = parts[3]
            # check match
            assert re.match(r"^\d+M$", cigar_val), cigar_val
            strand = pos_val[0]
            start = int(int(pos_val[1:]) + (-1 if strand == "-" else 0))
            length = int(cigar_val[:-1])
            end = int(start + length)
            assert (end - start) == length, (start, end, strand, length)
            parsed_data.append(
                dict(
                    chrom=chr_val,
                    start=start,
                    end=end,
                    strand=strand,
                    cigar=cigar_val,
                    NM=int(nm_val),
                )
            )
    return pd.DataFrame(parsed_data)


def get_extra_alignments(
    df1: pd.DataFrame,  # from read_bam
    genome: str,  ## to get sequences
    bed_path: str,
    alignments_max: int = 10,  ## how many imperfect alignments are considered, ootherwise the guide is dropped
    threads: int = 1,
) -> pd.DataFrame:
    """Get extra alignments

    Args:
        df1 (pd.DataFrame): input table
        alignments_max (int, optional): alignments max. Defaults to 10.
        threads (int, optional): threads. Defaults to 1.

    Returns:
        pd.DataFrame: output table

    TODOs:
        1. apply parallel processing to get_seq
    """
    ## get the alignment counts and filter
    df2 = (
        df1.assign(
            **{
                "alignments": lambda df: df["aligned XA"].apply(
                    lambda x: x.count(";") if isinstance(x, str) else 0
                ),
            }
        )
        .log("guide sequence")
        .log.query(expr=f"`alignments` <= {alignments_max}")
        .log("guide sequence")
    )
    # df2
    logging.info("separating extra alignments")
    cols_groupby = ["guide sequence", "aligned XA", "alignments"]
    df2_ = (
        df2.groupby(cols_groupby)
        .apply(lambda df: parse_XA(df.name[1]))
        .assign(
            **{
                "locus": lambda df: df.apply(
                    lambda x: to_locus("chrom", "start", "end", "strand", x=x), axis=1
                ),
            }
        )
        .reset_index()
        .rd.clean()
    )
    logging.info(f"{len(df2_)} extra alignments found")
    logging.info("could be a slow step depending on the genome size")

    df21_ = df2_.assign(
        # contig=lambda df: df.apply(lambda x: x['chrom'],axis=1),
        start=lambda df: df.apply(
            lambda x: x["start"] + (0 if x["strand"] == "+" else 1), axis=1
        ),
        end=lambda df: df.apply(
            lambda x: x["end"] - 1 + (0 if x["strand"] == "+" else 1), axis=1
        ),
        # strand=lambda df: df.apply(lambda x: x['strand'],axis=1),
    ).log.drop_duplicates()
    # logging.info(f"saving bed file with {len(df21_)} loci")
    df2_ = (
        get_sequences(
            df21_,
            # renames={
            # },
            src_path=None,
            genome_path=genome,
            p=bed_path,
            revcom=True,
        )
        .merge(
            right=df2_.drop(["chrom", "strand", "start", "end"], axis=1).reset_index(),
            on="locus",
            how="left",
            validate="m:m",
        )
        .set_index(cols_groupby)
        .add_prefix("aligned ")
        .reset_index()
    )

    assert not df2_["aligned sequence"].isnull().any()
    return pd.concat([df2, df2_], axis=0, ignore_index=True).sort_values(
        "guide sequence"
    )


def to_pam_coord(
    pam_pos: str,
    pam_len: int,
    align_start: int,
    align_end: int,
    strand: str,
) -> tuple:
    """Get PAM coords

    Args:
        pam_pos (str): PAM position
        pam_len (int): PAM length
        align_start (int): alignment start
        align_end (int): alignment end
        strand (str): strand

    Returns:
        tuple: start,end
    """
    if pam_pos == "up":
        # PAM---
        if strand == "+":
            # PAM---
            start, end = (align_start - pam_len) + 1, align_start
        elif strand == "-":
            # ---PAM
            start, end = align_end + 1, align_end + pam_len
        else:
            raise ValueError(strand)
    elif pam_pos == "down":
        # ---PAM
        if strand == "+":
            # ---PAM
            start, end = align_end + 1, align_end + pam_len
        elif strand == "-":
            # PAM---
            start, end = (align_start - pam_len) + 1, align_start
        else:
            raise ValueError(strand)
    else:
        raise ValueError(pam_pos)
    return start, end  # f"{chrom}:{int(start)}-{int(end)}"


def get_alignments(
    align_path: str,
    genome: str,  ## to get the sequence of the multi-alignments
    alignments_max: int,
    ## for getting the sequence of the aligned PAM
    pam_pos: str,
    pam_len: int,
    guide_len: int,
    pam_pattern: str,
    pam_bed_path: str,  # path with
    extra_bed_path: str,  # path with
    # fast=False,
    **kws_xa,
) -> pd.DataFrame:
    """Get alignments

    Args:
        align_path (str): alignement path
        genome (str): genome path
        pam_pos (str): PAM position
        pam_len (int): PAM length
        guide_len (int): sgRNA length
        pam_pattern (str): PAM pattern
        pam_bed_path (str): PAM bed path

    Returns:
        pd.DataFrame: output path
    """
    ## read the bam
    df1 = read_sam(align_path)
    ## get multi-alignments
    if "aligned strand" not in df1:
        logging.info("getting the `aligned strand`")
        # sequence is always fetched from + strand for the non XA alignments
        import pysam

        # pysam.sort("-o", cfg['sgRNAs']['url'], f"{output_dir_path}/04_offtargets/alignment.sam")
        # pysam.index(cfg['sgRNAs']['url'])
        # samfile = pysam.AlignmentFile(Path(align_path).with_suffix('.bam'), "rb")
        samfile = pysam.AlignmentFile(align_path, "rb")
        to_strand = {}
        for read in samfile.fetch():
            to_strand[read.qname] = "+" if read.is_forward else "-"
        from roux.lib.set import assert_overlaps_with

        assert_overlaps_with(df1["guide sequence"].tolist(), list(to_strand.keys()))
        df1 = df1.assign(
            **{
                "aligned strand": lambda df: df["guide sequence"].map(to_strand),
            }
        )
    if "aligned XA" not in df1:
        logging.info("no multiple alignments found.")
        df1 = df1.assign(alignments=0)
    else:
        logging.info(f"getting the `multi_alignments`, alignments_max={alignments_max}")
        # tmp
        from roux.lib.io import to_table

        to_table(df1, "test/extra_aligns.tsv")
        df1 = get_extra_alignments(
            df1,  # from read_bam
            genome=genome,  ## to get sequences
            alignments_max=alignments_max,
            bed_path=extra_bed_path,
            # fast=fast,
            **kws_xa,
        )
    ## get the aligned pam sequence
    import re

    df1 = df1.assign(
        **{
            "aligned PAM start,end": lambda df: df.apply(
                lambda x: to_pam_coord(
                    pam_pos=pam_pos,
                    pam_len=pam_len,
                    align_start=x["aligned start"],
                    align_end=x["aligned end"],
                    strand=x["aligned strand"],
                ),
                axis=1,
            ),
            "aligned PAM start": lambda df: df["aligned PAM start,end"].apply(
                pd.Series
            )[0],
            "aligned PAM end": lambda df: df["aligned PAM start,end"].apply(pd.Series)[
                1
            ],
            "aligned PAM location": lambda df: df.apply(
                lambda x: to_locus(
                    "aligned chrom",
                    "aligned PAM start",
                    "aligned PAM end",
                    "aligned strand",
                    x=x,
                ),
                axis=1,
            ),
        }
    )
    ## get pam sequence
    df2 = (
        get_sequences(
            df1,
            renames={
                "aligned chrom": "chrom",
                "aligned PAM start": "start",
                "aligned PAM end": "end",
                "aligned PAM location": "locus",
                "aligned strand": "strand",
            },
            src_path=None,
            genome_path=genome,
            p=pam_bed_path,
            revcom=True,
        )
        .rename(
            columns={"sequence": "aligned PAM"},
            errors="raise",
        )
        .loc[:, ["aligned PAM location", "aligned PAM"]]
        .merge(
            right=df1,
            on="aligned PAM location",
            how="left",
            validate="m:m",
        )
    )
    assert all(df2["aligned sequence"].apply(len) == guide_len), sum(
        df2["aligned sequence"] == guide_len
    )

    df2 = df2.assign(
        **{
            "aligned PAM valid": lambda df: df["aligned PAM"].apply(
                lambda x: bool(re.match(r"^{}$".format(pam_pattern), x))
            ),
        }
    ).log.drop(labels=["aligned PAM start,end"], axis=1)
    return df2.log.query(expr="`aligned PAM valid`==True"), df2.log.query(
        expr="`aligned PAM valid`==False"
    )


def get_penalties(
    aligns: pd.DataFrame,
    guides: pd.DataFrame,
    annots: pd.DataFrame,
) -> pd.DataFrame:
    """Get penalties

    Args:
        aligns (pd.DataFrame): alignements
        guides (pd.DataFrame): guides
        annots (pd.DataFrame): annotations

    Returns:
        pd.DataFrame: output table
    """
    ### Target is genic or intergenic region
    from beditor.lib.utils import align

    df3 = (
        aligns.merge(
            right=guides,
            # .loc[:,['guide sequence',
            # 'PAM position','guide strand' ##todo: remove unused
            # ]].drop_duplicates(),
            how="inner",
            on="guide sequence",
        )
        .astype(
            {
                "aligned start": int,
                "aligned end": int,
            }
        )
        .assign(
            **{
                "aligned genes": lambda df: df.apply(
                    lambda x: list(
                        annots.gene_ids_at_locus(
                            contig=x["aligned chrom"],
                            position=x["aligned start"],
                            end=x["aligned end"],
                        )
                    ),
                    axis=1,
                ),
                "aligned genic": lambda df: df["aligned genes"]
                .fillna("")
                .apply(
                    lambda x: True if isinstance(x, list) and len(x) != 0 else False
                ),
                "alignment": lambda df: df.apply(
                    lambda x: align(
                        x["guide sequence"],
                        x["aligned sequence"],
                    ),
                    axis=1,
                ),
                "alignment is perfect": lambda df: df["alignment"].apply(
                    lambda x: x == (len(x) * "|")
                ),
            }
        )
    )
    # print(df3.columns.tolist())
    df4 = (
        df3
        ## mark offtargets by comparing the location of guide and the aligned
        .assign(
            **{
                "aligned locus": lambda df: df.apply(
                    lambda x: to_locus(
                        "aligned chrom",
                        "aligned start",
                        "aligned end",
                        "aligned strand",
                        x=x,
                    ),
                    axis=1,
                ),
                "offtarget": lambda df: ~(
                    (df["guide locus"] == df["aligned locus"])
                    & df["alignment is perfect"]
                ),
            }
        )
    )
    # assert sum(~df4['offtarget'])!=0, "on target alignments are expected."

    ## check on targets
    df_ = df4.query("`offtarget`==False")
    assert df_["aligned PAM valid"].all(), df_["aligned PAM valid"].sum()
    assert (df_["aligned PAM"] == df_["PAM sequence"]).all(), (
        df_["aligned PAM"] == df_["PAM sequence"]
    ).sum()

    logging.info(df4["offtarget"].value_counts().to_string())
    return df4


def score_alignments(
    df4: pd.DataFrame,  # penalties
    # aligns,
    # annots,
    pam_len: int,
    pam_pos: str,
    # guide_len,
    pentalty_genic: float = 0.5,
    pentalty_intergenic: float = 0.9,
    pentalty_dist_from_pam: float = 0.1,
    verbose: bool = False,
) -> tuple:
    """score_alignments _summary_

    Args:
        df4 (pd.DataFrame): input table
        pam_pos (str): PAM position
        pentalty_genic (float, optional): penalty for offtarget in genic locus. Defaults to 0.5.
        pentalty_intergenic (float, optional): penalty for offtarget in intergenic locus. Defaults to 0.9.
        pentalty_dist_from_pam (float, optional): penalty for offtarget wrt distance from PAM. Defaults to 0.1.
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        tuple: tables

    Note:
        1. Low value corresponds to high penalty and vice versa, because values are multiplied.
        2. High penalty means consequential offtarget alignment and vice versa.
    """
    df5 = (
        df4.log("guide sequence")
        .groupby("guide sequence")
        .filter(lambda df: not df["offtarget"].all())
        .log("guide sequence")
    )
    assert len(df5) != 0, "on target alignments are expected."
    ## aligner did not report an on target alignment
    df5_ = (
        df4.log("guide sequence")
        .groupby("guide sequence")
        .filter(lambda df: df["offtarget"].all())
        .log("guide sequence")
    )

    ## Test if 1 alignment is an exact match to the designed guide
    assert (
        df5.groupby("guide sequence")
        .apply(lambda df: df["alignment is perfect"].any())
        .all()
    ), (
        df5.groupby("guide sequence")
        .apply(lambda df: df["alignment is perfect"].any())
        .sum()
    )
    # logging.debug(df5['aligned XT'].value_counts().to_string())
    # logging.warning("not all guides have at least one perfect alignment. potentially repeatitive sequence if aligned XT=R and mapq=0")
    # if not all(df4['aligned NM']==0):
    #     logging.debug(df5['aligned NM'].value_counts().to_string())

    ## get score
    from beditor.lib.get_scores import get_beditorscore_per_alignment

    df6 = df5.assign(
        **{
            "penalty mismatches": lambda df: df.apply(
                lambda x: get_beditorscore_per_alignment(
                    NM=x["aligned NM"],
                    alignment=x["alignment"],
                    pam_len=pam_len,
                    pam_pos=pam_pos,
                    pentalty_dist_from_pam=pentalty_dist_from_pam,
                    verbose=verbose,
                ),
                axis=1,
            ),
            "score per alignment": lambda df: df["penalty mismatches"]
            * df.apply(
                lambda x: (
                    1
                    if not x["offtarget"]
                    else pentalty_genic
                    if x["aligned genic"]
                    else pentalty_intergenic
                ),
                axis=1,
            ),
        },
    )

    assert all(
        df6.query("`offtarget`==False")["score per alignment"].unique() == np.array([1])
    ), df6.query("`offtarget`==False")["score per alignment"].unique()

    return df6, df5_


def score_guides(
    guides: pd.DataFrame,
    scores: pd.DataFrame,
    # penalty_activity_window=0.5,
    not_be: bool = False,
) -> pd.DataFrame:
    """Score guides

    Args:
        guides (pd.DataFrame): guides
        scores (pd.DataFrame): scores
        not_be (bool, optional): not a base editor. Defaults to False.

    Returns:
        pd.DataFrame: output table

    Changes:
        penalty_activity_window disabled as only the sgRNAs with target in the window are reported.
    """

    def agg_lists(ds):
        # try:
        return ";".join(
            list(
                set(
                    sum(
                        [
                            eval(s) if not isinstance(s, list) else s
                            for s in ds.tolist()
                        ],
                        [],
                    )
                )
            )
        )
        # except:
        #     print(ds.tolist())

    ## by number of alignmnets
    df1 = (
        scores.groupby("guide sequence")
        .agg(
            **{
                "score": ("score per alignment", "prod"),
                "offtargets": ("offtarget", "sum"),
                "aligned genes": ("aligned genes", lambda x: agg_lists(x)),
                "aligned locuss": (
                    "aligned locus",
                    lambda x: ";".join(list(set(x.tolist()))),
                ),
            },
        )
        .reset_index()
        .merge(
            right=guides.loc[:, ["guide locus", "guide sequence"]].drop_duplicates(),
            on="guide sequence",
            how="left",
        )
        .merge(
            right=scores.loc[:, ["guide sequence", "alignments"]]
            .drop_duplicates()
            .dropna()
            .rd.assert_no_dups(subset=["guide sequence"]),
            on="guide sequence",
            how="left",
        )
    )
    # if not_be:
    ## no penalty
    return df1
    # else:
    #     return (df1
    #         .assign(
    #             **{
    #                 ## whther the target is within window
    #                 'score':lambda df: df.apply(lambda x: x['score']*(penalty_activity_window if not x['target is in window'] else 1),axis=1),
    #             },
    #             )
    #            )
