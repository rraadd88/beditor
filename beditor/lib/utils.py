#!usr/bin/python
"""Utilities"""

import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

from roux.lib.sys import makedirs
import roux.lib.df as rd #noqa

## variables
cols_muts = ["chrom", "start", "end", "strand"]


## sys
def get_src_path() -> str:
    """Get the beditor source directory path.

    Returns:
        str: path
    """
    import site

    src_path = f"{site.getsitepackages()[0]}/beditor/"
    if not Path(src_path).exists():
        logging.warning("package is installed in the development mode")
        # import beditor
        # src_path=str(Path(beditor.__file__).parent)
        src_path = str(Path(os.path.realpath(__file__)).parent.parent)
    return src_path

# # EXT
def runbashcmd(cmd: str, test: bool = False, logf=None):
    """Run a bash command

    Args:
        cmd (str): command
        test (bool, optional): test-mode. Defaults to False.
        logf (optional): log file instance. Defaults to None.
    """
    import subprocess
    if test:
        print(cmd)
    err = subprocess.call(cmd, shell=True, stdout=logf, stderr=subprocess.STDOUT)
    if err != 0:
        import sys
        print("bash command error: {}\n{}\n".format(err, cmd))
        sys.exit(1)

def log_time_elapsed(start):
    """Log time elapsed.

    Args:
        start (datetime): start tile

    Returns:
        datetime: difference in time.
    """
    from datetime import datetime
    diff = datetime.now() - start
    return diff


## arrays
def rescale(a: np.array, mn: float = None) -> np.array:
    """Rescale a vector.

    Args:
        a (np.array): vector.
        mn (float, optional): minimum value. Defaults to None.

    Returns:
        np.array: output vector
    """
    a = (a - a.min()) / (a.max() - a.min())
    if mn is not None:
        a = 1 - a
        a = a * (1 - mn)
        a = 1 - a
    return a


## strings
## sequences
# # # regex
multint2reg = {
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ATGC]",
}
multint2regcomplement = {
    "R": "[TC]",
    "Y": "[GA]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[CA]",
    "M": "[TG]",
    "B": "[^T]",
    "D": "[^G]",
    "H": "[^C]",
    "V": "[^A]",
    "N": "[ATGC]",
}


def get_nt2complement():
    nt2complement = {
        "A": "T",
        "G": "C",
        "N": "N",
        "R": "Y",
        "S": "W",
        "K": "M",
        "B": "b",
        "D": "d",
        "H": "h",
        "V": "v",
        "N": "N",
    }
    nt2complement.update(dict(zip(nt2complement.values(), nt2complement.keys())))
    return nt2complement


def s2re(
    s: str,
    ss2re: dict,
) -> str:
    """String to regex patterns

    Args:
        s (str): string
        ss2re (dict): substrings to regex patterns.

    Returns:
        str: string with regex patterns.
    """
    for ss in ss2re:
        s = s.replace(ss, ss2re[ss])
    return s


def parse_locus(
    s: str,
    zero_based: bool = True,
) -> tuple:
    """parse_locus _summary_

    Args:
        s (str): location string.
        zero_based (bool, optional): zero-based coordinates. Defaults to True.

    Returns:
        tuple: chrom, start, end, strand

    Notes:
        beditor outputs (including bed files) use 0-based loci
        pyensembl and IGV use 1-based locations
    """
    loc = (
        s.split(":")[0],
        int(s.split(":")[1].split("-")[0]),
        int(s.split("-")[1].split("(")[0]),
        s.split("(")[1].split(")")[0],
    )
    if zero_based:
        return loc
    else:
        if loc[1] != loc[2]:
            return loc[0], loc[1] + 1, loc[2], loc[3]
        elif loc[1] == loc[2]:
            return loc[0], loc[1] + 1, loc[2] + 1, loc[3]


def get_pos(
    s: str,
    l: str,
    reverse: bool = True,
    zero_based: bool = True,
) -> pd.Series:
    """Expand locus to positions mapped to nucleotides.

    Args:
        s (str): sequence
        l (str): locus
        reverse (bool, optional): reverse the - strand. Defaults to True.
        zero_based (bool, optional): zero based coordinates. Defaults to True.

    Returns:
        pd.Series: output.
    """
    chrom, start, end, strand = parse_locus(l)
    if not zero_based:
        start -= 1
    if zero_based and (end - start) == len(s) - 1:
        logging.error(
            "maybe 1-based location provided, zero_based=False to resolve this error."
        )
    assert (end - start) == len(
        s
    ), f"0-based location expected ({l},length of sequence={len(s)})"
    if strand == "+":
        pos = dict(zip(range(start + 1, end + 1), list(s)))
    elif strand == "-":
        if reverse:
            # reverse complemented sequence
            pos = dict(zip(range(start + 1, end + 1), list(s)[::-1]))
        else:
            # because the sequence already from start to end
            pos = dict(zip(range(start + 1, end + 1), list(s)))
    return pd.Series(pos)

import gc

import numpy as np
import pandas as pd
from Bio.Seq import Seq as str2seq


def get_seq(
    genome: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    out_type: str = "str",
    verbose: bool = False,
) -> str:
    """Extract a sequence from a genome file based on start and end positions using streaming.

    Args:
        genome (str): The path to the genome file in FASTA format.
        contig (str): chrom
        start (int): start
        end (int): end
        strand (str): strand
        out_type (str, optional): type of the output. Defaults to 'str'.
        verbose (bool, optional): verbose. Defaults to False.

    Raises:
        ValueError: invalid strand.

    Returns:
        str: The extracted sequence.
    """
    if verbose:
        logging.debug(f"getting {contig}:{start}-{end}{strand}")
    # for r in genome:
    #     print(r.name)
    #     if contig==r.name:
    #         s=r.seq[start - 1:end]
    # logging.error(f"contig ({contig}) not found")
    contig = str(contig)
    if strand == "+":
        s = genome[contig][start:end].seq
    elif strand == "-":
        s = genome[contig][start:end].reverse_complement().seq
    else:
        raise ValueError(strand)
    # to avoid memory build up
    del genome
    gc.collect()
    if out_type == "str":
        return str(s)
    else:
        return s


### multiple seq fasta
def read_fasta(
    fap: str,
    key_type: str = "id",
    duplicates: bool = False,
    out_type="dict",
) -> dict:
    """Read fasta

    Args:
        fap (str): path
        key_type (str, optional): key type. Defaults to 'id'.
        duplicates (bool, optional): duplicates present. Defaults to False.

    Returns:
        dict: data.

    Notes:
        1. If `duplicates` key_type is set to `description` instead of `id`.
    """

    def post(id2seq, out_type):
        if out_type == "dict":
            return id2seq
        else:
            from roux.lib.str import get_bracket

            return (
                pd.DataFrame([{"locus": k, "sequence": v} for k, v in id2seq.items()])
                .assign(
                    **{
                        "strand": lambda df: df["locus"].apply(get_bracket),
                        # "locus":lambda df: df['locus'].apply(lambda x: x.split('(')[0]),
                        "chrom": lambda df: df["locus"].str.split(":", expand=True)[0],
                        "start": lambda df: df["locus"]
                        .apply(lambda x: x.split(":")[1].split("-")[0])
                        .astype(int),
                        "end": lambda df: df["locus"]
                        .apply(lambda x: x.split(":")[1].split("-")[1].split("(")[0])
                        .astype(int),
                    }
                )
                .loc[:, ["chrom", "start", "end", "locus", "strand", "sequence"]]
            )

    from Bio import SeqIO

    if (not duplicates) or key_type == "id":
        try:
            id2seq = SeqIO.to_dict(SeqIO.parse(fap, format="fasta"))
            id2seq = {k: str(id2seq[k].seq) for k in id2seq}
            return post(id2seq, out_type)
        except:
            duplicates = True
    if duplicates or key_type == "description":
        id2seq = {}
        for seq_record in SeqIO.parse(fap, "fasta"):
            id2seq[getattr(seq_record, key_type)] = str(seq_record.seq)
        return post(id2seq, out_type)


def format_coords(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Format coordinates

    Args:
        df (pd.DataFrame): table

    Returns:
        pd.DataFrame: formated table
    """
    return df.astype(
        dict(
            chrom=str,
            start=int,
            end=int,
            strand=str,
        )
    )


def fetch_sequences_bp(
    p: str,
    genome: str,
) -> pd.DataFrame:
    """Fetch sequences using biopython.

    Args:
        p (str): path to the bed file.
        genome (str): genome path.

    Returns:
        pd.DataFrame: sequences.
    """
    if isinstance(genome, str):
        from beditor.lib.io import read_genome

        genome = read_genome(genome, fast=True)
    from beditor.lib.io import read_bed

    return (
        read_bed(p)
        .pipe(format_coords)
        .assign(
            **{
                "sequence": lambda df: df.apply(
                    lambda x: get_seq(
                        # genome=genome,
                        genome=genome,
                        contig=x["chrom"],
                        start=x["start"],
                        end=x["end"],
                        strand=x["strand"],
                    ),
                    axis=1,
                ).apply(str)
            }
        )
    )


def fetch_sequences(
    p: str,
    genome_path: str,
    outp: str = None,
    src_path: str = None,
    revcom: bool = True,
    method="2bit",
    out_type="df",
) -> pd.DataFrame:
    """Fetch sequences

    Args:
        p (str): path to the bed file
        genome_path (str): genome path
        outp (str, optional): output path for fasta file. Defaults to None.
        src_path (str, optional): source path. Defaults to None.
        revcom (bool, optional): reverse-complement. Defaults to True.
        method (str, optional): method name. Defaults to '2bit'.
        out_type (str, optional): type of the output. Defaults to 'df'.

    Returns:
        pd.DataFrame: sequences.
    """
    # todo revcom=True for fetching strand-aware sequences i.e. revcom for - strand
    if outp is None:
        outp = Path(p).with_suffix(".fa").as_posix()
    if not isinstance(genome_path, str):
        df1 = fetch_sequences_bp(
            p,
            genome_path,
        )
    else:
        from beditor.lib.utils import runbashcmd
        
        if src_path is None:
            import site
            src_path=Path(site.getsitepackages()[0]).parent.parent.parent / 'bin/'
            #todo change 
            # src_path = "micromamba run -n twobit twoBitToFa"

        genome_2bit_path = f"{genome_path.split('.fa')[0]}.2bit"
        if method == "2bit" and Path(genome_2bit_path).exists() and Path(f"{str(src_path)}/twoBitToFa").exists():
            logging.info("using 2bit reference")
            cmd = f"{src_path} -bed={p} {genome_2bit_path} {outp}"
            runbashcmd(cmd)
            print(outp)
            df1 = read_fasta(outp, out_type=out_type)
        elif (
            method == "bedtools" or not Path(genome_2bit_path).exists()
        ) and method != "biopython" and Path(f"{str(src_path)}/bedtools").exists():
            if method != "bedtools":
                logging.warning(f"falling back to bedtools (set method={method})")
            else:
                logging.info("using bedtools")
            cmd = f"{src_path}/bedtools getfasta {'-s' if revcom else ''} -name -fi {genome_path} -bed {p} -fo {outp}"
            runbashcmd(cmd)
            df1 = read_fasta(outp, out_type=out_type)
            df1 = df1.assign(
                locus=lambda df: df["locus"].str.rsplit("(", n=1, expand=True)[0]
            )
        else:
            if method != "biopython":
                logging.warning(f"falling back to biopython (set method={method})")
            else:
                logging.warning("using a slow method to fetch sequences")
            df1 = fetch_sequences_bp(
                p,
                genome_path,
            )
            if (df1["score"] == 0).all():
                df1 = df1.drop(["score"], axis=1)
    return df1


def get_sequences(
    df1: pd.DataFrame,
    p: str,  # bed path
    genome_path: str,
    outp: str = None,
    src_path: str = None,
    revcom: bool = True,
    out_type: str = "df",
    renames: dict = {},  # original to bed
    **kws_fetch_sequences,
) -> pd.DataFrame:
    """Get sequences for the loci in a table

    Args:
        df1 (pd.DataFrame): input table
        p (str): path to the beb file
        outp (str, optional): output path. Defaults to None.
        src_path (str, optional): source path. Defaults to None.
        revcom (bool, optional): reverse complement. Defaults to True.
        out_type (str, optional): output type. Defaults to 'df'.
        renames (dict, optional): renames. Defaults to {}.

    Returns:
        pd.DataFrame: output sequences

    Notes:
        Input is 1-based
        Output is 0-based
        Saves bed file and gets the sequences
    """
    from beditor.lib.io import to_bed

    to_bed(
        (
            df1.drop(labels=list(renames.values()), axis=1, errors="ignore")
            .rename(columns=renames, errors="raise")
            .assign(start=lambda df: df["start"] - 1)
        ),  ## 1-based to 0-based
        makedirs(p),
    )
    df2 = fetch_sequences(
        p=p,
        genome_path=genome_path,
        outp=outp,
        src_path=src_path,
        revcom=revcom,
        out_type=out_type,
        **kws_fetch_sequences,
    )

    if len(renames) != 0:
        from roux.lib.dict import flip_dict

        df2 = df2.rename(columns=flip_dict(renames), errors="raise")
    return df2


def to_locus(
    chrom: str = "chrom",
    start: str = "start",
    end: str = "end",
    strand: str = "strand",
    x: pd.Series = None,
) -> str:
    """To locus

    Args:
        chrom (str, optional): chrom. Defaults to 'chrom'.
        start (str, optional): strart. Defaults to 'start'.
        end (str, optional): end. Defaults to 'end'.
        strand (str, optional): strand. Defaults to 'strand'.
        x (pd.Series, optional): row of the dataframe. Defaults to None.

    Returns:
        str: locus
    """
    if x is None:
        return f"{chrom}:{start}-{end}({strand})"
    else:
        return f"{x[chrom]}:{x[start]}-{x[end]}({x[strand]})"


def get_flanking_seqs(
    df1: pd.DataFrame,
    targets_path: str,
    flanks_path: str,
    genome: str = None,
    search_window: list = None,
) -> pd.DataFrame:
    """Get flanking sequences

    Args:
        df1 (pd.DataFrame): input table
        targets_path (str): target sequences path
        flanks_path (str): flank sequences path
        genome (str, optional): genome path. Defaults to None.
        search_window (list, optional): search window around the target. Defaults to None.

    Returns:
        pd.DataFrame: output table with sequences
    """
    df1 = df1.log.drop_duplicates(subset=["target location"]).assign(
        **{
            "start flanking": lambda df: df["start"] - (search_window // 2),
            "end flanking": lambda df: df["end"] + (search_window // 2),
        }
    )
    ## get sequences
    df2 = get_sequences(
        df1,
        renames={
            "target location": "locus",  ## for merging
        },
        src_path=None,
        genome_path=genome,
        p=targets_path,
        revcom=True,
    )
    df3 = get_sequences(
        df1,
        renames={
            "target location": "locus",
            "start flanking": "start",
            "end flanking": "end",
        },
        src_path=None,
        genome_path=genome,
        p=flanks_path,
        revcom=True,
    ).drop(["start flanking", "end flanking"], axis=1)
    df4 = (
        df2.set_index(["target location", "chrom", "strand", "start", "end"])
        .add_suffix(" target")
        .reset_index()
        .log.merge(
            right=df3.set_index(["target location", "chrom", "strand"])
            .add_suffix(" flanking")
            .reset_index(),
            how="inner",
            on=["target location", "chrom", "strand"],
            validate="1:1",
            validate_equal_length=True,
            # suffixes=[' target',' flanking'],
        )
        .log.merge(
            right=df1.loc[
                :, ["target location", "strand", "start flanking", "end flanking"]
            ],
            how="inner",
            on="target location",
            validate="1:1",
            validate_equal_length=True,
            suffixes=["", " original"],
        )
    )
    assert all(df4["strand"] == df4["strand original"]), sum(
        ~(df4["strand"] == df4["strand original"])
    )
    return df4.drop(["strand original"], axis=1)


def get_strand(
    genome,
    df1: pd.DataFrame,
    col_start: str,
    col_end: str,
    col_chrom: str,
    col_strand: str,
    col_seq: str,
) -> pd.DataFrame:
    """
    Get strand by comparing the aligned and fetched sequence

    Args:
        genome: genome instance
        df1 (pd.DataFrame): input table.
        col_start (str): start
        col_end (str): end
        col_chrom (str): chrom
        col_strand (str): strand
        col_seq (str): sequences

    Returns:
        pd.DataFrame: output table

    Notes:
        used for tests.
    """
    from beditor.lib.utils import get_seq

    df1_ = df1.assign(
        **{
            "_fetched": lambda df: df.apply(
                lambda x: get_seq(
                    genome=genome,
                    contig=x[col_chrom],
                    start=x[col_start] + 1,
                    end=x[col_end],
                    strand="+",
                ),
                axis=1,
            ),
            col_strand: lambda df: df.apply(
                lambda x: "+" if x[col_seq] == x["_fetched"] else "-", axis=1
            ),
        }
    )
    assert (df1_[col_seq].apply(len) == df1_["_fetched"].apply(len)).all(), (
        df1_[col_seq].apply(len) == df1_["_fetched"].apply(len)
    ).sum()
    return df1


def reverse_complement_multintseq(
    seq: str,
    nt2complement: dict,
) -> str:
    """Reverse complement multi-nucleotide sequence

    Args:
        seq (str): sequence
        nt2complement (dict): nucleotide to complement

    Returns:
        str: sequence
    """
    complement = []
    for s in list(seq):
        for ss in nt2complement:
            if ss == s:
                #                 print(nt2complement[s],s)
                complement.append(nt2complement[s])
                break
    return "".join(complement[::-1])


def reverse_complement_multintseqreg(
    seq: str,
    multint2regcomplement: dict,
    nt2complement: dict,
) -> str:
    """Reverse complement multi-nucleotide regex patterns

    Args:
        seq (str): _description_
        multint2regcomplement (dict): mapping.
        nt2complement (dict): nucleotide to complement

    Returns:
        str: regex pattern
    """
    complement = []
    for s in list(seq):
        if s in multint2regcomplement.keys():
            for ss in multint2regcomplement:
                if ss == s:
                    #                 print(nt2complement[s],s)
                    complement.append(multint2regcomplement[s])
                    break
        elif s in nt2complement.keys():
            for ss in nt2complement:
                if ss == s:
                    complement.append(nt2complement[s])
                    break
        else:
            logging.error(f"odd character {s} in seq {seq}")

    return "".join(complement[::-1])


def hamming_distance(
    s1: str,
    s2: str,
) -> int:
    """Return the Hamming distance between equal-length sequences

    Args:
        s1 (str): sequence #1
        s2 (str): sequence #2

    Raises:
        ValueError: Undefined for sequences of unequal length

    Returns:
        int: distance.
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1.upper(), s2.upper()))


def align(
    q: str,  # query
    s: str,
    test: bool = False,
    psm: float = 2,
    pmm: float = 0.5,
    pgo: float = -3,
    pge: float = -1,
) -> str:
    """Creates pairwise local alignment between seqeunces.

    Args:
        q (str): query
        s (str): subject
        test (bool, optional): test-mode. Defaults to False.

    Returns:
        str: alignment with symbols.

    Notes:
        REF: http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
        The match parameters are:

        CODE  DESCRIPTION
        x     No parameters. Identical characters have score of 1, otherwise 0.
        m     A match score is the score of identical chars, otherwise mismatch
            score.
        d     A dictionary returns the score of any pair of characters.
        c     A callback function returns scores.
        The gap penalty parameters are:

        CODE  DESCRIPTION
        x     No gap penalties.
        s     Same open and extend gap penalties for both sequences.
        d     The sequences have different open and extend gap penalties.
        c     A callback function returns the gap penalties.
    """

    import itertools
    import operator

    from Bio import pairwise2

    def format_alignment(a):
        alignstr = pairwise2.format_alignment(*a)
        alignsymb = alignstr.split("\n")[1]
        return alignsymb.replace(" ", "-")

    assert len(q) == len(s), (len(q), len(s))
    q, s = q.upper(), s.upper()
    aligns = {}
    for i, (q_, s_) in enumerate(
        itertools.product(
            [q, str2seq(q).reverse_complement()], [s, str2seq(s).reverse_complement()]
        )
    ):
        if any([p is None for p in [psm, pmm, pgo, pge]]):
            alignments = pairwise2.align.localxx(q_, s_)
        else:
            alignments = pairwise2.align.localms(q_, s_, psm, pmm, pgo, pge)
        ## sort
        a = sorted(alignments, key=operator.itemgetter(2), reverse=True)[
            0
        ]  # descending order as opposed to stringent settings
        out = format_alignment(a)  # symbols
        if out == ("|" * len(a)):
            # perfect alignment
            return out
        else:
            aligns[i] = a
    ## no perfect alignment found
    a = sorted(list(aligns.values()), key=operator.itemgetter(2), reverse=True)[0]
    for i, a_ in aligns.items():
        if a_.score == a.score:
            out = format_alignment(a)  # symbols
            if i < 2:
                return out
            else:
                ## because query was rev com.d
                return out[::-1]


def get_orep(seq: str) -> int:
    """
    Get the overrepresentation
    """
    seq = seq.upper()
    return pd.Series({k: seq.count(k) for k in set(list(seq))}).max()


def get_polyt_length(s: str) -> int:
    """
    Counts the length of the longest polyT stretch (RNA pol3 terminator) in sequence

    :param s: sequence in string format
    """
    if isinstance(s, str):
        from itertools import groupby

        groups = groupby(s)
        result = [(label, sum(1 for _ in group)) for label, group in groups]
        df = pd.DataFrame(result, columns=["nt", "count"])
        df = df.loc[(df["nt"] == "T"), :]
        if len(df) > 0:
            return df.sort_values(by="count", ascending=False).iloc[0, 1]
        else:
            return np.nan
    else:
        return np.nan


## annots
def get_annots_installed() -> pd.DataFrame:
    """Get a list of annotations installed.

    Returns:
        pd.DataFrame: output.
    """
    import pandas as pd
    from pyensembl.shell import collect_all_installed_ensembl_releases

    genomes = collect_all_installed_ensembl_releases()
    return pd.DataFrame(
        [
            {
                "species": [g.species.latin_name] + g.species.synonyms,
                "assembly": g.reference_name,
                "source": g.annotation_name,
                "release": g.release,
            }
            for g in genomes
        ]
    ).explode("species")


def get_annots(
    species_name: str = None,
    release: int = None,
    gtf_path: str = None,
    transcript_path: str = None,
    protein_path: str = None,
    reference_name: str = "assembly",
    annotation_name: str = "source",
    verbose: bool = False,
    **kws_Genome,
):
    """Get pyensembl annotation instance

    Args:
        species_name (str, optional): species name. Defaults to None.
        release (int, optional): release number. Defaults to None.
        gtf_path (str, optional): GTF path. Defaults to None.
        transcript_path (str, optional): transcripts path. Defaults to None.
        protein_path (str, optional): protein path. Defaults to None.
        reference_name (str, optional): reference name. Defaults to 'assembly'.
        annotation_name (str, optional): annotation name. Defaults to 'source'.
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        pyensembl annotation instance
    """
    _loglevel = logging.getLevelName(logging.root.getEffectiveLevel())
    if species_name is None:
        from pyensembl import Genome

        annots = Genome(
            reference_name=reference_name,
            annotation_name=annotation_name,
            gtf_path_or_url=gtf_path,
            transcript_fasta_paths_or_urls=transcript_path,
            protein_fasta_paths_or_urls=protein_path,
            **kws_Genome,
        )
        # parse GTF and construct database of genomic features
        annots.index()
    else:
        if (
            get_annots_installed()
            .query(expr=f"`species`=='{species_name}' and `release`=={release}")
            .shape[0]
            == 0
        ):
            from beditor.lib.io import download_annots

            logging.getLogger().disabled = not verbose
            download_annots(
                species_name=species_name,
                release=release,
            )
            logging.getLogger().disabled = False
        from pyensembl import EnsemblRelease
        from pyensembl.species import find_species_by_name

        annots = EnsemblRelease(
            release=release,
            species=find_species_by_name(species_name),
        )
    logging.basicConfig(force=True, level=getattr(logging, _loglevel))
    return annots


def to_pid(
    annots,
    gid: str,
) -> str:
    """To protein ID

    Args:
        annots: pyensembl annotation instance
        gid (str): gene ID

    Returns:
        str: protein ID
    """
    ## get the longest transcript
    g = annots.gene_by_id(gid)
    assert g.is_protein_coding, g.is_protein_coding

    lens = {}
    for t in g.transcripts:
        if t.is_protein_coding and t.contains_start_codon and t.contains_stop_codon:
            lens[t.protein_id] = len(t.protein_sequence)
    if len(lens) != 0:
        return pd.Series(lens).sort_values().index.tolist()[-1]


## viz genome
def to_one_based_coordinates(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """To one based coordinates

    Args:
        df (pd.DataFrame): input table

    Returns:
        pd.DataFrame: output table.
    """
    for c in df.filter(like="start"):
        logging.warning(f"converting to 1-based:{c}")
        df = df.assign(
            **{
                c: lambda df: df[c] + 1,
            }
        )
    return df
