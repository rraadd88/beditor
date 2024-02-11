#!usr/bin/python
"""Input/Output"""
import logging
from pathlib import Path

import pandas as pd
from roux.lib.io import read_dict, read_table, to_dict
from roux.lib.sys import makedirs
import roux.lib.df as rd #noqa

from beditor.lib.utils import format_coords


def download_annots(
    species_name: str,
    release: int,
) -> bool:
    """Download annotations using pyensembl

    Args:
        species_name (str): species name
        release (int): release number

    Returns:
        bool: whether annotation is downloaded or not
    """

    # species_name,assembly,release
    class DictAsObject:
        def __init__(self, dictionary):
            self.__dict__.update(dictionary)

    # Example usage:
    args = DictAsObject(
        dict(species=[species_name], release=[release], custom_mirror=None)
    )
    from pyensembl.shell import all_combinations_of_ensembl_genomes

    genome = all_combinations_of_ensembl_genomes(args)[0]
    logging.info(
        "Downloading the genome files if using the genome for the first time. Depending on the size of the genome, this step can take a few minutes to finish."
    )
    genome.download(overwrite=False)
    genome.index(overwrite=False)
    return True


## local (cache) directory
from os.path import exists, join, split

import datacache


def cache_subdirectory(
    reference_name: str = None,
    annotation_name: str = None,
    annotation_version: int = None,
    CACHE_BASE_SUBDIR: str = "beditor",
) -> str:
    """
    Which cache subdirectory to use for a given annotation database
    over a particular reference. All arguments can be omitted to just get
    the base subdirectory for all pyensembl cached datasets.

    Args:
        reference_name (str, optional): reference name. Defaults to None.
        annotation_name (str, optional): annotation name. Defaults to None.
        annotation_version (int, optional): annotation version. Defaults to None.
        CACHE_BASE_SUBDIR (str, optional): cache path. Defaults to 'beditor'.

    Returns:
        str: output path
    """
    if reference_name is None:
        reference_name = ""
    if annotation_name is None:
        annotation_name = ""
    if annotation_version is None:
        annotation_version = ""
    reference_dir = join(CACHE_BASE_SUBDIR, reference_name)
    annotation_dir = "%s%s" % (annotation_name, annotation_version)
    return join(reference_dir, annotation_dir)


def cached_path(
    path_or_url: str,
    cache_directory_path: str,
    # decompress_on_download=False,
):
    """
    When downloading remote files, the default behavior is to name local
    files the same as their remote counterparts.
    """

    def is_url_format(path_or_url):
        """
        Is the given string a URL?

        Parameters
        ----------
        path_or_url : str

        Returns
        -------
        bool
        """
        if path_or_url is None or path_or_url == "":
            raise ValueError("Expected non-empty string for path_or_url")
        return "://" in path_or_url

    if path_or_url is None or path_or_url == "":
        raise ValueError("Expected non-empty string for path_or_url")
    remote_filename = split(path_or_url)[1]
    if is_url_format(path_or_url):
        # passing `decompress=False` since there is logic below
        # for stripping decompression extensions for both local
        # and remote files
        local_filename = datacache.build_local_filename(
            download_url=path_or_url, filename=remote_filename, decompress=False
        )
    else:
        local_filename = remote_filename

    # if we expect the download function to decompress this file then
    # we should use its name without the compression extension
    # if decompress_on_download:
    #     local_filename = _remove_compression_suffix_if_present(local_filename)

    if len(local_filename) == 0:
        raise ValueError("Can't determine local filename for %s" % (path_or_url,))

    return join(cache_directory_path, local_filename)


def to_downloaded_cached_path(
    url: str,
    annots=None,
    # or
    reference_name: str = None,
    annotation_name: str = "ensembl",
    ensembl_release: str = None,
    CACHE_BASE_SUBDIR: str = "pyensembl",
) -> str:
    """To downloaded cached path

    Args:
        url (str): URL
        annots (optional): pyensembl annotation. Defaults to None.
        reference_name (str, optional): reference name. Defaults to None.
        annotation_name (str, optional): annotation name. Defaults to 'ensembl'.
        ensembl_release (str, optional): ensembl release. Defaults to None.
        CACHE_BASE_SUBDIR (str, optional): cache path. Defaults to 'pyensembl'.

    Returns:
        str: output path
    """
    if annots is not None:
        reference_name = annots.reference_name
        annotation_name = annots.annotation_name
        ensembl_release = annots.release
    cache_subdirectory_path = cache_subdirectory(
        reference_name=reference_name,
        annotation_name=annotation_name,
        annotation_version=ensembl_release,
        CACHE_BASE_SUBDIR=CACHE_BASE_SUBDIR,
    )
    # If `CACHE_DIR_ENV_KEY` is set, the cache will be saved there
    cache_directory_path = datacache.get_data_dir(
        subdir=cache_subdirectory_path,
        envkey="PYENSEMBL_CACHE_DIR",
    )

    downloaded_cached_path = cached_path(
        path_or_url=url,
        cache_directory_path=cache_directory_path,
    )
    return downloaded_cached_path


def download_genome(
    species: str,
    ensembl_release: int,
    force: bool = False,
    verbose: bool = False,
) -> str:
    """Download genome

    Args:
        species (str): species name
        ensembl_release (int): release
        force (bool, optional): force. Defaults to False.
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        str: output path
    """
    ## sci-names
    from pyensembl.ensembl_url_templates import normalize_release_properties

    ensembl_release, species, reference_name = normalize_release_properties(
        ensembl_release, species
    )
    if verbose:
        logging.info(f"{ensembl_release}, {species}, {reference_name}")

    ## url
    from pyensembl.ensembl_url_templates import make_fasta_url

    # .dna.alt.:         These files contain alternate sequences or alternative haplotypes for regions in the genome.
    # .dna.toplevel.:    These files include the primary assembly sequences for the entire genome.
    # .dna_rm.alt.:      Similar to .dna.alt., but specifically representing repeat-masked alternate sequences.
    # .dna_sm.alt.:      These files contain soft-masked alternate sequences, where repetitive regions are marked in lowercase letters.
    # .dna_sm.toplevel.: Similar to .dna.toplevel., but with soft-masked sequences for repetitive regions.
    url = make_fasta_url(ensembl_release, species, "dna").replace(
        ".dna.all.", ".dna.toplevel."
    )
    if verbose:
        logging.info(url)
    downloaded_cached_path = to_downloaded_cached_path(
        url=url,
        annots=None,
        # or
        reference_name=reference_name,
        # annotation_name=annotation_name,
        ensembl_release=ensembl_release,
        CACHE_BASE_SUBDIR="beditor",
    )
    if verbose:
        logging.info(downloaded_cached_path)
    ## downnload
    if not exists(downloaded_cached_path) or force:
        if verbose:
            print("Fetching %s from URL %s", downloaded_cached_path, url)
            print(
                "Note: Depending on the size of the genome, this step can take a while."
            )
            print(
                "Note: Downloading is necessary only for the first-time use of a genome."
            )
            print(
                "Note: At later use of the genome, beditor will use a downloaded local file."
            )
        datacache.ensure_dir(str(Path(downloaded_cached_path).parent))
        datacache.download._download_and_decompress_if_necessary(
            full_path=downloaded_cached_path, download_url=url, timeout=3600
        )
        assert exists(
            downloaded_cached_path
        ), f"Error in downloading file {downloaded_cached_path} from URL {url}"
    return downloaded_cached_path


## common functions
def read_genome(
    genome_path: str,
    fast=True,
):
    """Read genome

    Args:
        genome_path (str): genome path
        fast (bool, optional): fast mode. Defaults to True.
    """
    from Bio import SeqIO

    if fast:
        return SeqIO.index(genome_path, "fasta")
    else:
        import gzip

        return SeqIO.parse(
            open(genome_path, "r")
            if not genome_path.endswith(".gz")
            else gzip.open(genome_path, "rt"),
            "fasta",
        )


def to_fasta(
    sequences: dict,
    output_path: str,
    molecule_type: str,
    force: bool = True,
    **kws_SeqRecord,
) -> str:
    """Save fasta file.

    Args:
        sequences (dict): dictionary mapping the sequence name to the sequence.
        output_path (str): path of the fasta file.
        force (bool): overwrite if file exists.

    Returns:
        output_path (str): path of the fasta file
    """
    from Bio import Seq, SeqIO, SeqRecord

    assert len(sequences) != 0
    from roux.lib.sys import exists, makedirs

    if exists(output_path) and not force:
        logging.warning("file exists.")
        return
    makedirs(output_path)
    assert molecule_type in ["Protein", "RNA", "DNA"], molecule_type
    seqs = (
        SeqRecord.SeqRecord(
            Seq.Seq(sequences[k]),
            id=k,
            annotations=dict(molecule_type=molecule_type),
            **kws_SeqRecord,
        )
        for k in sequences
    )
    SeqIO.write(seqs, output_path, "fasta")
    return output_path


def to_2bit(
    genome_path: str,
    src_path: str = None,  # path to bin
    force: bool = False,
    verbose: bool = False,
) -> str:
    """To 2bit

    Args:
        genome_path (str): genome path
        src_path (str, optional): source path. Defaults to None.
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        str: output path
    """
    from beditor.lib.utils import runbashcmd

    if genome_path.endswith("gz"):
        outp = Path(genome_path).with_suffix("")
        if not Path(outp).exists() or force:
            com = f"gunzip {genome_path}"
            if verbose:
                print(com)
            runbashcmd(com)
        genome_path = outp
    genome_2bit_path = Path(genome_path).with_suffix(".2bit")
    if not Path(genome_2bit_path).exists() or force:
        if src_path is None:
#             import site

#             src_path = (
#                 Path(site.getsitepackages()[0]).parent.parent.parent.as_posix() + "/bin/faToTwoBit"
#             )
#             logging.warning(f"inferred src_path={src_path}")
            src_path="micromamba run -n twobit faToTwoBit"
        com = f"{src_path} {genome_path} {genome_2bit_path}"
        if verbose:
            print(com)
        runbashcmd(com)

        com = f"{src_path} {genome_2bit_path} stdout | sort -k2rn > {Path(genome_2bit_path).with_suffix('.chrom.sizes')}"
        if verbose:
            print(com)
        runbashcmd(com)
    return genome_2bit_path

def to_fasta_index(
    genome_path: str,
    bgzip: bool=False,
    bgzip_path: str = None,
    threads: int = 1,
    verbose: bool = True,
    force: bool = False,
    indexed: bool = False,
) -> str:
    """To fasta index

    Args:
        genome_path (str): genome path
        bgzip_path (str, optional): bgzip path. Defaults to None.
        threads (int, optional): threads. Defaults to 1.
        verbose (bool, optional): verbose. Defaults to True.
        force (bool, optional): force. Defaults to False.
        indexed (bool, optional): indexed or not. Defaults to False.

    Returns:
        str: output path
    """
    from pyfaidx import Faidx

    if Path(genome_path).suffix == ".gz":
        ## check if gzip'ed or bgzip'ed
        try:
            out = Faidx(genome_path, rebuild=force)
            genome_fmt = "bgzip"
            indexed = True
        except:
            genome_fmt = "gzip"
            logging.warning(f"genome_fmt={genome_fmt}")
            indexed = False
    else:
        genome_fmt=genome_path.split('.')[-1]
    from beditor.lib.utils import runbashcmd
    if genome_fmt == "gzip":
        if not Path(genome_path[:-3]).exists():
            runbashcmd(f"gunzip -k {genome_path}")
        genome_path=genome_path[:-3]
        if bgzip:
            if bgzip_path is None:
                import site
                bgzip_path = site.getsitepackages()[0] + "/pysam/include/htslib/bgzip"
            if not Path(bgzip_path).exists() and Path(bgzip_path).parent.exists():
                if verbose:
                    print("Compiling bgzip")
                cmd = f"cd {Path(bgzip_path).parent};./configure --disable-bz2 --disable-lzma;make"
                if verbose:
                    print(cmd)
                runbashcmd(cmd)
            if verbose:
                print(f"Using bgzip compiled at {bgzip_path}")

            # outp = genome_path.replace(".gz", ".bgz")
            # if not Path(outp).exists():
            if verbose:
                print("Converting gzip to bgzip")
                print(
                    "Note: Depending on the size of the genome, this step can take a while."
                )
            # cmd=f"zcat {genome_path} | {bgzip_path} -c --threads {threads} > {outp}"
            cmd = f"{bgzip_path} -i -f --threads {threads} {genome_path}"
            if verbose:
                print(cmd)
            runbashcmd(cmd)
            genome_path=genome_path+'.gzip'
    print(f"genome_path={genome_path}")
    if not indexed:
        out = Faidx(genome_path, rebuild=force)
    return genome_path


def to_bed(
    df: pd.DataFrame,
    outp: str,
    cols: list = ["chrom", "start", "end", "locus", "score", "strand"],
) -> str:
    """To bed path

    Args:
        df (pd.DataFrame): input table
        outp (str): output path
        cols (list, optional): columns. Defaults to ['chrom','start','end','locus','score','strand'].

    Returns:
        str: output path
    """
    if "locus" not in df:
        from beditor.lib.utils import to_locus

        df = df.assign(
            locus=lambda df: df.apply(
                lambda x: to_locus("chrom", "start", "end", "strand", x=x), axis=1
            )
        )
    if "score" not in df:
        df = df.assign(score=0)  ## dummy score needed for UCSC-style bed
    else:
        df = df.assign(
            score=lambda df: df["score"].fillna(0)
        )  ## dummy score needed for UCSC-style bed
    makedirs(outp)
    (
        df.loc[:, cols]
        .rd.assert_no_na()
        .sort_values(cols[:3])
        .drop_duplicates()
        .to_csv(
            outp,
            index=False,
            header=False,
            sep="\t",
        )
    )
    return outp


def read_bed(
    p: str,
    cols: list = ["chrom", "start", "end", "locus", "score", "strand"],
) -> pd.DataFrame:
    """Read bed file

    Args:
        p (str): path
        cols (list, optional): columns. Defaults to ['chrom','start','end','locus','score','strand'].

    Returns:
        pd.DataFrame: output table
    """
    return pd.read_csv(p, sep="\t", header=None, names=cols)


def to_viz_inputs(
    gtf_path: str,
    genome_path: str,
    output_dir_path: str,
    output_ext: str = "tsv",
    threads: int = 1,
    force: bool = False,
) -> dict:
    """To viz inputs for the IGV

    Args:
        gtf_path (str): GTF path
        genome_path (str): genome path
        output_dir_path (str): output directory path
        output_ext (str, optional): output extension. Defaults to 'tsv'.
        threads (int, optional): threads. Defaults to 1.
        force (bool, optional): force. Defaults to False.

    Returns:
        dict: configuration
    """
    import pysam

    makedirs(f"{output_dir_path}/06_viz/")

    cfg = {
        "reference": {
            "fastaURL": genome_path,
            "indexURL": f"{genome_path}.fai",
        },
        # tracks
        "Targets": {
            "url": f"{output_dir_path}/01_sequences/targets.bed",
        },
        "On-target sgRNAs": {
            "url": f"{output_dir_path}/06_viz/ontarget_sgrnas.bed",
        },
        "Off-target sgRNAs": {
            "url": f"{output_dir_path}/06_viz/offtarget_sgrnas.bed",
        },
        "Annotations": {
            "url": f"{gtf_path}.sorted.gtf.gz",
            "indexURL": f"{gtf_path}.sorted.gtf.gz.tbi",
        },
        "sgRNA alignments": {
            "url": f"{output_dir_path}/04_offtargets/alignment.bam",
            "indexURL": f"{output_dir_path}/04_offtargets/alignment.bam.bai",
        },
    }
    ### Reference sequence
    if genome_path.endswith(".gz"):
        cfg["reference"]["compressedIndexURL"] = f"{genome_path}.gzi"

    if not Path(cfg["reference"]["indexURL"]).exists():
        print("indexing fasta")
        ## index fasta
        to_fasta_index(
            genome_path=genome_path,
            bgzip_path=None,
            dbug=True,
            threads=threads,
        )

    #### Targets
    # if not Path(f"{output_dir_path}/06_viz/input.bed").exists() or force:
    # print("track for Targets")
    # p=glob(f'{output_dir_path}/01_sequences.*')[0]
    # df_=(read_table(p)
    # .assign(
    #     start=lambda df: df.apply(lambda x: x['start']-1,axis=1), ## 1-based to 0-based
    #     # label=lambda df: df.apply(lambda x: f"{'#'+x['aa pos'] if 'aa pos' in x else str(x['start'])+':'+str(x['end'])}({x['sequence target']})",axis=1),
    #     label=lambda df: df['sequence target'],
    #     score=None,#lambda df: df['sequence target'],
    # )
    # )
    # x=df_.sort_values(['chrom','start']).iloc[0,:]
    # cfg['locus']=f"chr{x['chrom']}:{x['start']}-{x['end']}"
    # to_bed(df_,outp=cfg['Targets']['url'])
    # del df_,x
    with open(cfg["Targets"]["url"], "r") as f:
        first_line = f.readline()
    c, s, e = first_line.split("\t")[:3]
    cfg["locus"] = "{c}:{s}-{e}"

    ## sgRNA on and offtargets
    df3 = read_table(
        f"{output_dir_path}/04_offtargets/01_mappedby_alignments.{output_ext}"
    ).log("offtarget")
    # df3.head(1)
    df3 = df3.assign(
        score=None,
        label=lambda df: df.apply(
            lambda x: f"{x['guide sequence']} ({x['guide locus']})", axis=1
        ),
    )
    # %run ../beditor/lib/io.py
    for k, df in df3.groupby("offtarget"):
        if len(df) != 0:
            # outp=f"{output_dir_path}/06_viz/"+('on' if not k else 'off')+'target_sgrnas'+".bed"
            to_bed(
                df.loc[
                    :,
                    [
                        "aligned chrom",
                        "aligned start",
                        "aligned end",
                        "label",
                        "score",
                        "aligned strand",
                    ],
                ].rd.renameby_replace({"aligned ": ""}, ignore=True),
                cfg[("On" if not k else "Off") + "-target sgRNAs"]["url"],
            )
            # del outp

    ### Annotations
    _outp = f"{gtf_path}.sorted.gtf.gz.tbi"
    if not Path(_outp).exists() or force:
        _outp1 = f"{gtf_path}.sorted.gtf"
        if not Path(_outp).exists() or force:
            print("sorting Annotations")
            read_table(
                gtf_path,
                params=dict(sep="\t", header=None, comment="#", low_memory=False),
            ).sort_values([0, 3]).to_csv(
                _outp1,
                index=False,
                header=False,
                sep="\t",
            )
        if not Path(_outp).exists() or force:
            print("indexing Annotations")
            pysam.tabix_index(_outp1, preset="gff", keep_original=True, force=force)

    ### sgRNAs alignments
    if not Path(cfg["sgRNA alignments"]["indexURL"]).exists() or force:
        print("indexing sgRNAs alignments")
        pysam.sort(
            "-o",
            cfg["sgRNA alignments"]["url"],
            f"{output_dir_path}/04_offtargets/alignment.sam",
        )
        pysam.index(cfg["sgRNA alignments"]["url"])

    ## set symlinks withing the output folder
    from roux.lib.sys import create_symlink

    outp = f"{output_dir_path}/06_viz/{Path(cfg['reference']['fastaURL']).name}"
    if not Path(outp).exists() or force:
        cfg["reference"]["fastaURL"] = create_symlink(
            cfg["reference"]["fastaURL"],
            outp,
        )
    cfg["reference"]["fastaURL"] = outp

    outp = f"{output_dir_path}/06_viz/{Path(cfg['reference']['indexURL']).name}"
    if not Path(outp).exists() or force:
        cfg["reference"]["indexURL"] = create_symlink(
            cfg["reference"]["indexURL"],
            outp,
        )
    cfg["reference"]["indexURL"] = outp

    # ds=[]
    # for d in cfg['tracks']:
    #     if d['name']=='Annotations':
    outp = f"{output_dir_path}/06_viz/{Path(cfg['Annotations']['url']).name}"
    if not Path(outp).exists() or force:
        cfg["Annotations"]["url"] = create_symlink(
            cfg["Annotations"]["url"],
            outp,
        )
    cfg["Annotations"]["url"] = outp

    outp = f"{output_dir_path}/06_viz/{Path(cfg['Annotations']['indexURL']).name}"
    if not Path(outp).exists() or force:
        cfg["Annotations"]["indexURL"] = create_symlink(
            cfg["Annotations"]["indexURL"],
            outp,
        )
    cfg["Annotations"]["indexURL"] = outp
    #     ds.append(d)
    # cfg['tracks']=ds

    return cfg


def to_igv_path_prefix() -> str:
    """Get IGV path prefix

    Returns:
        str: URL
    """
    import subprocess

    r = subprocess.run(
        "jupyter server list",
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        url, dir_path = r.stdout.split("\n")[1].split(" :: ")
        if Path.cwd().as_posix().startswith("/mnt"):  # wsl for external storage
            dir_path = "/".join(Path.cwd().as_posix().split("/")[3:])
    except:
        print(r.stdout)
        url = "https://localhost:port/"
        dir_path = "path/to/wd"
    return f"{url}files/{dir_path}/"


def to_session_path(
    p: str,
    path_prefix: str = None,
    outp: str = None,
) -> str:
    """To session path

    Args:
        p (str): session configuration path
        path_prefix (str, optional): path prefix. Defaults to None.
        outp (str, optional): output path. Defaults to None.

    Returns:
        str: output path
    """
    if path_prefix is None:
        path_prefix = to_igv_path_prefix()
    if isinstance(p, str):
        cfg = read_dict(p)
    elif isinstance(p, dict):
        cfg = p.copy()
    for k in cfg:
        if isinstance(cfg[k], dict):
            for k1 in cfg[k]:
                if isinstance(cfg[k][k1], str):
                    if Path(cfg[k][k1]).exists():
                        cfg[k][k1] = path_prefix + cfg[k][k1]
        if isinstance(cfg[k], list):
            ds = []
            for d in cfg[k]:
                for k1 in d:
                    if isinstance(d[k1], str):
                        if Path(d[k1]).exists():
                            d[k1] = path_prefix + d[k1]
                ds.append(d)
            cfg[k] = ds
    # cfg
    if outp is None:
        import tempfile

        outp = tempfile.NamedTemporaryFile().name + ".json"
    return to_dict(cfg, outp)


# viz
def read_cytobands(
    cytobands_path: str,
    col_chrom: str = "chromosome",  # Todo change to chrom
    remove_prefix: str = "chr",  # remove chr
) -> pd.DataFrame:
    """Read cytobands

    Args:
        cytobands_path (str): path
        col_chrom (str, optional): column with contig. Defaults to 'chromosome'.

    Returns:
        pd.DataFrame: output table
    """
    from beditor.lib.utils import to_one_based_coordinates

    df02 = pd.read_csv(
        cytobands_path,
        sep="\t",
        compression="gzip" if cytobands_path.endswith(".gz") else None,
        header=None,
        names=[col_chrom, "start", "end", "cytoband", "cytoband type"],
    )

    ## clean
    df0 = df02.pipe(to_one_based_coordinates)
    # df0

    ## subset
    df1 = (
        df0.query("`cytoband type`=='acen'")
        .assign(
            **{
                "arm type": lambda df: df["cytoband"].str[0],
            }
        )
        .pivot(index=col_chrom, columns=["arm type"], values=["start", "end"])
        .swaplevel(axis=1)
        .rd.flatten_columns()
        .reset_index()
    )
    # df1.head(1)

    assert all(df1["p end"] + 1 == df1["q start"])
    df1 = df1.drop(["p end", "q start"], axis=1)
    assert all(df1["p start"] < df1["q end"])

    # ## lengths
    # df01=(
    #     pd.read_csv(genome_path+'.fai',sep='\t',names=['chrom','length','byte index','bases per line','bytes per line'])
    #     .loc[:,['chrom','length']]
    #     .astype({'chrom':str,'length': int})
    # )
    # ## centromeres
    # df2=(df1
    #     .log.merge(right=df01,
    #               on='chrom',
    #              how='inner',
    #              validate="1:1")
    #     )

    df4 = df1.log.merge(
        right=df0,
        how="inner",
        on=col_chrom,
        validate="1:m",
    )
    if remove_prefix is not None:
        df4 = df4.assign(
            **{
                col_chrom: lambda df: df[col_chrom].apply(
                    lambda x: x[(len(remove_prefix)) :]
                    if str(x).startswith(remove_prefix)
                    else x
                ),
            },
        )

    return (
        df4.assign(
            **{
                "arm": lambda df: df["cytoband"].str[0],
                f"{col_chrom} arm": lambda df: "chr" + df[col_chrom] + df["arm"],
            },
        )
        # .rename(
        #     columns={
        #         'chrom':col_chrom,
        #         # "chrom arm":f"{col_chrom} arm",
        #     },
        #     )
    )


def to_output(
    inputs: pd.DataFrame,
    guides: pd.DataFrame,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    """To output table

    Args:
        inputs (pd.DataFrame): inputs
        guides (pd.DataFrame): guides
        scores (pd.DataFrame): scores

    Returns:
        pd.DataFrame: output table
    """
    ### Map back to the mutations
    # types={
    #     'chrom':str,
    #     'start':int,
    #     'end':int,
    #     'strand':str,
    # }
    cols_merge_on = ["target location"] + (["mutation"] if "mutation" in inputs else [])
    cols_right = list(set(guides.columns.tolist()) - set(inputs.columns.tolist()))
    df0 = (
        inputs
        .pipe(format_coords)
        .log.merge(
            right=guides.loc[:, cols_merge_on + cols_right],
            how="left",
            on=cols_merge_on,
        )
        .log.drop_duplicates()
        )
    cols_merge_on = [
        # 'PAM position',
        "guide locus",
        "guide sequence",
    ]
    cols_right = list(set(scores.columns.tolist()) - set(df0.columns.tolist()))
    df1 = (df0
        .log.merge(
            right=scores.loc[:, cols_merge_on + cols_right],
            how="inner",
            on=cols_merge_on,
            suffixes=["", " aligned"],
            )
        .log.drop_duplicates()
          )
    ### Post-processing
    from beditor.lib.utils import get_polyt_length

    return df1.assign(
        **{
            "aligned genes count": lambda df: df["aligned genes"].apply(
                lambda x: x.count(";") + 1
            ),
            "polyT stretch length": lambda df: df["guide sequence"]
            .apply(lambda x: get_polyt_length(x))
            .fillna(0),
        }
    )
