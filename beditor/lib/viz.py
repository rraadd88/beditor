#!usr/bin/python
"""Visualizations."""
import logging
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def to_igv(
    cfg: dict = None,
    gtf_path: str = None,
    genome_path: str = None,
    output_dir_path: str = None,
    threads: int = 1,
    output_ext: str = None,
    force: bool = False,
) -> str:
    """To IGV session file.

    Args:
        cfg (dict, optional): configuration of the run. Defaults to None.
        gtf_path (str, optional): path to the gtf file. Defaults to None.
        genome_path (str, optional): path to the genome file. Defaults to None.
        output_dir_path (str, optional): path to the output directory. Defaults to None.
        threads (int, optional): threads. Defaults to 1.
        output_ext (str, optional): extension of the output. Defaults to None.
        force (bool, optional): force. Defaults to False.

    Returns:
        str: path to the session file.
    """
    if output_ext is None:
        output_ext = "tsv"
    # ## to ensure relative paths do not interfere with igv display
    # from os.path import abspath
    # gtf_path=abspath(gtf_path)
    # genome_path=abspath(genome_path)
    # output_dir_path=abspath(output_dir_path)

    if cfg is None:
        from beditor.lib.io import to_viz_inputs

        cfg = to_viz_inputs(
            gtf_path=gtf_path,
            genome_path=genome_path,
            output_dir_path=output_dir_path,
            threads=threads,
            output_ext=output_ext,
            force=force,
        )
    cfg_base = {
        "reference": {
            "id": "Genome sequence",
            "name": "Genome sequence",
            "fastaURL": "genome_path",
            "indexURL": "{genome_path}.fai",
        },
        # "genome": "dm6",
        # "locus": "chr4:521,504-538,551",
        "tracks": [
            ## Targets (input mutations)
            {
                "name": "Targets",
                "url": "{output_dir_path}/06_viz/input.bed",
                "indexed": False,
            },
            ## sgRNAs
            {
                "name": "On-target sgRNAs",
                "url": "{output_dir_path}/06_viz/ontarget_sgrnas.bed",
                "indexed": False,
            },
            {
                "name": "Off-target sgRNAs",
                "url": "{output_dir_path}/06_viz/offtarget_sgrnas.bed",
                "indexed": False,
            },
            ## Annotations from gtf
            {
                "name": "Annotations",
                "type": "annotation",
                "format": "gtf",
                "displayMode": "expanded",
                "height": 200,
                "url": "{output_dir_path}/06_viz/ann_sorted.gtf.gz",
                "indexURL": "{output_dir_path}/06_viz/ann_sorted.gtf.gz.tbi",
                "colorBy": "biotype",
                # 'colorTable': {
                # },
            },
            ## sgRNAs
            {
                "name": "sgRNA alignments",
                "url": "{output_dir_path}/04_offtargets/alignment.bam",
                "indexURL": "{output_dir_path}/04_offtargets/alignment.bam.bai",
                "format": "bam",
                "type": "alignment",
                "height": 200,
            },
        ],
    }
    ## filll with inputs
    cfg_base["reference"] = {**cfg_base["reference"], **cfg["reference"]}
    if "locus" in cfg:
        cfg_base["locus"] = cfg["locus"]
    for i, d in enumerate(cfg_base["tracks"]):
        cfg_base["tracks"][i] = {
            **cfg_base["tracks"][i],
            **cfg[cfg_base["tracks"][i]["name"]],
        }

    outp1 = f"{output_dir_path}/06_viz/igv_session_paths.json"
    from roux.lib.io import to_dict

    to_dict(cfg_base, outp1)

    if Path(gtf_path).name.startswith("Homo_sapiens.GRCh38"):
        logging.warning("using cloud-based files for hg38.")
        # e.g. 'Homo_sapiens.GRCh38.110.gtf.gz.sorted.gtf.gz'
        ## for faster and robust visualization
        cfg_base["reference"] = {
            "id": "hg38",
            "name": "Human (GRCh38/hg38)",
            "fastaURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg38/hg38.fa",
            "indexURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/seq/hg38/hg38.fa.fai",
            "cytobandURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/hg38/cytoBandIdeo.txt.gz",
            "aliasURL": "https://igv-genepattern-org.s3.amazonaws.com/genomes/hg38/hg38_alias.tab",
            "chromSizesURL": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
            "twoBitURL": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit",
            "chromosomeOrder": [
                "chr1",
                "chr2",
                "chr3",
                "chr4",
                "chr5",
                "chr6",
                "chr7",
                "chr8",
                "chr9",
                "chr10",
                "chr11",
                "chr12",
                "chr13",
                "chr14",
                "chr15",
                "chr16",
                "chr17",
                "chr18",
                "chr19",
                "chr20",
                "chr21",
                "chr22",
                "chrX",
                "chrY",
            ],
        }

    from beditor.lib.io import to_session_path

    return to_session_path(
        p=cfg_base, outp=f"{output_dir_path}/06_viz/igv_session_urls.json"
    )


def get_nt_composition(
    seqs: list,
) -> pd.DataFrame:
    """Get nt composition.

    Args:
        seqs (list): list of sequences

    Returns:
        pd.DataFrame: table with the frequencies of the nucleotides.
    """
    import pandas as pd
    from Bio import motifs
    from beditor.lib.utils import str2seq

    instances = [str2seq(s) for s in seqs]
    motif = motifs.create(instances)
    d = pd.DataFrame(motif.counts)
    d.index = range(1, d.shape[0] + 1)
    return d


def plot_ntcompos(
    seqs: list,
    pam_pos: str,
    pam_len: int,
    window: list = None,
    ax: plt.Axes = None,
    color_pam: str = "lime",
    color_window: str = "gold",
) -> plt.Axes:
    """Plot nucleotide composition

    Args:
        seqs (list): list of sequences.
        pam_pos (str): PAM position.
        pam_len (int): PAM length.
        window (list, optional): activity window bounds. Defaults to None.
        ax (plt.Axes, optional): subplot. Defaults to None.
        color_pam (str, optional): color of the PAM. Defaults to 'lime'.
        color_window (str, optional): color of the wnindow. Defaults to 'gold'.

    Returns:
        plt.Axes: subplot
    """
    df1 = get_nt_composition(seqs)
    ax = df1.loc[:, list("ATGC")].plot(
        grid=True,
        # xticks=dp.index,
        ax=ax,
        lw=3,
        color=["r", "b", "orange", "g"],
        label=None,
        # cmap='Spectral'
        legend=False,
    )

    if pam_pos == "up":
        df1["PAM"] = ([True] * pam_len) + ([False] * (len(df1) - pam_len))
        df1["pos"] = list(range(-1 * (pam_len), 0, 1)) + list(
            range(
                1,
                (len(df1) - pam_len) + 1,
            )
        )
    elif pam_pos == "down":
        df1["PAM"] = ([False] * (len(df1) - pam_len)) + ([True] * pam_len)
        df1["pos"] = (
            list(
                range(
                    1,
                    (len(df1) - pam_len) + 1,
                )
            )[::-1]
            + list(range(-1 * (pam_len), 0, 1))[::-1]
        )

    start, end = (
        df1.query("`PAM`==True").index.min(),
        df1.query("`PAM`==True").index.max(),
    )
    ax.axvspan(
        start - 0.5,
        end + 0.5,
        lw=4,
        color=color_pam,
        zorder=-1,
        label="PAM",
    )

    if window is not None:
        guidepam_seq_len = len(df1)
        if pam_pos == "up":
            _wstart, _wend = pam_len + window[0], pam_len + window[1]
        elif pam_pos == "down":
            _wstart, _wend = (
                (guidepam_seq_len - pam_len) - (window[1] - 1),
                (guidepam_seq_len - pam_len) - (window[0] - 1),
            )
        ax.axvspan(
            _wstart - 0.5,
            _wend + 0.5,
            lw=4,
            color=color_window,
            zorder=-1,
            label="Activity window",
        )
    ax.set(
        xticks=df1.index,
        xlim=[df1.index.min(), df1.index.max()],
        xticklabels=df1["pos"].tolist(),
        xlabel="Position",
        ylabel="Nucleotide\ncomposition",
    )
    plt.grid(linewidth=0.5, color="lightgray", alpha=0.8)
    ax.legend(loc=2, bbox_to_anchor=[1, 1])
    return ax


def plot_ontarget(
    guide_loc: str,
    # guide_strand: str,
    pam_pos: str,
    pam_len: int,
    guidepam_seq: str,
    window: list = None,
    show_title: bool = False,
    figsize: list = [10, 2],
    verbose: bool = False,
    kws_sg: dict = {},
) -> plt.Axes:
    """plot_ontarget _summary_

    Args:
        guide_loc (str): sgRNA locus
        pam_pos (str): PAM position
        pam_len (int): PAM length
        guidepam_seq (str): sgRNA and PAM sequence
        window (list, optional): activity window bounds. Defaults to None.
        show_title (bool, optional): show the title. Defaults to False.
        figsize (list, optional): figure size. Defaults to [10,2].
        verbose (bool, optional): verbose. Defaults to False.
        kws_sg (dict, optional): keyword arguments to plot the sgRNA. Defaults to {}.

    Returns:
        plt.Axes: subplot

    TODOs:
        1. convert to 1-based coordinates
        2. features from the GTF file
    """
    # guide_start=int(guide_loc.split(':')[1].split('-')[0])
    # guide_end=int(guide_loc.split(':')[1].split('-')[1].split('(')[0])
    # chrom=guide_loc.split(':')[0]
    from beditor.lib.utils import parse_locus

    chrom, guide_start, guide_end, guide_strand = parse_locus(guide_loc)

    if window is not None:
        wlen = window[1] - window[0]
    if (pam_pos == "up" and guide_strand == "+") or (
        pam_pos == "down" and guide_strand == "-"
    ):
        _guide_start, _guide_end = pam_len, len(guidepam_seq)
        start, end = guide_start - pam_len, guide_end + 1
        if window is not None:
            _wstart, _wend = pam_len + window[0] - 1, pam_len + window[1] - 1
    elif (pam_pos == "up" and guide_strand == "-") or (
        pam_pos == "down" and guide_strand == "+"
    ):
        _guide_start, _guide_end = 0, len(guidepam_seq) - pam_len
        start, end = guide_start, guide_end + pam_len + 1
        if window is not None:
            _wstart, _wend = (
                (len(guidepam_seq) - pam_len) - window[1],
                (len(guidepam_seq) - pam_len) - window[0],
            )

    if verbose:
        print(f"guide_strand={guide_strand}")
        print(_guide_start, _guide_end)
        print(start, end, end - start)

    from dna_features_viewer import GraphicFeature, GraphicRecord

    features = [
        GraphicFeature(
            start=_guide_start,
            end=_guide_end,
            strand=1 if guide_strand == "+" else -1,
            # color="#ffd700",
            label="sgRNA",
            **kws_sg,
        ),
    ]

    from beditor.lib.utils import str2seq

    record = GraphicRecord(
        sequence=guidepam_seq
        if guide_strand == "+"
        else str(str2seq(guidepam_seq).complement())[::-1],
        sequence_length=len(guidepam_seq),
        features=features,
        ticks_resolution=1,
        first_index=0,
        feature_level_height=1,
    )
    # cropped_record = record.crop((guide_start, guide_end))
    ax, _ = record.plot(
        plot_sequence=True,
        x_lim=None,
        figure_width=figsize[0],
        figure_height=figsize[1],
    )
    _ = ax.set_xticklabels(range(start, end), rotation=90)
    if window is not None:
        ax.axvspan(xmin=_wstart - 0.5, xmax=_wend + 0.5, color="gold", zorder=-1)
    if show_title:
        ax.set_title(f"{chrom}:{start}-{end}", loc="left")
    return ax


def get_plot_inputs(
    df2: pd.DataFrame,
) -> list:
    """Get plot inputs.

    Args:
        df2 (pd.DataFrame): table.

    Returns:
        list: list of tables.
    """
    from beditor.lib.utils import cols_muts

    dfs = {}
    dfs["Off-target alignments"] = (
        df2.groupby("guide sequence")["offtargets"].sum().value_counts()
    )
    dfs["Off-target alignments"].name = "Off-target alignments"

    dfs["Genes targeted by sgRNA"] = (
        df2.groupby("aligned genes count")["guide sequence"].nunique().sort_index()
    )
    dfs["Genes targeted by sgRNA"].name = "Genes targeted by sgRNA"

    dfs["PolyT stretch"] = (
        df2.groupby("polyT stretch length")["guide sequence"].nunique().sort_index()
    )
    dfs["PolyT stretch"].name = "PolyT stretch"
    dfs["PolyT stretch"].index = dfs["PolyT stretch"].index.astype(int)

    dfs["beditor score"] = df2.drop_duplicates(
        subset=["guide+PAM sequence"] + cols_muts
    )["score"]
    return dfs


def plot_library_stats(
    dfs: list,
    palette: dict = {True: "b", False: "lightgray"},
    cutoffs: dict = None,
    not_be: bool = True,
    dbug: bool = False,
    figsize: list = [10, 2.5],
) -> list:
    """Plot library stats

    Args:
        dfs (list): list of tables.
        palette (_type_, optional): color palette. Defaults to {True:'b',False:'lightgray'}.
        cutoffs (dict, optional): cutoffs to be applied. Defaults to None.
        not_be (bool, optional): not a base editor. Defaults to True.
        dbug (bool, optional): debug mode. Defaults to False.
        figsize (list, optional): figure size. Defaults to [10,2.5].

    Returns:
        list: list of subplots.
    """
    fig, axs = plt.subplots(1, 4, figsize=figsize, sharey=True)
    
    k = "beditor score"
    ax = axs[0]
    ax = (dfs[k] * 100).hist(color=palette[True], ax=ax)
    if cutoffs is not None:
        ax.axvspan(
            0,
            cutoffs[k][0] * 100,
            color=palette[False],
            zorder=0,
        )
    ax.grid(False)
    ax.set(
        xlim=[0, 100],
        xlabel=k + ("*" if not_be else ""),
        ylabel="sgRNAs",
    )

    for i, k in enumerate(
        [
            "Off-target alignments",
            "Genes targeted by sgRNA",
            "PolyT stretch",
        ]
    ):
        # k='Off-target alignments'
        data = pd.DataFrame(dfs[k])
        if cutoffs is not None:
            data["color"] = [
                palette[b]
                for b in (
                    pd.Series(data.index >= cutoffs[k][0])
                    & pd.Series(data.index <= cutoffs[k][1])
                )
            ]
        else:
            data["color"] = [palette[True] for b in range(len(data))]

        ax = data[k].plot.bar(
            color=data["color"],
            ax=axs[i+1],
        )
        ax.set(xlabel=k)

    if dbug:
        print(min(data.index.tolist()), max(data.index.tolist()))

    plt.tight_layout()
    return axs
