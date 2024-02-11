#!usr/bin/python
"""Mutation co-ordinates using pyensembl"""

## import required libraries
import logging
from pathlib import Path

import pandas as pd
from roux.lib.io import read_table, to_table
import roux.lib.df as rd #noqa

def get_protein_cds_coords(
    annots,
    protein_id: str,
) -> pd.DataFrame:
    """Get protein CDS coordinates

    Args:
        annots: pyensembl annotations
        protein_id (str): protein ID

    Returns:
        pd.DataFrame: output table
    """
    t = annots.transcript_by_protein_id(protein_id)
    # g=annots.gene_by_protein_id(protein_id)
    if not (t.is_protein_coding and t.contains_start_codon and t.contains_stop_codon):
        logging.error(
            f"Excluded protein id {protein_id} because either it is not protein-coding or does not contain start and stop codons."
        )
        return None, None, None
    assert (
        len(t.protein_sequence) * 3
        == len(t.coding_sequence) - 3
        == sum([i[1] - i[0] + 1 for i in t.coding_sequence_position_ranges])
    ), "CDS length is not compatible with the protein sequence"

    from roux.lib.set import flatten

    return (
        pd.DataFrame(
            {
                "pos": flatten(
                    [
                        list(range(i[0], i[1] + 1))
                        for i in t.coding_sequence_position_ranges
                    ]
                ),
            }
        )
        .sort_values("pos", ascending=t.strand == "+")
        .assign(
            **{
                "chrom": t.contig,
                "strand": t.strand,
                "base": list(t.coding_sequence[:-3]),
                "aa": flatten([[s] * 3 for s in list(t.protein_sequence)]),
                "aa pos": flatten(
                    [[i + 1] * 3 for i, s in enumerate(t.protein_sequence)]
                ),
            }
        )
        .astype(
            {
                "chrom": str,
                "strand": str,
                "pos": int,
                "base": str,
                "aa": str,
                "aa pos": int,
            }
        )
    )


def get_protein_mutation_coords(
    data: pd.DataFrame,
    aapos: int,
    test=False,
) -> tuple:
    """Get protein mutation coordinates

    Args:
        data (pd.DataFrame): input table
        aapos (int): amino acid position
        test (bool, optional): test-mode. Defaults to False.

    Raises:
        ValueError: invalid positions

    Returns:
        tuple: aapos,start,end,seq
    """
    # if isinstance(aapos,str):
    #     aapos=int(aapos[1:-1])
    ## identify strand
    if data.iloc[0, :]["pos"] < data.iloc[-1, :]["pos"]:
        strand = "+"
    elif data.iloc[0, :]["pos"] > data.iloc[-1, :]["pos"]:
        strand = "-"
    else:
        raise ValueError(
            f"data.iloc[0,:]['pos']({data.iloc[0,:]['pos']})=data.iloc[-1,:]['pos']({data.iloc[-1,:]['pos']})?"
        )
    ##
    if (aapos * 3) - 1 > len(data):
        logging.error(
            f"aa pos {aapos} ({(aapos*3)-1}/3) exeeds the length of the protein ({len(data)}/3)."
        )
    start, end = (
        data.iloc[(aapos - 1) * 3, :]["pos"],
        data.iloc[(aapos * 3) - 1, :]["pos"],
    )
    # seq=data.iloc[start:end,:]['base'].sum()
    seq = data.set_index("pos").loc[start:end, "base"].sum()
    if test:
        logging.debug(f"{aapos},{start},{end}")
    # data.iloc[range((aapos-1)*3,(aapos)*3),:]
    # pos=df2['pos'].tolist()[1]
    # print(int(mutation[1:-1]),pos)
    if strand == "-":
        start, end = end, start
    return aapos, start, end, seq


def map_coords(
    df_: pd.DataFrame,  # input amino acid positions
    df1_: pd.DataFrame,  # all aa positions
    verbose: bool = False,
) -> pd.DataFrame:
    """Map coordinates

    Args:
        df_ (pd.DataFrame): input table

    Returns:
        pd.DataFrame: output table
    """
    df2_ = (
        df_["aa pos"]
        .apply(
            lambda x: get_protein_mutation_coords(
                df1_,
                aapos=x,
                # search_window=search_window,
                test=verbose,
            )
        )
        .apply(pd.Series)
        .rename(
            columns={
                0: "aa pos",
                1: "start",
                2: "end",
                3: "sequence target codon",
                4: "chrom",
                5: "start flanking",
                6: "end flanking",
            }
        )
    )
    return (
        df_
        ## get coords
        .merge(
            right=df2_,
            on="aa pos",
            how="left",
            # validate="1:1", # turned off because there could be mutations on left side and aa to nt position mapping on right
        )
    )


def get_mutation_coords_protein(
    df0: pd.DataFrame,
    annots,
    search_window: int,
    outd: str = None,
    force: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    """Get mutation coordinates for protein

    Args:
        df0 (pd.DataFrame): input table
        annots (_type_): pyensembl annotations
        search_window (int): search window length on either side of the target
        outd (str, optional): output directory path. Defaults to None.
        force (bool, optional): force. Defaults to False.
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        pd.DataFrame: output table
    """
    dfs = {}
    for protein_id, df_ in df0.groupby("protein id"):
        outp = f"{outd}/{protein_id}.tsv"
        if not Path(outp).exists() or force:
            df1_ = get_protein_cds_coords(
                annots,
                protein_id,
            )
            chrom, strand = df1_["chrom"].tolist()[0], df1_["strand"].tolist()[0]
            if df1_ is None:
                continue
            if outd is not None:
                logging.info("saved protein positions: " + to_table(df1_, outp))
        else:
            df1_ = read_table(outp)
        dfs[protein_id] = df_.reset_index(drop=True).assign(  ## for join
            chrom=chrom,
            strand=strand,
        )
        # if 'mutation' in dfs[protein_id]:
        #     df2_=(dfs[protein_id]['aa pos'].apply(lambda x: get_protein_mutation_coords(
        #                 df1_,
        #                 aapos=x,
        #                 # search_window=search_window,
        #                 test=verbose,
        #             )).apply(pd.Series).rename(columns={
        #                 0:'aa pos',1:'start', 2:'end',
        #                 3:'sequence target check',
        #                 4:'chrom',5:'start flanking',6: 'end flanking',
        #         }))
        #     dfs[protein_id]=(
        #         dfs[protein_id]
        #         ## get coords
        #         .merge(df2_)
        #     )
        #     del df2_
        if "aa pos" in dfs[protein_id]:
            dfs[protein_id] = map_coords(
                dfs[protein_id],  # input amino acid positions
                df1_,  # all aa positions
                verbose=verbose,
            )
        elif "aa start" in dfs[protein_id] and "aa end" in dfs[protein_id]:
            dfs[protein_id] = pd.concat(
                [
                    ## aa start
                    (
                        map_coords(
                            dfs[protein_id].rename(
                                columns={"aa start": "aa pos"}, errors="raise"
                            ),  # input amino acid positions
                            df1_,  # all aa positions
                            verbose=verbose,
                        )
                        # .drop(['aa pos'],axis=1)
                        # .add_prefix('aa start ')
                        .rename(
                            columns={
                                "aa pos": "aa start",
                                "start": "aa start start",
                                "end": "aa start end",
                            },
                            errors="raise",
                        )
                    ),
                    ## aa end
                    (
                        map_coords(
                            dfs[protein_id]
                            .loc[:, ["aa end"]]
                            .rename(
                                columns={"aa end": "aa pos"}, errors="raise"
                            ),  # input amino acid positions
                            df1_,  # all aa positions
                            verbose=verbose,
                        )
                        .drop(["aa pos"], axis=1)
                        .add_prefix("aa end ")
                        # .rename(columns={'aa pos':'aa end'},errors='raise'),
                    ),
                ],
                axis=1,
                # names=['start','end'],
            )
            dfs[protein_id] = (
                dfs[protein_id]
                .assign(
                    **{
                        "codons": lambda df: df.apply(
                            lambda x: f"from {x['sequence target codon']} to {x['aa end sequence target codon']}",
                            axis=1,
                        ),
                        "start": lambda df: df.apply(
                            lambda x: x[
                                [
                                    "aa start start",
                                    "aa start end",
                                    "aa end start",
                                    "aa end end",
                                ]
                            ].min(),
                            axis=1,
                        ),
                        "end": lambda df: df.apply(
                            lambda x: x[
                                [
                                    "aa start start",
                                    "aa start end",
                                    "aa end start",
                                    "aa end end",
                                ]
                            ].max(),
                            axis=1,
                        ),
                    }
                )
                .drop(["sequence target codon", "aa end sequence target codon"], axis=1)
            )
    df1 = pd.concat(dfs.values(), axis=0).log.dropna()
    if df0["protein id"].nunique() != df1["protein id"].nunique():
        logging.warning(
            f"protein ids filtered out: {set(df0['protein id'].unique())-set(df1['protein id'].unique())}"
        )
    return df1


def get_mutation_coords(
    df0: pd.DataFrame,
    annots,
    search_window: int,
    verbose: bool = False,
    **kws_protein,
) -> pd.DataFrame:
    """Get mutation coordinates

    Args:
        df0 (pd.DataFrame): input table
        annots (_type_): pyensembl annotation
        search_window (int): search window length on either side of the target
        verbose (bool, optional): verbose. Defaults to False.

    Returns:
        pd.DataFrame: output table
    """
    if "pos" in df0:
        df1 = df0.assign(
            start=lambda df: df["pos"],  # todo check 1 based or 0 based
            end=lambda df: df["start"],  # todo check 1 based or 0 based
            strand="+",
        )
    elif all([c in df0 for c in ["chrom", "start", "end"]]) and "strand" not in df0:
        df1 = df0.assign(strand="+")
    elif "protein id" in df0 and "start" not in df0:
        df1 = get_mutation_coords_protein(
            df0,
            annots,
            search_window,
            verbose=False,
            **kws_protein,
        )

    df1 = df1.astype({"chrom": str, "start": int, "end": int})
    ## validate
    assert "chrom" in df1
    assert "strand" in df1
    assert all(df1["start"] <= df1["end"])
    from beditor.lib.utils import to_locus

    return df1.assign(
        **{
            "target location": lambda df: df.apply(
                lambda x: to_locus("chrom", "start", "end", "strand", x=x), axis=1
            ),
        }
    )
