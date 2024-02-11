import logging
import pandas as pd
import numpy as np


from beditor.lib.utils import (
    s2re,
    reverse_complement_multintseq,
    reverse_complement_multintseqreg,
    multint2reg,
    multint2regcomplement,
    get_nt2complement,
)


def dpam2dpam_strands(
    dpam: pd.DataFrame,
    pams: list,
) -> pd.DataFrame:
    """Duplicates dpam dataframe to be compatible for searching PAMs on - strand

    Args:
        dpam (pd.DataFrame): dataframe with pam information
        pams (list): pams to be used for actual designing of guides.

    Returns:
        pd.DataFrame: table
    """
    nt2complement = get_nt2complement()

    dpam = dpam.rd.clean()
    dpam["rPAM"] = dpam.apply(lambda x: s2re(x["PAM"], multint2reg), axis=1)
    dpam = dpam.set_index("PAM")
    if "strand" in dpam:
        all(dpam["strand"] == "+")
    else:
        dpam["strand"] = "+"
    dpamr = pd.DataFrame(columns=dpam.columns)
    dpam.loc[:, "reverse complement"] = np.nan
    dpam.loc[:, "original"] = np.nan
    for pam in dpam.index:
        pamr = reverse_complement_multintseq(pam, nt2complement)
        dpam.loc[pam, "reverse complement"] = pamr
        dpam.loc[pam, "original"] = pam
        dpamr.loc[pamr, "original"] = pam
        dpam.loc[pam, "original position"] = dpam.loc[pam, "PAM position"]
        dpamr.loc[pamr, "original position"] = dpam.loc[pam, "PAM position"]
        dpamr.loc[pamr, ["PAM position", "guide length"]] = dpam.loc[
            pam, ["PAM position", "guide length"]
        ]
        dpamr.loc[pamr, ["rPAM"]] = reverse_complement_multintseqreg(
            pam, multint2regcomplement, nt2complement
        )
    dpamr["PAM position"] = dpamr.apply(
        lambda x: "up" if x["PAM position"] == "down" else "down", axis=1
    )
    dpamr["strand"] = "-"
    # dpam_strands=dpam.append(dpamr,sort=True)
    dpam_strands = pd.concat([dpam, dpamr], ignore_index=True)
    dpam_strands.index.name = "PAM"
    dpam_strands.loc[:, "is a reverse complement"] = pd.isnull(
        dpam_strands.loc[:, "reverse complement"]
    )
    # print(dpam_strands)
    # print(pams)
    # pams_strands=pams+dpam_strands.loc[pams,'reverse complement'].dropna().tolist()
    # dpam_strands=dpam_strands.loc[pams_strands,:]
    return dpam_strands

def get_be2dpam(
    din: pd.DataFrame,
    methods: list = None,
    test: bool = False,
    cols_dpam: list = ["PAM", "PAM position", "guide length"],
) -> dict:
    """
    Make BE to dpam mapping i.e. dict

    Args:
        din (pd.DataFrame): table with BE and PAM info all cols_dpam needed
        methods (list, optional): method names. Defaults to None.
        test (bool, optional): test-mode. Defaults to False.
        cols_dpam (list, optional): columns to be used. Defaults to ['PAM', 'PAM position', 'guide length'].

    Returns:
        dict: output dictionary.
    """
    be2dpam = {}
    be2pam = (
        din.loc[:, ["method", "PAM"]]
        .drop_duplicates()
        .set_index("method")
        .to_dict()["PAM"]
    )
    if methods is None:
        methods = be2pam.keys()
    if test:
        logging.info(methods)
    for be in be2pam:
        if be in methods:
            pam = be2pam[be]
            dpam = din.loc[((din["PAM"] == pam) & (din["method"] == be)), cols_dpam]
            dpam_strands = dpam2dpam_strands(dpam, pams=[pam])
            if "PAM" in dpam_strands:
                dpam_strands = dpam_strands.set_index("PAM")
            be2dpam[be] = dpam_strands
    return be2dpam
