import pandas as pd

def df2unstack(df,coln='columns',idxn='index',col='value'):
    if df.columns.name is None:
        df.columns.name=coln
    if df.index.name is None:
        df.index.name=idxn
    df=df.unstack()
    df.name=col
    return pd.DataFrame(df).reset_index()