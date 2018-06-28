import pandas as pd

def df2unstack(df,coln='columns',idxn='index',col='value'):
    if df.columns.name is None:
        df.columns.name=coln
    if df.index.name is None:
        df.index.name=idxn
    df=df.unstack()
    df.name=col
    return pd.DataFrame(df).reset_index()

def df2info(df):
    if len(df.columns)>5:
        print('**COLS**: ',df.columns.tolist())
    print('**HEAD**: ',df.loc[:,df.columns[:5]].head())
    print('**SHAPE**: ',df.shape)