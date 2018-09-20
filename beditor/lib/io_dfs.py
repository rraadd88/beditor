#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_dfs``
================================
"""
from os.path import basename,exists
import pandas as pd
import numpy as np
from beditor.lib.io_nums import is_numeric

import logging

def set_index(data,col_index):
    """
    Sets the index if the index is not present

    :param data: pandas table 
    :param col_index: column name which will be assigned as a index
    """
    if col_index in data:
        data=data.reset_index().set_index(col_index)
        if 'index' in data:
            del data['index']
        return data
    elif data.index.name==col_index:
        return data
    else:
        logging.error("something's wrong with the df")
        df2info(data)

# dfs
def concat_cols(df1,df2,idx_col,df1_cols,df2_cols,
                df1_suffix,df2_suffix,wc_cols=[],suffix_all=False):
    """
    Concatenates two pandas tables 

    :param df1: dataframe 1
    :param df2: dataframe 2
    :param idx_col: column name which will be used as a common index 
    """

    df1=df1.set_index(idx_col)
    df2=df2.set_index(idx_col)    
    if not len(wc_cols)==0:
        for wc in wc_cols:
            df1_cols=df1_cols+[c for c in df1.columns if wc in c]
            df2_cols=df2_cols+[c for c in df2.columns if wc in c]
    combo=pd.concat([df1.loc[:,df1_cols],df2.loc[:,df2_cols]],axis=1)
    # find common columns and rename them
    # print df1_cols
    # print df2_cols    
    if suffix_all:
        df1_cols=["%s%s" % (c,df1_suffix) for c in df1_cols]
        df2_cols=["%s%s" % (c,df2_suffix) for c in df2_cols]
        # df1_cols[df1_cols.index(col)]="%s%s" % (col,df1_suffix)
        # df2_cols[df2_cols.index(col)]="%s%s" % (col,df2_suffix)
    else:
        common_cols=[col for col in df1_cols if col in df2_cols]
        for col in common_cols:
            df1_cols[df1_cols.index(col)]="%s%s" % (col,df1_suffix)
            df2_cols[df2_cols.index(col)]="%s%s" % (col,df2_suffix)
    combo.columns=df1_cols+df2_cols
    combo.index.name=idx_col
    return combo

def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)
     

def get_colmin(data):
    """
    Get rowwise column names with minimum values

    :param data: pandas dataframe
    """
    data=data.T
    colmins=[]
    for col in data:
        colmins.append(data[col].idxmin())
    return colmins

def fhs2data_combo(fhs,cols,index,labels=None,col_sep=': '):
    """
    Collates data from multiple csv files

    :param fhs: list of paths to csv files
    :param cols: list of column names to concatenate
    :param index: name of the column name to be used as the common index of the output pandas table 
    """

    if labels is None:
        labels=[basename(fh) for fh in fhs]
    if len(fhs)>0:
        for fhi,fh in enumerate(fhs):
            label=labels[fhi]
            data=pd.read_csv(fh).set_index(index)
            if fhi==0:
                data_combo=pd.DataFrame(index=data.index)
                for col in cols:
                    data_combo.loc[:,'%s%s%s' % (label,col_sep,col)]=data.loc[:,col]
            else:
                for col in cols:
                    data_combo.loc[:,'%s%s%s' % (label,col_sep,col)]=data.loc[:,col]    
        return data_combo
    else:
        logging.error('no fhs found: len(fhs)=0')

def fhs2data_combo_appended(fhs, cols=None,labels=None,labels_coln='labels',sep=','):
    """
    Collates data from multiple csv files vertically

    :param fhs: list of paths to csv files
    :param cols: list of column names to concatenate
    """    
    if labels is None:
        labels=[basename(fh) for fh in fhs]
    if len(fhs)>0:
        data_all=pd.DataFrame(columns=cols)
        for fhi,fh in enumerate(fhs):
            label=labels[fhi]
            data=pd.read_csv(fh,sep=sep)
            if len(data)!=0:
                data.loc[:,labels_coln]=label
                if not cols is None:
                    data=data.loc[:,cols]
                data_all=data_all.append(data,sort=True)
        return data_all

def rename_cols(df,names,renames=None,prefix=None,suffix=None):
    """
    rename columns of a pandas table

    :param df: pandas dataframe
    :param names: list of new column names
    """
    if not prefix is None:
        renames=[ "%s%s" % (prefix,s) for s in names]
    if not suffix is None:    
        renames=[ "%s%s" % (s,suffix) for s in names]
    if not renames is None:
        for i,name in enumerate(names):
#             names=[renames[i] if s==names[i] else s for s in names]    
            rename=renames[i]    
            df.loc[:,rename]=df.loc[:,name]
        df=df.drop(names,axis=1)
        return df 


def reorderbydf(df2,df1):
    """
    Reorder rows of a dataframe by other dataframe

    :param df2: input dataframe
    :param df1: template dataframe 
    """
    df3=pd.DataFrame()
    for idx,row in df1.iterrows():
        df3=df3.append(df2.loc[idx,:])
    return df3

def df2unstack(df,coln='columns',idxn='index',col='value'):
    if df.columns.name is None:
        df.columns.name=coln
    if df.index.name is None:
        df.index.name=idxn
    df=df.unstack()
    df.name=col
    return pd.DataFrame(df).reset_index()

def df2info(df,col_searches=None):
    if len(df.columns)>5:
        print('**COLS**: ',df.columns.tolist())
    print('**HEAD**: ',df.loc[:,df.columns[:5]].head())
    print('**SHAPE**: ',df.shape)
    if not col_searches is None:
        cols_searched=[c2 for c1 in col_searches for c2 in df if c1 in c2]
        print('**SEARCHEDCOLS**:\n',cols_searched)
        print('**HEAD**: ',df.loc[:,cols_searched].head())
    
def lambda2cols(df,lambdaf,in_coln,to_colns):
    df_=df.apply(lambda x: lambdaf(x[in_coln]),
                 axis=1).apply(pd.Series)
    df_.columns=to_colns
    df=df.join(df_)        
    return df

def df2chucks(din,chunksize,outd,fn,return_fmt='\t',force=False):
    """
    :param return_fmt: '\t': tab-sep file, lly, '.', 'list': returns a list
    """
    from os.path import exists#,splitext,dirname,splitext,basename,realpath
    from os import makedirs
    din.index=range(0,len(din),1)

    chunkrange=list(np.arange(0,len(din),chunksize))
    chunkrange=list(zip([c+1 if ci!=0 else 0 for ci,c in enumerate(chunkrange)],chunkrange[1:]+[len(din)-1]))

    chunk2range={}
    for ri,r in enumerate(chunkrange):    
        chunk2range[ri+1]=r

    if not exists(outd):
        makedirs(outd)
    chunks=[]
    chunkps=[]
    for chunk in chunk2range:
        chunkp='{}/{}_chunk{:08d}.tsv'.format(outd,fn,chunk)
        rnge=chunk2range[chunk]
        din_=din.loc[rnge[0]:rnge[1],:]
        if not exists(chunkp) or force:
            if return_fmt=='list':
                chunks.append(din_)
            else:
                din_.to_csv(chunkp,sep=return_fmt)
            del din_
        chunkps.append(chunkp)
    if return_fmt=='list':
        return chunks
    else:
        return chunkps        

def filldiagonal(df,cols,filler=None):
    try:
        d=df.loc[cols,cols]
    except:
        logging.error('cols should be in cols and idxs')
    if filler is None:
        filler=np.nan
    for r,c in zip(cols,cols):
        df.loc[r,c]=filler
    return df

def df2submap(df,col,idx,aggfunc='sum',binary=False,binaryby='nan'):
    df['#']=1
    dsubmap=df.pivot_table(columns=col,index=idx,values='#',
                       aggfunc=aggfunc)
    if binary:
        if binaryby=='nan':
            dsubmap=~pd.isnull(dsubmap)
        else:
            dsubmap=dsubmap!=binaryby
    return dsubmap

def completesubmap(dsubmap,fmt,
    fmt2vals={'aminoacid':["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"], 
    'aminoacid_3letter':['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'],
    'codon':["TTT",    "TTC",    "TTA",  "TTG",  "TCT",  "TCC",  "TCA",  "TCG",  "TAT",  "TAC",  "TAA",  "TAG",  "TGT",  "TGC",  "TGA",  "TGG",  "CTT",  "CTC",  "CTA",  "CTG",  "CCT",  "CCC",  "CCA",  "CCG",  "CAT",  "CAC",  "CAA",  "CAG",  "CGT",  "CGC",  "CGA",  "CGG",  "ATT",  "ATC",  "ATA",  "ATG",  "ACT",  "ACC",  "ACA",  "ACG",  "AAT",  "AAC",  "AAA",  "AAG",  "AGT",  "AGC",  "AGA",  "AGG",  "GTT",  "GTC",  "GTA",  "GTG",  "GCT",  "GCC",  "GCA",  "GCG",  "GAT",  "GAC",  "GAA",  "GAG",  "GGT",  "GGC",  "GGA",  "GGG"],
    'nucleotide': ['A','T','G','C'],}):
    
    vals=fmt2vals[fmt]
    for v in vals: 
        if not v in dsubmap.columns:            
            dsubmap[v]=np.nan
        if not v in dsubmap.index:
            dsubmap.loc[v,:]=np.nan
    return dsubmap.loc[vals,vals]

def dfswapcols(df,cols):
    df[f"_{cols[0]}"]=df[cols[0]].copy()
    df[cols[0]]=df[cols[1]].copy()
    df[cols[1]]=df[f"_{cols[0]}"].copy()
    df=df.drop([f"_{cols[0]}"],axis=1)
    return df
