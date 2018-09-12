#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_nums``
================================
"""
import numpy as np
import pandas as pd 

def is_numeric(obj):
    """
    This detects whether an input object is numeric or not.

    :param obj: object to be tested.
    """
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        return True
    
def str2num(x):
    """
    This extracts numbers from strings. eg. 114 from M114R.

    :param x: string
    """
    return int(''.join(ele for ele in x if ele.isdigit()))

def str2numorstr(x,method=int):
    """
    This extracts numbers from strings. eg. 114 from M114R.

    :param x: string
    """
    try:
        x=method(x)
        return x
    except:
        return x

def plog(x,p = 0.5):
    """
    psudo-log

    :param x: number
    :param p: number added befor logarithm 
    """

    return np.log2(x+p)

def glog(x,l = 2):
    """
    Generalised logarithm

    :param x: number
    :param p: number added befor logarithm 

    """
    return np.log((x+np.sqrt(x**2+l**2))/2)/np.log(l)

def float2int(x):
    """
    converts floats to int when only float() is not enough.

    :param x: float
    """
    if not pd.isnull(x):
        if is_numeric(x):
            x=int(x)
    return x    

def rescale(a,mn=None):
    a=(a-a.min())/(a.max()-a.min())
    if not mn is None:
        a=1-a
        a=a*(1-mn)
        a=1-a        
    return a