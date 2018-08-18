#!usr/bin/python

# Copyright 2018, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_strs``
================================
"""
import logging

def s2re(s,ss2re):
    for ss in ss2re:
        s=s.replace(ss,ss2re[ss])
    return s

import logging
import os.path

def get_datetime():
    import datetime
    return make_pathable_string(str(datetime.datetime.now())).replace('-','_')

def get_logger(program='program',argv=None,level=None,dp=None):
# def initialize_logger(output_dir):
    cmd='_'.join([str(s) for s in argv]).replace('/','_')
    if dp is None:
        dp=''
    else:
        dp=dp+'/'
    date=get_datetime()
    logp=f"{dp}.log_{program}_{date}_{cmd}.log"
    log_format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s'
    #'[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..):%(lineno)d: %(message)s'
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # create console handler and set level to info
    handler = logging.StreamHandler()
    handler.setLevel(level)
    formatter = logging.Formatter(log_format)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

#     # create error file handler and set level to error
#     handler = logging.FileHandler(os.path.join(output_dir, "error.log"),"w", encoding=None, delay="true")
#     handler.setLevel(logging.ERROR)
#     formatter = logging.Formatter(log_format)
#     handler.setFormatter(formatter)
#     logger.addHandler(handler)

    # create debug file handler and set level to debug
    handler = logging.FileHandler(logp)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(log_format)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logp

def isstrallowed(s,form):
    """
    Checks is input string conforms to input regex (`form`).

    :param s: input string.
    :param form: eg. for hdf5: `"^[a-zA-Z_][a-zA-Z0-9_]*$"`
    """
    import re
    match = re.match(form,s)
    return match is not None

def convertstr2format(col,form):
    """
    Convert input string to input regex (`form`).
    
    :param col: input string.
    :param form: eg. for hdf5: `"^[a-zA-Z_][a-zA-Z0-9_]*$"`
    """
    if not isstrallowed(col,form):
        col=col.replace(" ","_") 
        if not isstrallowed(col,form):
            chars_disallowed=[char for char in col if not isstrallowed(char,form)]
            for char in chars_disallowed:
                col=col.replace(char,"_")
    return col

def make_pathable_string(s,replacewith='_'):
    """
    Removes symbols from a string to be compatible with directory structure.

    :param s: string
    """
    import re
    return re.sub('\W+',replacewith, s )

def linebreaker(l,break_pt=16):
    """
    used for adding labels in plots.

    :param l: list of strings
    :param break_pt: number, insert new line after this many letters 
    """

    l_out=[]
    for i in l:
        if len(i)>break_pt:
            i_words=i.split(' ')
            i_out=''
            line_len=0
            for w in i_words:
                line_len+=len(w)+1
                if i_words.index(w)==0:
                    i_out=w
                elif line_len>break_pt:
                    line_len=0
                    i_out="%s\n%s" % (i_out,w)
                else:
                    i_out="%s %s" % (i_out,w)
            l_out.append(i_out)    
#             l_out.append("%s\n%s" % (i[:break_pt],i[break_pt:]))
        else:
            l_out.append(i)
    return l_out

def splitlabel(label,splitby=' ',ctrl='__'):
    """
    used for adding labels in plots.

    :param label: string
    :param splitby: string split the label by this character/string
    :param ctrl: string, marker that denotes a control condition  
    """
    splits=label.split(splitby)
    if len(splits)==2:
        return splits
    elif len(splits)==1:

        return splits+[ctrl]

def get_time():
    """
    Gets current time in a form of a formated string. Used in logger function.

    """
    import datetime
    time=make_pathable_string('%s' % datetime.datetime.now())
    return time.replace('-','_').replace(':','_').replace('.','_')