"""
Module to augment pandas functionality. 
"""
import pandas as pd
import numpy as np
from cStringIO import StringIO

def latex_strip(tab,path):
    """
    Strip off tab[4:-3] and write to path
    """
    
    tab = tab.getvalue()
    tab = tab.split('\n')[4:-3]
    tab = [s+'\n' for s in tab]
    with open(path,'w') as f:
        f.writelines(tab)

def LittleEndian(r):
    names = r.dtype.names
    data = {}
    for n in names:
        if r[n].dtype.byteorder=='>':
            data[n] = r[n].byteswap().newbyteorder() 
        else:
            data[n] = r[n] 
    q = pd.DataFrame(data,columns=names)
    return np.array(q.to_records(index=False))

def df_to_ndarray(df):
    """
    Convert Pandas DataFrame to ndarray
    
    If there are objects in the array, convert them to a string type.

    Parameters
    ----------
    df  : DataFrame
    res : numpy ndarray
    """

    arrayList = []
    for c in df.columns:
        if df[c].dtype==np.dtype('O'):
            arr = np.array(df[c]).astype(str)
        else:
            arr = np.array(df[c])
        arrayList += [arr]

    res = np.rec.fromarrays(arrayList,names=list(df.columns))
    res = np.array(res)
    return res  
    



