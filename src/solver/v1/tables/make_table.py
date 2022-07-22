import sys, os

#pat=os.getcwd()
pat = os.path.dirname(__file__)
sys.path.append(pat) # Carlos edit

import tables as cb
from importlib import reload # Carlos add  
reload(cb)
import numpy as np
import os.path
import config

import solver.v1.prec_float.prec_float as PF

def get_filename(eps, nmax, mmax):
    return config.tabledir+"/CBeps"+str(eps)+"n"+str(nmax)+"m"+str(mmax)+".txt"


def get_table(eps, nmax, mmax, make_only=False, lmax=40):
    """Returns a CB_Table object for a table with the given params assumed to be
    in the default table dir.  If no table exists one is created."""
    filename=get_filename(eps,nmax,mmax)
    if os.path.isfile(filename):
        if make_only:
            return
        else:
            tab = cb.CB_Table(FILE = filename)
            return tab
    else:
        print ("File %s does not exist, creating."%filename)
    lowstep= [(0.01,2),(0.05,10),(0.25,10),(1,lmax)]
    highstep= [(0.01,2),(0.05,10),(0.25,10),(1,lmax)]
    steps={}
    # 2nd arg to range should be lmax+2
    for l in range(0,lmax+2,2):
        if l < 0:
            steps[l]=lowstep
        else:
            steps[l]=highstep
    dlist = cb.make_dlist1(eps,lmax, lmax+6, steps, eps+1e-10)
    tab = cb.CB_Table(eps = eps, dlist = dlist, nmax = nmax, mmax = mmax)
    FILE = tab.save(filename)
    if make_only:
        return
    return tab
