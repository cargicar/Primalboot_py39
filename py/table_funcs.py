import os
import datetime

import setpath

import solver.v1.tables.lp_table as lp_table1
import solver.v1.tables.make_table as mt
import solver.v1.tables.tables as cb
import utils.scan.tabholder as tabh
import utils.logging as log
import solver.v2.tables.lp_table as lp_table2

from version import *

import probspec as ps

import config # Carlos edit

if version==1:
    tabfilename=mt.get_filename(float(ps.dim-2)/2, config.nmax, config.mmax)
else:
    tabfilename='../tables/eps'+str(float(ps.dim-2)/2)+'n'\
            +str(config.nmax)+'m'+str(config.mmax)+'.txt'

def checktable():
    if version==1:
        #tabfilename=mt.get_filename(float(ps.dim-2)/2, config.nmax, config.mmax)
        if not os.path.isfile(tabfilename):
            mt.get_table(float(ps.dim-2)/2, config.nmax,config.mmax,make_only=False)

def loadtable():
    if version==1:
        #tabfilename=mt.get_filename(float(dim-2)/2, nmax, mmax)
        if not os.path.isfile(tabfilename):
            mt.get_table(float(dim-2)/2, nmax,mmax,make_only=False)
        tab = cb.CB_Table(FILE = tabfilename)

def gettable(sig):
    if version==1:
        #tabfilename=mt.get_filename(float(dim-2)/2, nmax, mmax)
        cbtab = tabh.get_tabfile(tabfilename)
        sigmatab = cb.Sigma_Table(ds=sig, cbtab = cbtab)
        cblen = sigmatab.CBlen
        lptab =lp_table1.LP_Table(sigmatab)
        return cblen,lptab
    else:
        #lptab = lp_table2.LP_Table("../tables/ds-0.5182-n3-m1.txt") 
        t0 = datetime.datetime.now()
        lptab = lp_table2.LP_Table(tabfilename, pfwrap(sig))
        t1 = datetime.datetime.now()
        tm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
        print ('time loading/convolving table: %s ' %tm )
        return lptab.CBlen,lptab
        
