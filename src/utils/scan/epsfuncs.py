import pp
import numpy as np
import scipy 
import sqlite3
import matplotlib.pyplot as plt


import lin_prog.lp_problem as lp_problem
import tables.tables as cb
import utils.logging as log
import bisect
import config
#import rundata
reload(bisect)
reload(lp_problem)

global splrep,splev
splrep=scipy.interpolate.splrep
splev=scipy.interpolate.splev

from scipy.interpolate import interp1d

def epslambda(const_epsdata):
    """return a constant function giving eps. useful for defining trivial
    epsdata"""
    e0=lambda x: const_epsdata[0]
    e1=lambda y: const_epsdata[1]
    emin=lambda z: const_epsdata[2]
    return (e0,e1,emin)

def min_max_eps(results, absmin, absmax):
    # drop failed points
    resclean=[r for r in results if r is not None and len(r.keys()) > 0]
    # sort by sigma
    sres=sorted(resclean, key=lambda r: r.keys()[0])
    #sres=results
    sig=[r.keys()[0][0] for r in sres]
    min_eps=[max([absmin] + [x[1] for x in r if r[x][0]==lp_problem.STATUS_AUX_ELIMINATED])
            for r in sres]
    max_eps=[min([absmax]+[x[1] for x in r if r[x][0]==lp_problem.STATUS_COST_MINIMIZED])
            for r in sres]
    return sig,min_eps,max_eps

def new_epsdata(sig, eps_min,eps_max, eps_res,n):
    """return new eps data and sigma by adding 'n' new points between each two
    sigma points.
    sig- full sorted list of sigma's checked so far
    epsmin/max - full set of eps min/max data
    n - number of new points at each iteration
    returns: 
    new delta-sigma, new sigma points to check, eps min/max interpolation functions"""
    s0=np.array(sig)
    #offset=np.array(eps_max)-np.array(eps_min)
    offset=eps_res(sig)
    e0=np.array(eps_min)-4*offset
    e1=np.array(eps_max)+4*offset
    kval=min(1, len(eps_min))
    log.debug("interpolate (s): "+str(s0))
    log.debug("interpolate (e0): "+str(e0))
    log.debug("interpolate (e1): "+str(e1))
    tckmin = splrep(s0,e0,k=kval)
    tckmax = splrep(s0,e1,k=kval)
    # in some version splev returns array so we convert it to an array anyway
    # for consistency
    emin=lambda s: np.array([splev(s, tckmin,der=0)])[0]
    emax=lambda s: np.array([splev(s, tckmax,der=0)])[0]
    #emin=interp1d(s0, e0)
    #emax=interp1d(s0, e1)
    signew=[]
    ds=float(sig[1]-sig[0])/n
    for i in range(len(sig)-1):
        signew+=[sig[i] + j*ds for j in range(1,n)]
    #print emin(signew)
    #print emax(signew)
    #print emax(signew)-emin(signew)
    newvals=''
    for x in signew:
        newvals+= (str((x,emin(x),emax(x))) +',')
    log.debug('new bounds: '+newvals)
    return ds,signew,emin,emax

def new_epsdata2(sig, eps_min,eps_max, n):
    """return new eps data and sigma by adding 'n' new points between each two
    sigma points.
    sig- full sorted list of sigma's checked so far
    epsmin/max - full set of eps min/max data
    n - number of new points at each iteration
    returns: 
    new delta-sigma, new sigma points to check, eps min/max interpolation functions"""
    s0=np.array(sig)
    emin=[]
    emax=[]
    signew=[]
    ds=float(sig[1]-sig[0])/n
    for i in range(len(sig)-1):
        mn=min(sig[i],sig[i+1])
        mx=max(sig[i],sig[i+1])
        signew+=[sig[i] + j*ds for j in range(1,n)]
        emin+=[mn for j in range(1,n)]
        emax+=[mx for j in range(1,n)]
    #print emin(signew)
    #print emax(signew)
    #print emax(signew)-emin(signew)
    newvals=''
    for i in range(len(signew)):
        newvals+= (str((signew[i],emin[i],emax[i])) +',')
    log.debug('new bounds (uninterpolated): '+newvals)
    return ds,signew,emin,emax


