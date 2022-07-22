import pp
import numpy as np
import scipy 
import sqlite3
import matplotlib.pyplot as plt


import lin_prog.lp_problem as lp_problem
import tables.tables as cb
import utils.logging as log
import bounds.bisect as bisect
import config
import scan.rundata 
import scan.epsfuncs as epsfuncs
import scan.tabholder as tabh
import bounds.check_point as check_point 

reload(bisect)
reload(lp_problem)
reload(check_point)
reload(tabh)

global splrep,splev
splrep=scipy.interpolate.splrep
splev=scipy.interpolate.splev

from scipy.interpolate import interp1d

from jug import TaskGenerator, barrier, value

@TaskGenerator
def run_check_point2(sigma, epsilon, bresults, specfunc, tabfilename, hotstart=None, lp=None):
    lptab=tabh.get_lptab(sigma, tabfilename)
    status,hotstart,lp=check_point.check_point(sigma, epsilon, bresults, lptab, specfunc, hotstart, lp,
            tabfilename=tabfilename)
    if status == lp_problem.STATUS_AUX_ELIMINATED:
        e0=epsilon
        log.info("found sol:"+str((sigma, epsilon)))
    elif status == lp_problem.STATUS_COST_MINIMIZED:
        e1=epsilon
        log.info("no sol:"+str((sigma, epsilon)))
    elif status == lp_problem.STATUS_STALLED:
        e1=epsilon
        log.warning("stalled on LP ("+str((sigma,epsilon))+
                "), treating as no solution (cost minimzed)")
        log.info("stalled:"+str((sigma, epsilon)))
    else:
        log.critical('unexpected status returned from lp for ('+
                          str((sigma,epsilon))+'): '+status)
    return status,bresults

@TaskGenerator
def run_bisect(s, edata, absbounds, tabfilename, specfunc):
    job=bisect.bisect_sigma(s,edata,absbounds,tabfilename,specfunc)
    return job

@TaskGenerator
def join(partials):
    return list(chain(*partials))

@TaskGenerator
def create_run(description, tabfilename, sigdata, init_eps):
    rundata=scan.rundata.rundata(description)
    log.info('called with sigdata: '+str(sigdata))
    log.info('called with initeps: '+str(init_eps))
    if config.precision:
        rundata.prec="mpfr"
    else:
        rundata.prec="machine"
    # load the tables to check cblen/dim
    tab = cb.CB_Table(FILE = tabfilename)
    rundata.spacetimedim=2*tab.eps+2
    if config.dbfile == None:
        return
    conn = sqlite3.connect(config.dbfile)
    c = conn.cursor()
    c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?,?,?,?, -1, 0.0)",
                (rundata.desc, rundata.prec, rundata.spacetimedim, 
                    sigdata[0],sigdata[1],sigdata[2], -1))
    runid=c.lastrowid
    # Save (commit) the changes
    conn.commit()
    conn.close()
    return runid

#@TaskGenerator
#def sigscan(description, tabfilename, specfunc, sigdata, init_eps, absbounds):
#    log.info('called with sigdata: '+str(sigdata))
#    log.info('called with initeps: '+str(init_eps))
#    runid=create_run(description, tabfilename)
#    for s in sigdata:
#        res=run_bisect(s,edata,absbounds,tabfilename,specfunc)
#    return res

@TaskGenerator
def run_save_bisection(res,runid, edata,absmin,absmax):
    save_bisection(res,runid,edata,absmin,absmax)


def save_bisection(res,runid,edata,absmin,absmax):
    if config.dbfile == None:
        return
    conn = sqlite3.connect(config.dbfile)
    c = conn.cursor()
    # mme is in different sort order!!
    mme=epsfuncs.min_max_eps([res],absmin,absmax)
    #mmedict=mmedict[mme[0][0]]=(mme[1][0],mme[2][0])
    #resclean=[r for r in results if r is not None and len(r.keys()) > 0]
    #for res in resclean:
    #    if res is None:
    #        log.warning("Empty entry in results list")
    #        continue
    s=res.keys()[0][0]
    epsfinmin,epsfinmax=mme[1][0],mme[2][0]
    #edata = [e(s) for e in epsdata]
    c.execute("INSERT INTO bisections VALUES (datetime('now'),\
                ?,?,?,?,?,?,?,?)",(runid, s, edata[0],edata[1],edata[2],
                epsfinmin,epsfinmax,0.0))
    bisectid=c.lastrowid
    for k in res.keys():
        status,tm,lpdata,cblen=res[k]
        c.execute("INSERT INTO points VALUES (datetime('now'),\
                    ?,?,?,?,?,?,?,?,?,?,?,?)",(bisectid, s, k[1], status,
                        lpdata[0], lpdata[1], str(lpdata[2]), str(lpdata[3]),
                        lpdata[4], lpdata[5], cblen, str(tm)))
    # Save (commit) the changes
    conn.commit()
    conn.close()


######################################################################################################


@TaskGenerator
def sigscan_old(description, tabfilename, specfunc, sigdata, epslist, absbounds):
    """Generatea  sigma scan via bisection:
    description - run descriptoin, stored in DB
    tabfile - file name of the table
    sigdata - (init sigma, end sigma, init res, sig bisections)
    epsdata - (min eps func, max eps func, bisec resolution)
                epsdata are functions of sigma!
    specfunc - a function of sigma, epsilon that returns a Spectrum
    object"""
    print "Running sigscan."
    epsdata=epsfuncs.epslambda(epslist)
    rundata=scan.rundata.rundata(description)
    log.info('called with sigdata: '+str(sigdata))
    log.info('called with epsdata: '+str(epsdata))
    if config.precision:
        rundata.prec="mpfr"
    else:
        rundata.prec="machine"
    # load the tables to check cblen/dim
    tab = cb.CB_Table(FILE = tabfilename)
    rundata.spacetimedim=2*tab.eps+2
    (s0,s1,ds0,nbis)=sigdata
    (efunc0,efunc1,eres)=epsdata
    rng=np.arange(s0,s1,ds0)
    #absmin=0.8*min([efunc0(s) for s in rng])
    #absmax=1.2*max([efunc1(s) for s in rng])
    absmin=absbounds[0]
    absmax=absbounds[1]
    edata=epsdata
    results=[]
    joblist=[]
    for i in range(nbis):
        print 'Checking sigmas [',len(rng),']:', rng
        log.info('Checking sigmas ['+str(len(rng))+']:'+str(rng))
        #joblist=queuejobs(rng,edata,absbounds,tabfilename,specfunc)
        for s in sigdata:
            edata = [e(s) for e in epsdata]
            #epsstr+=(str((s,edata))+', ')
            joblist+=[run_bisect(s,edata,absbounds,tabfilename,specfunc)]
        results=join(joblist)
        #barrier()
        #results=value(joblist)
        #results=joblist
        print 'got results: ',results
        #for j in self.joblist:
        #    res=j.result()
        #    log.info('job done: '+str(res)[0:80])
        #    self.results+=[res]
        s,e0,e1=epsfuncs.min_max_eps(results,absmin,absmax)
        log.debug("results (s): "+str(s0))
        log.debug("results (e0): "+str(e0))
        log.debug("results (e1): "+str(e1))
        if len(s) > 1:
            ds,rng,ef0,ef1=epsfuncs.new_epsdata(s,e0,e1,eres,2)
            edata=(ef0,ef1,epsdata[2])
    save_results(results, rundata, sigdata,epsdata,absmin,absmax)
