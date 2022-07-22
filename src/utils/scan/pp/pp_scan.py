import pp
import numpy as np
import scipy 
import sqlite3
import matplotlib.pyplot as plt


import lin_prog.lp_problem as lp_problem
import tables.tables as cb
import utils.logging as log
import config

import bounds.bisect as bisect
import scan.rundata as rundata
import scan.epsfuncs as epsfuncs
reload(bisect)
reload(lp_problem)
reload(epsfuncs)

global splrep,splev
splrep=scipy.interpolate.splrep
splev=scipy.interpolate.splev

from scipy.interpolate import interp1d



class scan:
#--------------------------------------------------------------------------
    def __init__(self, cluster=config.pp_cluster_mode):
        #self.results={}
        self.results=[]
        self.joblist=[]
        self.poolSize=4
        if cluster:
            ppservers = ("*",)
            self.job_server = pp.Server(ncpus=0,ppservers=ppservers)
        else:
            self.job_server = pp.Server(ncpus=config.pp_num_cpus,secret="hey")
#--------------------------------------------------------------------------
    def sigscan(self, description, tabfilename, specfunc, sigdata, epsdata, absbounds):
        """Generatea  sigma scan via bisection:
        description - run descriptoin, stored in DB
        tabfile - file name of the table
        sigdata - (init sigma, end sigma, init res, sig bisections)
        epsdata - (min eps func, max eps func, bisec resolution)
                    epsdata are functions of sigma!
        specfunc - a function of sigma, epsilon that returns a Spectrum
        object"""
        self.rundata=rundata.rundata(description)
        log.info('called with sigdata: '+str(sigdata))
        log.info('called with epsdata: '+str(epsdata))
        if config.precision:
            self.rundata.prec="mpfr"
        else:
            self.rundata.prec="machine"
        # load the tables to check cblen/dim
        tab = cb.CB_Table(FILE = tabfilename)
        self.rundata.spacetimedim=2*tab.eps+2
        (s0,s1,ds0,nbis)=sigdata
        (efunc0,efunc1,eres)=epsdata
        rng=np.arange(s0,s1,ds0)
        #absmin=0.8*min([efunc0(s) for s in rng])
        #absmax=1.2*max([efunc1(s) for s in rng])
        absmin=absbounds[0]
        absmax=absbounds[1]
        edata=epsdata
        for i in range(nbis):
            print 'Checking sigmas [',len(rng),']:', rng
            log.info('Checking sigmas ['+str(len(rng))+']:'+str(rng))
            self.joblist=self.queuejobs(rng,edata,absbounds,tabfilename,specfunc)
            for j in self.joblist:
                job=j()
                log.info('job done: '+str(job)[0:80])
                self.results+=[job]
            s,e0,e1=epsfuncs.min_max_eps(self.results,absmin,absmax)
            log.debug("results (s): "+str(s0))
            log.debug("results (e0): "+str(e0))
            log.debug("results (e1): "+str(e1))
            if len(s) > 1:
                ds,rng,ef0,ef1=epsfuncs.new_epsdata(s,e0,e1,eres,2)
                edata=(ef0,ef1,epsdata[2])
        self.save_results(sigdata,epsdata,absmin,absmax)
#--------------------------------------------------------------------------
    def queuejobs(self,srange,epsdata,absbounds,tabfilename,specfunc):
        joblist=[]
        epsstr='epsdata: '
        for s in srange:
            edata = [e(s) for e in epsdata]
            epsstr+=(str((s,edata))+', ')
            joblist+=[self.job_server.submit(bisect.bisect_sigma, 
                                    #(s,edata,absbounds,tabfilename,specfunc),
                                    (s,edata,absbounds,tabfilename,specfunc),
                                    depfuncs=(bisect.check_point,
                                        bisect.check_bounds,bisect.get_tabfile,
                                        specfunc),
                                    modules=("setpath",
                                             "lin_prog.lp_problem as lp_problem",
                                             "utils.timing as timing", 
                                             "tables.tables as cb","datetime", 
                                             "tables.lp_table as lp_table",
                                             "utils.logging as log",
                                             "import config",
                                             "import traceback",
                                             "import scan.tabholder as tabh"))]
        print epsstr
        log.debug(epsstr)
        return joblist
#--------------------------------------------------------------------------
    def save_results(self,sigdata,epsdata,absmin,absmax):
        if config.dbfile == None:
            return
        conn = sqlite3.connect(config.dbfile)
        c = conn.cursor()
        #c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?, -1, 0.0)",
        #             (sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
        c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?,?,?,?, -1, 0.0)",
                    (self.rundata.desc, self.rundata.prec, self.rundata.spacetimedim, 
                        sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
        runid=c.lastrowid
        # mme is in different sort order!!
        mme=epsfuncs.min_max_eps(self.results,absmin,absmax)
        mmedict={}
        for i in range(len(mme[0])):
            mmedict[mme[0][i]]=(mme[1][i],mme[2][i])
        resclean=[r for r in self.results if r is not None and len(r.keys()) > 0]
        #for res in self.results:
        for res in resclean:
            if res is None:
                log.warning("Empty entry in results list")
                continue
            s=res.keys()[0][0]
            epsfinmin,epsfinmax=mmedict[s]
            edata = [e(s) for e in epsdata]
            c.execute("INSERT INTO bisections VALUES (datetime('now'),\
                        ?,?,?,?,?,?,?,?)",(runid, s, edata[0],edata[1],edata[2],
                        epsfinmin,epsfinmax,0.0))
            bisectid=c.lastrowid
            for k in res.keys():
                status,tm,lpdata,cblen=res[k]
                #print lpdata
                c.execute("INSERT INTO points VALUES (datetime('now'),\
                            ?,?,?,?,?,?,?,?,?,?,?,?)",(bisectid, s, k[1], status,
                                lpdata[0], lpdata[1], str(lpdata[2]), str(lpdata[3]),
                                lpdata[4], lpdata[5], cblen, str(tm)))
        # Save (commit) the changes
        conn.commit()
        # We can also close the connection if we are done with it.
        #Just be sure any changes have been committed or they will be lost.
        conn.close()
#--------------------------------------------------------------------------
    def test(self):
        self.problem_module="bounds.problem_spec"
        j=self.job_server.submit(bisect_sigma, (0.115, (0.7,1.1,1e-2)),
                          depfuncs=(check_point,),
                          modules=("lin_prog.lp_problem as lp_problem",
                                   "utils.timing as timing", 
                                   "tables.tables as cb", "datetime",
                                   "tables.lp_table as lp_table",
                                   "utils.logging as log",
                                    "import config",
                                    "import traceback",
                                   self.problem_module + " as pm"))
        print j()










#--------------------------------------------------------------------------

#def epslambda(const_epsdata):
#    """return a constant function giving eps. useful for defining trivial
#    epsdata"""
#    e0=lambda x: const_epsdata[0]
#    e1=lambda y: const_epsdata[1]
#    emin=lambda z: const_epsdata[2]
#    return (e0,e1,emin)
#
#def min_max_eps(results, absmin, absmax):
#    # drop failed points
#    resclean=[r for r in results if r is not None and len(r.keys()) > 0]
#    # sort by sigma
#    sres=sorted(resclean, key=lambda r: r.keys()[0])
#    #sres=results
#    sig=[r.keys()[0][0] for r in sres]
#    min_eps=[max([absmin] + [x[1] for x in r if r[x][0]==lp_problem.STATUS_AUX_ELIMINATED])
#            for r in sres]
#    max_eps=[min([absmax]+[x[1] for x in r if r[x][0]==lp_problem.STATUS_COST_MINIMIZED])
#            for r in sres]
#    return sig,min_eps,max_eps
#
#def new_epsdata(sig, eps_min,eps_max, eps_res,n):
#    """return new eps data and sigma by adding 'n' new points between each two
#    sigma points.
#    sig- full sorted list of sigma's checked so far
#    epsmin/max - full set of eps min/max data
#    n - number of new points at each iteration
#    returns: 
#    new delta-sigma, new sigma points to check, eps min/max interpolation functions"""
#    s0=np.array(sig)
#    #offset=np.array(eps_max)-np.array(eps_min)
#    offset=eps_res(sig)
#    e0=np.array(eps_min)-4*offset
#    e1=np.array(eps_max)+4*offset
#    kval=min(1, len(eps_min))
#    log.debug("interpolate (s): "+str(s0))
#    log.debug("interpolate (e0): "+str(e0))
#    log.debug("interpolate (e1): "+str(e1))
#    tckmin = splrep(s0,e0,k=kval)
#    tckmax = splrep(s0,e1,k=kval)
#    # in some version splev returns array so we convert it to an array anyway
#    # for consistency
#    emin=lambda s: np.array([splev(s, tckmin,der=0)])[0]
#    emax=lambda s: np.array([splev(s, tckmax,der=0)])[0]
#    #emin=interp1d(s0, e0)
#    #emax=interp1d(s0, e1)
#    signew=[]
#    ds=float(sig[1]-sig[0])/n
#    for i in range(len(sig)-1):
#        signew+=[sig[i] + j*ds for j in range(1,n)]
#    #print emin(signew)
#    #print emax(signew)
#    #print emax(signew)-emin(signew)
#    newvals=''
#    for x in signew:
#        newvals+= (str((x,emin(x),emax(x))) +',')
#    log.debug('new bounds: '+newvals)
#    return ds,signew,emin,emax
#
#def new_epsdata2(sig, eps_min,eps_max, n):
#    """return new eps data and sigma by adding 'n' new points between each two
#    sigma points.
#    sig- full sorted list of sigma's checked so far
#    epsmin/max - full set of eps min/max data
#    n - number of new points at each iteration
#    returns: 
#    new delta-sigma, new sigma points to check, eps min/max interpolation functions"""
#    s0=np.array(sig)
#    emin=[]
#    emax=[]
#    signew=[]
#    ds=float(sig[1]-sig[0])/n
#    for i in range(len(sig)-1):
#        mn=min(sig[i],sig[i+1])
#        mx=max(sig[i],sig[i+1])
#        signew+=[sig[i] + j*ds for j in range(1,n)]
#        emin+=[mn for j in range(1,n)]
#        emax+=[mx for j in range(1,n)]
#    #print emin(signew)
#    #print emax(signew)
#    #print emax(signew)-emin(signew)
#    newvals=''
#    for i in range(len(signew)):
#        newvals+= (str((signew[i],emin[i],emax[i])) +',')
#    log.debug('new bounds (uninterpolated): '+newvals)
#    return ds,signew,emin,emax
#
#    
#def save_results(scan, sigdata,epsdata,absmin,absmax):
#    if config.dbfile == None:
#        return
#    conn = sqlite3.connect(config.dbfile)
#    c = conn.cursor()
#    #c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?, -1, 0.0)",
#    #             (sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
#    c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?,?,?,?, -1, 0.0)",
#                (scan.rundata.desc, scan.rundata.prec, scan.rundata.spacetimedim, 
#                    sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
#    runid=c.lastrowid
#    # mme is in different sort order!!
#    mme=min_max_eps(scan.results,absmin,absmax)
#    mmedict={}
#    for i in range(len(mme[0])):
#        mmedict[mme[0][i]]=(mme[1][i],mme[2][i])
#    resclean=[r for r in scan.results if r is not None and len(r.keys()) > 0]
#    #for res in self.results:
#    for res in resclean:
#        if res is None:
#            log.warning("Empty entry in results list")
#            continue
#        s=res.keys()[0][0]
#        epsfinmin,epsfinmax=mmedict[s]
#        edata = [e(s) for e in epsdata]
#        c.execute("INSERT INTO bisections VALUES (datetime('now'),\
#                    ?,?,?,?,?,?,?,?)",(runid, s, edata[0],edata[1],edata[2],
#                    epsfinmin,epsfinmax,0.0))
#        bisectid=c.lastrowid
#        for k in res.keys():
#            status,tm,lpdata,cblen=res[k]
#            #print lpdata
#            c.execute("INSERT INTO points VALUES (datetime('now'),\
#                        ?,?,?,?,?,?,?,?,?,?,?,?)",(bisectid, s, k[1], status,
#                            lpdata[0], lpdata[1], str(lpdata[2]), str(lpdata[3]),
#                            lpdata[4], lpdata[5], cblen, str(tm)))
#    # Save (commit) the changes
#    conn.commit()
#    # We can also close the connection if we are done with it.
#    #Just be sure any changes have been committed or they will be lost.
#    conn.close()
#

