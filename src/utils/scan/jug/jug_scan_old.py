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
reload(bisect)
reload(lp_problem)

global splrep,splev
splrep=scipy.interpolate.splrep
splev=scipy.interpolate.splev

from scipy.interpolate import interp1d

from jug import TaskGenerator, barrier, value

@TaskGenerator
def run_bisect(s, edata, absbounds, tabfilename, specfunc):
    job=bisect.bisect_sigma(s,edata,absbounds,tabfilename,specfunc)
    return job

@TaskGenerator
def join(partials):
    return list(chain(*partials))

@TaskGenerator
def sigscan(description, tabfilename, specfunc, sigdata, epslist, absbounds):
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
#-----------------------------------------------------------------------
def queuejobs(srange,epsdata,absbounds,tabfilename,specfunc):
    joblist=[]
    epsstr='epsdata: '
    for s in srange:
        edata = [e(s) for e in epsdata]
        epsstr+=(str((s,edata))+', ')
        joblist+=[run_bisect(s,edata,absbounds,tabfilename,specfunc)]
    print epsstr
    log.debug(epsstr)
    return joblist
#-----------------------------------------------------------------------
def save_results(results, rundata, sigdata,epsdata,absmin,absmax):
    if config.dbfile == None:
        return
    conn = sqlite3.connect(config.dbfile)
    c = conn.cursor()
    #c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?, -1, 0.0)",
    #             (sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
    c.execute("INSERT INTO runs VALUES (datetime('now'), ?,?,?,?,?,?,?, -1, 0.0)",
                (rundata.desc, rundata.prec, rundata.spacetimedim, 
                    sigdata[0],sigdata[1],sigdata[2], sigdata[3]))
    runid=c.lastrowid
    # mme is in different sort order!!
    mme=epsfuncs.min_max_eps(results,absmin,absmax)
    mmedict={}
    for i in range(len(mme[0])):
        mmedict[mme[0][i]]=(mme[1][i],mme[2][i])
    resclean=[r for r in results if r is not None and len(r.keys()) > 0]
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


class scanold:
#--------------------------------------------------------------------------
    def __init__(self, cluster=config.pp_cluster_mode):
        #self.results={}
        self.results=[]
        self.joblist=[]
        self.poolSize=4
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
        print "Running sigscan."
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
            barrier()
            self.results=value(self.joblist)
            print 'got results: ',self.results
            #for j in self.joblist:
            #    res=j.result()
            #    log.info('job done: '+str(res)[0:80])
            #    self.results+=[res]
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
            joblist+=[run_bisect(s,edata,absbounds,tabfilename,specfunc)]
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


