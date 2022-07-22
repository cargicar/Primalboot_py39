# same modules as passed to ppservers (used when calling functions directly)
import datetime
import traceback
import numpy as np

import tables.lp_table as lp_table
import lin_prog.lp_problem as lp_problem
import tables.tables as cb
import utils.logging as log
import scan.tabholder as tabh
import config
import check_point 
import scan.jug.jug_tasks as jt


reload(cb)
reload(lp_problem)
reload(tabh)
reload(check_point)
reload(jt)



def scan_eps(sigma,epsdata,specfunc,tabfilename, use_jug=False):
    log.proc_id=str(sigma)
    t0 = datetime.datetime.now()
    cbtab = tabh.get_tabfile(tabfilename)
    #cbtab = mt.get_table(float(dim-2)/2, nmax, mmax)
    sigmatab = cb.Sigma_Table(ds=sigma, cbtab = cbtab)
    lptab = lp_table.LP_Table(sigmatab)
    log.proc_id=str(sigma)+' d='+str(lptab.spacetimedim)
    t0 = datetime.datetime.now()
    e0=epsdata[0]
    e1=epsdata[1]
    deltaeps=epsdata[2]
    if e0 > e1:
        log.warning("min eps bound greater than max, resetting: "+
                    str((e0,e1)))
        e1=e0+2*deltaeps
    log.info('scanning: ('+str((sigma,e0,e1,deltaeps))+')')
    bresults={}
    hotstart=None
    lp=None
    try:
        for eps in np.arange(e0, e1, deltaeps):
            log.debug('checking: (%f,%f))'%(sigma,eps))
            if use_jug:
                print 'adding point (%d,%d)'%(sigma,eps)
                print sigma
                print eps
                print bresults
                print tabfilename
                status,hotstart,lp=jt.run_check_point2(sigma,eps,bresults,specfunc,tabfilename)
                #jt.run_check_point2(sigma,eps,bresults,specfunc,tabfilename)
            else:
                status,hotstart,lp=check_point.check_point(sigma, eps,bresults, lptab,specfunc,hotstart,lp)
            hotstart=None
            lp=None
            if status == lp_problem.STATUS_AUX_ELIMINATED:
                e0=eps
                log.info("found sol:"+str((sigma, eps)))
            elif status == lp_problem.STATUS_COST_MINIMIZED:
                e1=eps
                log.info("no sol:"+str((sigma, eps)))
            elif status == lp_problem.STATUS_STALLED:
                e1=eps
                log.warning("stalled on LP ("+str((sigma,eps))+
                        "), treating as no solution (cost minimzed)")
                log.info("stalled:"+str((sigma, eps)))
            else:
                raise UserWarning('unexpected status returned from lp for ('+
                                  str((sigma,eps))+'): '+status)
            #log.info("cur bounds:"+str((sigma, e0, e1, e1-e0)))
    except Exception as inst:
        log.warning('Exception checking bounds: '+str((traceback.format_exc())))
        raise
    log.debug("results: "+str(bresults)[0:80])
    t1 = datetime.datetime.now()
    tottm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    log.debug('total time: '+str(tottm))
    return bresults

