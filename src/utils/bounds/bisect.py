# same modules as passed to ppservers (used when calling functions directly)
import datetime
import traceback

import config

import scan.tabholder as tabh
import tables.lp_table as lp_table
import lin_prog.lp_problem as lp_problem
import tables.tables as cb
import utils.logging as log
import check_point as cp

reload(cb)
reload(lp_problem)
reload(config)
reload(cp)

def check_bounds(sig, eps0, eps1,absbounds,bresults,lptab,specfunc,hotstart,lp):
    status=''
    mineps=eps0
    maxeps=eps1
    avgeps=(eps0+eps1)/2
    difeps=(eps1-eps0)/2
    absmin=absbounds[0]
    absmax=absbounds[1]
    tries=0
    maxtries=20
    while tries < maxtries and maxeps <= absmax:
        log.info('checking upper bound [try %d of %d]: (%f, %f)' %
                (tries,maxtries,sig,maxeps))
        status,hotstart,lp=cp.check_point(sig,maxeps,bresults,lptab,specfunc,hotstart,lp)
        if status == lp_problem.STATUS_COST_MINIMIZED:
            log.info('found good upper bound: '+str((sig,maxeps)))
            break
        if status == lp_problem.STATUS_STALLED:
            log.info('LP stalled upper bound (this is ok): '+str((sig,maxeps)))
            break
        elif status == lp_problem.STATUS_AUX_ELIMINATED:
            log.info('Not an upper bound (solution exists): '+str((sig,maxeps)))
            mineps = maxeps
        else:
            raise UserWarning('unexpected status returned from lp for ('+
                    str((sig,maxeps))+'): '+status)
        tries+=1
        # if we tried absmax once and failed we break
        if maxeps == absmax:
            break
        # if we overshoot absmax just try it once: never go above it
        maxeps = min(absmax, avgeps + 4**tries * difeps)
    if status == lp_problem.STATUS_AUX_ELIMINATED:
        log.warning('could not find upper bound for '+str((sig,maxeps))+
                 ' after '+str(maxtries)+' tries. Failing!')
        raise Exception('Could not find upper eps bound for sigma '
                +str((sig,maxeps)))
    if maxeps >= absmax:
        log.warning('eps exceeding/equaled max allowed value: '+str((maxeps, absmax)))
        maxeps=absmax
        log.warning('setting eps bound to max allowed value: '
                +str((sig,maxeps)))
    # if tries > 0 we already know a lower bound so we're done
    if tries > 0:
        log.info('found good bounds (s,e0,e1): '+str((sig,mineps,maxeps)))
        return (mineps,maxeps,lp,hotstart)
    tries=0
    while tries < maxtries and mineps >= absmin:
        log.info('checking lower bound [%d of %d]: (%f, %f)' %
                (tries,maxtries,sig,mineps))
        status,hotstart,lp=cp.check_point(sig,mineps,bresults,lptab,specfunc,hotstart,lp)
        if status == lp_problem.STATUS_AUX_ELIMINATED:
            log.info('found good lower bound: '+str((sig,mineps)))
            break
        elif status == lp_problem.STATUS_COST_MINIMIZED:
            log.info('not a lower bound (no solution): '+str((sig,mineps)))
            maxeps = mineps
        elif status == lp_problem.STATUS_STALLED:
            log.info('LP stalled on lower bound (treating as if minimzed): '+str((sig,mineps)))
            maxeps = mineps
        else:
            log.warning('unexpected status returned from lp for ('+
                    str((sig,mineps))+'): '+status)
            raise UserWarning('unexpected status returned from lp for ('+
                    str((sig,mineps))+'): '+status)
        tries+=1
        # if we tried absmax once and failed we break
        if mineps == absmin:
            break
        mineps = max(absmin, avgeps - 4**tries * difeps)
    if status == lp_problem.STATUS_COST_MINIMIZED or\
        status == lp_problem.STATUS_STALLED:
        log.warning('could not find lower bound for '+str((sig,mineps))+
                 ' after '+str(maxtries)+' tries. Failing!')
        raise Exception('Could not find lower eps bound for sigma '
                +str((sig,mineps)))
    if mineps <= absmin:
        log.warning('eps below/equal to min allowed value: '+str((mineps, absmin)))
        mineps=absmin
        log.warning('setting eps bound to min allowed value:'
                     +str((sig,mineps)))
    # if tries > 1 we already know a lower bound so we're done
    log.info('found good bounds (s,e0,e1): '+str((sig,mineps,maxeps)))
    return (mineps,maxeps,lp,hotstart)


def bisect_sigma(sigma,epsdata,absbounds,tabfilename,specfunc,checkbounds=True):
#def bisect_sigma(sigma, epsdata,absbounds,tabfilename,checkbounds=True):
    log.proc_id=str(sigma)
    t0 = datetime.datetime.now()
    cbtab = tabh.get_tabfile(tabfilename)
    #sigmatab = cb.Sigma_Table(ds=sigma, cbtab = pm.tab)
    sigmatab = cb.Sigma_Table(ds=sigma, cbtab = cbtab)
    lptab = lp_table.LP_Table(sigmatab)
    log.proc_id=str(sigma)+' d='+str(lptab.spacetimedim)
    t0 = datetime.datetime.now()
    e0=epsdata[0]
    e1=epsdata[1]
    minres=epsdata[2]
    log.info("Using CB table "+tabfilename+" with "+str(lptab.CBlen)+" derivs.")
    if e0 > e1:
        log.warning("min eps bound greater than max, resetting: "+
                    str((e0,e1)))
        e1=e0+2*minres
    if e0 < absbounds[0]:
        log.warning("epsmin %f below allowed min %f. setting to min."
                    %(e0,absbounds[0]))
        e0=absbounds[0]
    if e1 > absbounds[1]:
        log.warning("epsmax %f above allowed max %f. setting to max."
                    %(e1,absbounds[1]))
        e1=absbounds[1]
    log.info('bisecting: ('+str((sigma,e0,e1,minres))+')')
    bresults={}
    hotstart=None
    lp=None
    try:
        if checkbounds:
            (e0,e1,lp,hotstart)=check_bounds(sigma,e0,e1,absbounds,bresults,lptab,specfunc,hotstart,lp) 
        eps=(e0+e1)/2
        while abs(e1-e0) > minres:
            status,hotstart,lp=cp.check_point(sigma, eps,bresults, lptab,specfunc,hotstart,lp)
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
            log.info("cur bounds:"+str((sigma, e0, e1, e1-e0)))
            eps=(e0+e1)/2
    except Exception as inst:
        log.warning('Exception checking bounds: '+str((traceback.format_exc())))
    log.info("final bounds:"+str((sigma, e0, e1, e1-e0)))
    log.debug("results: "+str(bresults)[0:80])
    t1 = datetime.datetime.now()
    tottm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    log.debug('total time: '+str(tottm))
    return bresults

