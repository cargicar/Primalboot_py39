import datetime
import traceback

import utils.logging as log
import utils.timing as timing
import lin_prog.lp_problem as lp_problem
import config

reload(lp_problem)

def check_point(sigma, epsilon, bresults, lptab, specfunc, hotstart=None,
        lp=None, tabfilename=None):
    """Solve the LP for a pair (sigma, epsilon) which are passed to the spectrum
    generating function in problem_spec."""
    spectrum = specfunc(sigma, epsilon, lptab.spacetimedim)
    #log.debug('spectrum: ' + str(spectrum))
    hotstarted=False
    if lp==None:
        lp = lp_problem.LP_Problem(spectrum, lptab, parallel =
                config.point_parallel, pool_size = config.point_poolsize, 
                useLU = False)
    if hotstart is not None:
        log.info('hotstarting')
        lp.hotstart(spectrum, hotstart)
        hotstarted=True
    else:
        log.info('resetting lp')
        lp.spectrum=spectrum
        lp.reset()
    timing.t0()
    costs=[]
    status=None
    while lp.is_still_running() and lp.iter < config.max_steps:
        # check if lp has stalled
        if len(costs) > 4 and costs[-1]==costs[-2] and costs[-2]==costs[-3]:
            status=lp_problem.STATUS_STALLED
            break
        lp.fewsteps(500, status_after=config.show_status_after)
        log.debug("status after "+str(lp.iter)+" iterations ("
                +str((sigma,epsilon)) + "):"+lp.STATUS+
                " (cost: "+str(lp.curCost()) + ", vars: " + str(lp.toElim) + ")")
        costs+=[lp.curCost()]
        #print lp.MRC
    tm=timing.deltat()
    if status is None:
        status = lp.STATUS
    if lp.iter >= config.max_steps:
        log.warning('lp (%f,%f) took more max iterations (%d).'%
                    (sigma,epsilon,config.max_steps))
    log.info("time taken: "+str(tm))
    hotstart_data = hotstart
    if status == lp_problem.STATUS_COST_MINIMIZED:
        hotstart_data=lp_problem.Hotstart_Data(lp)
    lpdata =(lp_problem.to_double(lp.curCost()), lp.toElim, 
            lp.Xb, lp_problem.to_double(lp.curSol()), lp.iter,hotstarted)
    #bresults[(sigma,epsilon)]=(lp.STATUS,tm,lpdata,lptab.CBlen)
    bresults[(sigma,epsilon)]=(status,tm,lpdata,lptab.CBlen)
    return  status,hotstart_data,lp


