import setpath
import numpy as np
import os
import sys
import datetime
import pickle


import utils.logging as log


import config
import table_funcs as tf

# defines lp_problem, mpfr_array, prec_float, etc...specific to version
from version import  Spectrum, lp_problem, NP_PREC, mpfr_array, prec_float, to_double, TD

from importlib import reload
reload(lp_problem) 

from probspec import dim, opedim, opespin, s0, s1, sdelta, threshold, epsmin

def pf(x):
    return prec_float.prec_float(x,prec=212)
 
def reap_fewsteps(sig, n):
    cblen,lptab=tf.gettable(sig)
    spectrum = Spectrum(spacetimedim = pf(lptab.spacetimedim),
                    deltamax = pf(lptab.lmax+1.5),
                    lmax = lptab.lmax)
    lp = lp_problem.LP_Problem(spectrum, lptab, parallel =
                    config.point_parallel, pool_size = config.point_poolsize,
                useLU = False)
    while lp.is_still_running():
            #t0()    
            #lp.fewsteps(500)
        lp.fewsteps(n)     
    hotstartdata=lp_problem.Hotstart_Data(lp)
    lp.reset()
    lp.hotstart(spectrum, hotstartdata)
    lp.opebound(opedim,opespin)
    xp_at_n_1=lp.fewsteps(n)
    xp_at_n=lp.fewsteps(4)
    return [xp_at_n_1, xp_at_n]

def setlogfile(logfile):
    config.logfile=logfile

niter=config.iters # Carlos edit

#output_file=f"xbs_{config.nmax}" # Carlos edit
def findbound(sig):
    #os.remove("xbs")
    # emin=epsmin
    emin=pf(epsmin) # Carlos edit
    log.proc_id=str(sig)
    log.debug2('entered findbound')
    t0 = datetime.datetime.now()
    #tabfilename=mt.get_filename(float(dim-2)/2, nmax, mmax)
    #cbtab = tabh.get_tabfile(tabfilename)
    #sigmatab = cb.Sigma_Table(ds=sig, cbtab = cbtab)
    #cblen = sigmatab.CBlen
    #lptab = tables.lp_table.LP_Table(sigmatab)
    log.stats('before get table')
    cblen,lptab=tf.gettable(sig)
    log.debug2('Tables loaded, cblen: %d'%cblen)
    log.stats('after get table')
    log.debug2('starting feasible vector search')
    log.stats('starting feasable vector search')
    while emin > pf(0): # Carlos edit
        spectrum = Spectrum(spacetimedim = pf(lptab.spacetimedim),
                    deltamax = pf(lptab.lmax+1.5),
                    lmax = lptab.lmax,
                    de = pf(emin))
        log.debug2('about to call LP_problem')
        lp = lp_problem.LP_Problem(spectrum, lptab, parallel =
                config.point_parallel, pool_size = config.point_poolsize,
                useLU = False)               
        # find a feasable vector
        while lp.is_still_running():
            #t0()    
            lp.fewsteps(niter)
            #lp.fewsteps(300)                         # Carlos edit
            log.debug2('done %s iterations' %(niter)  )
            #print deltat()
        log.debug("status: %s"%lp.STATUS)
        log.stats('feasable search (%d iterations)'%lp.iter)
        if lp.STATUS != lp_problem.STATUS_AUX_ELIMINATED:
            emin = pf(0.9) * emin
            log.debug("Could not find feasible solution, reducing epsmin.")
            log.debug("resetting emin to: "+str(emin))
        else:
            log.debug("Found feasable solution.")
            log.debug("leading scalar (feasible sol): "+str(lp.Xb[0][2]))
            break
    log.debug2('found feasible')
    #log.debug("solution:"+str(lp.curSol()))
    #print '\n\n\nRerunning with OPE bound.'
    oldbnd=None
    bnd=None
    t1 = datetime.datetime.now()
    tm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    log.debug('found feasable vector after %d iterations in %f sec.'%(lp.iter,
        tm))
    log.stats('after found feasible vector')
    hotstartdata=lp_problem.Hotstart_Data(lp)
    lp.reset()
    lp.hotstart(spectrum, hotstartdata)
    #print ('***XB: %s' %(lp.Xb))
    lp.opebound(opedim,opespin)
    log.debug2('starting optimization')
    #Xbs=[]                              # Carlos edit
    count_step=0
    while lp.is_still_running():
        #t0()    
        count_step+=1
        #print("COUNT_STEP %s" %count_step)
        lp.fewsteps(niter)                         # Carlos edit
        log.debug2('done %s iterations' %(niter)  )
    #    print lp.MRC
        #print deltat()
        bind=[i for i in range(len(lp.Xb)) if lp.Xb[i][0]=='BOUND']
        if len(bind) > 0:
            bnd= lp.curSol()[bind[0]]
            log.debug('CURRENT BOUND (%s) : %s ' %(sig, to_double(bnd))) # Carlos edit
            log.debug('CURRENT Cb: %s' %( to_double(lp.Cb[bind[0]]) ))
            log.debug('OLD BOUND: %s' %(to_double(oldbnd) ))
        #print to_double(lp.Cb)
        if oldbnd is not None and\
            abs(to_double(oldbnd-bnd))/to_double(bnd) < config.threshold: # Carlos edit
            log.debug('~~~~~~~~bound not varying much.  STOPPING')
            break
        if oldbnd is not None:
            log.debug('Bound change ratio: '+str(abs(to_double(oldbnd)-to_double(bnd))/to_double(bnd)))
        if bnd is not None:
            log.debug2('copying NP_PREC')
            oldbnd=NP_PREC.copy(bnd)
        
        #Xbs.append(lp.Xb) # Carlos edit
#         with open(output_file, "ab") as f:  # Carlos edit
#             #f.write("%r\n" %lp.Xb)
#             pickle.dump(lp.Xb, f)
        log.stats('bound optimization (%d iterations)'%lp.iter)
        #print
    log.debug("status:"+str(lp.STATUS))
    
    #log.debug("solution:"+str(lp.curSol()))
    #epsdim+=[lp.Xb[0][2]]
    log.debug("leading scalar (ope maxed): "+str(lp.Xb[0][2]))
    t2 = datetime.datetime.now()
    tm=(t2-t1).seconds+((t2-t1).microseconds*1e-6)
    log.debug('done maximizing bound after %d iterations in %f sec.'%(lp.iter,
        tm))
    bind=[i for i in range(len(lp.Xb)) if lp.Xb[i][0]=='BOUND']
    #specs += [ lp.Xb ]
    if len(bind) > 0:
        #boundlist+= [ (sig, lp.curSol()[bind[0]], lp.STATUS, lp.curCost()) ]
        #cTlist+= [  (sig, (sig**2 * 4**dim * dim)/((dim-1)*TD(lp.curSol()[bind[0]])) ) ]
        bound=(sig, to_double(lp.curSol()[bind[0]]), lp.STATUS,
                to_double(lp.curCost())) 
        cT=(sig, (sig**2 * 4**dim * dim)/((dim-1)*TD(lp.curSol()[bind[0]])) ) 
        #print 'BOUND: ', boundlist[-1]
        log.debug('BOUND: %s' %str(bound ))
        log.debug('Cb: %s'%to_double((lp.Cb[bind[0]])))
    else:
        bound=(sig, 0)
        cT=(sig, 0)
        log.debug('***NO BOUND FOUND!!!')
    #return bound,cT,to_double(lp.Xb),to_double(lp.Xb[0][2])
    #spec=[(s[0], s[1], to_double(s[2])) for s in lp.Xb]
    t3 = datetime.datetime.now()
    tm=(t3-t0).seconds+((t3-t0).microseconds*1e-6)
    log.debug2('about to leave findbounds')
    res=bound,cT,lp.Xb[0][2],lp.Xb,to_double(lp.curSol()),tm,cblen,lp.iter
    
   # lp.__finalize__() # Carlos edit
    #return Xbs

