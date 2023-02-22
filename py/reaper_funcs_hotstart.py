import setpath
import numpy as np
import os
import sys
import datetime
import pickle

import config
import utils.logging as log
import table_funcs as tf

from solver.v2.mpfr_array.mpfr_array import to_pf_array
from solver.v2.mpfr_array.mpfr_array import array
from solver.v2.mpfr_array.mpfr_array import copy as copy_mpfr_array


# defines lp_problem, mpfr_array, prec_float, etc...specific to version
from version import  Spectrum, lp_problem, NP_PREC, mpfr_array, prec_float, to_double, TD
from importlib import reload # Carlos add  
reload(lp_problem) 

from probspec import dim, opedim, opespin, s0, s1, sdelta, threshold, epsmin


# def setlogfile(logfile):
#     config.logfile=logfile

# Some tag-holders for files
resumefile="resumefile" 
initspecfile = 'initspecfile' 
statusfile = 'statusfile' 
nmax=config.nmax
niter=config.iters 


def pf(x):
    """ Converts a string or a float into a big float 
    input=x: String or Float
    output= prec_float : mpfr float with precision prec """
    
    return prec_float.prec_float(x,prec=212)

# remove the log of previous run if running multiples runs the same day
#os.remove(log.logfile+ str(datetime.date.today())) 

# write the current LP status to a log file keeping a single previous copy
def write_status(lp, jobdir):
    if os.path.isfile(jobdir+statusfile):
        os.rename(jobdir+statusfile, jobdir+statusfile+'.prev')
    #statfile=open(jobdir+statusfile, 'wb', 0)
    statfile=open(jobdir+statusfile, 'w') # Carlos edit
    statfile.write(lp.get_status())
    statfile.close() 
    
def write_resume(jobdir,lp,spectrum):
    log.debug('Writing resume file')
    hd_pf = Hotstart_Data_to_pf(lp)
    pickle.dump((hd_pf,spectrum), open(jobdir+resumefile,'wb'))


def get_resume(jobdir):
    if os.path.isfile(jobdir+resumefile):
        log.debug('Found resume file.')
        (hd_pf_resumed,spectrum) = pickle.load(open(jobdir+resumefile,'rb'))
        hotstartdata = Hotstart_Data_to_ndarray(hd_pf_resumed)
        return hotstartdata,spectrum
    else:
        return None
    
class Hotstart_Data_to_pf:
    """ convert class Hotstart_Data to prec_float """
    def __init__(self, lp):
        self.Cb = to_pf_array(lp.Cb)
        self.Xb = [x for x in lp.Xb]
        self.toElim = lp.toElim
        # THIS SHOULD BE IMPROVED. WHERE ARE THE ARRAY DIM STORED?
        #self.AbT = [to_pf_array(lp.AbT[i]) for i in range(lp.AbT.size/lp.AbT[0].size)]
        self.AbT = [to_pf_array(lp.AbT[i]) for i in range(int(lp.AbT.size/lp.AbT[0].size))] # Carlos edit

#
# Convert back pickled hostart data to proper hotstart data.
#
class Hotstart_Data_to_ndarray:
    """convert class Hotstart_Data_pf from prec_float back to mpfr_array """
    def __init__(self, hotstart):
        self.Cb = array(hotstart.Cb)
        self.Xb = [x for x in hotstart.Xb]
        self.toElim = hotstart.toElim
        self.AbT = array(hotstart.AbT)
                                        

def find_feasible(emin,jobdir, lptab, sig):
    log.debug2('starting feasible vector search')
    log.stats('starting feasable vector search')
    t0 = datetime.datetime.now()
    # emin=epsmin
    #emin=pf(epsmin)
    while emin > pf(0):  
        spectrum = Spectrum(spacetimedim = pf(lptab.spacetimedim),
                    deltamax = pf(lptab.lmax+1.5),
                    lmax = lptab.lmax,
                    de = pf(emin))
        log.debug2('about to call LP_problem')
        lp = lp_problem.LP_Problem(spectrum, lptab, parallel =
                config.point_parallel, pool_size = config.point_poolsize,
                useLU = True)
        
        #set trial spectrum
        trial_spec = GFF_spectrum(sig,int((nmax+1)*(nmax+2)/2),lptab.lmax,
                                  lambda l,delta:delta-0.0*l)

        # if we find an initial spectrum we overwrite as much of the GFF as we
        # can with data from the initial spectrum
        if os.path.isfile(jobdir+initspecfile):
            log.info('Found initial spec file.')
            spec0 = pickle.load(open(jobdir+initspecfile,'rb'))
            initspec=[(s[0], pf(s[1])) for s in spec0]
            GFFpad=True
            if GFFpad:
                usespec=replace_GFF(trial_spec, initspec,
                                  lambda l,delta:delta-pf(0.0*l))
            else:
                usespec=initspec
        else:
            log.info('No spec file -- hotstarting using GFF.')
            usespec=trial_spec
        # SR: the relative coefficient between delta and l can be played with.
        # Say delta-0.5*l gives more weight to operators near the unitarity bound
        # however I checked that 0.0 gives the smallest number of resulting negative coefficients
        lp.set_trial_spectrum(usespec) 
        #lp.set_trial_spectrum_old(usespec) # Carlos temporal
        log.debug('Initialized Xb vector to: %s'%(lp.Xb))
        
        # find a feasable vector
        while lp.is_still_running():
            print(f"lp.is_still_running")
            lp.fewsteps(niter, status_after=True)
            log.debug2(f'done {niter} iterations. (mrc,inv) times ({lp.mrctime}, {lp.invtime})')
            write_status(lp, jobdir)
        log.debug("status: %s"%lp.STATUS)
        log.debug2('Xb: '+str(lp.Xb))
        log.debug2('Cb: '+str(lp.Cb))
        log.debug2('toElim: %d'%lp.toElim)
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
    t1 = datetime.datetime.now()
    tm=(t1-t0).seconds+((t1-t0).microseconds*1e-6)
    log.debug('found feasable vector after %d iterations in %f sec.'%(lp.iter,
        tm))
    log.stats('after found feasible vector')
    return lp,spectrum
#################################################################################

# Run simplex algorithm until the cost function stop changing less than a given threshold
def findbound(sig, jobdir):
    # emin=epsmin
    emin=pf(epsmin) 
    log.proc_id=str(sig)
    log.debug2('entered findbound')
    t0 = datetime.datetime.now()
    log.stats('before get table')
    cblen,lptab=tf.gettable(sig)
    log.debug2('Tables loaded, cblen: %d'%cblen)
    log.stats('after get table')
    # ------ feasible vector    
    log.debug2('starting feasible vector search')
    log.stats('starting feasable vector search')
    resdata=get_resume(jobdir)
    # if there is no resume data (so this is a fresh job) we look for a feasible
    # vector
    if resdata is None:
        log.info('No resume data found: about to search for feasible vector')
        (lp,spectrum)=find_feasible(emin,jobdir,lptab,sig)
        hotstartdata=lp_problem.Hotstart_Data(lp)
        write_resume(jobdir,lp,spectrum)
    else:
        log.info('Found resume data: skipping feasible vector search.')
        log.info('Hotstarting from previous run')
        (hotstartdata,spectrum)=resdata
        lp = lp_problem.LP_Problem(spectrum, lptab, parallel =
                config.point_parallel, pool_size = config.point_poolsize,
                useLU = True)

    lp.reset()
    lp.hotstart(spectrum, hotstartdata)
    log.debug('Value for XB (feasibility or previous run): %s'%str(lp.Xb))
    lp.opebound(opedim,opespin)
    log.debug2('starting optimization')
    oldbnd=None
    bnd=None
    t1 = datetime.datetime.now()
    while lp.is_still_running():
        lp.fewsteps(niter, status_after=False)
        log.debug2(f'done {niter} iterations. (mrc,inv) times ({lp.mrctime}, {lp.invtime})')
       # log.simlog(jobdir+'/simlog',lp.debuginfo) # Carlos edit
        # save status and pause file
        write_status(lp, jobdir)
        write_resume(jobdir,lp,spectrum)
        t2 = datetime.datetime.now()
        tm=(t2-t0).seconds+((t2-t0).microseconds*1e-6)
        # if we have a maxruntime set check and respect it
        if config.maxruntime > 0 and tm > config.maxruntime:
            log.warning('Job exceeded max runtime. About to exit')
            log.info('Removing runfile')
            os.remove(jobdir+runfile)
            log.info('Exiting.')
            sys.exit(0)
        bind=[i for i in range(len(lp.Xb)) if lp.Xb[i][0]=='BOUND']
        if len(bind) > 0:
            bnd= lp.curSol()[bind[0]]
            log.debug('CURRENT BOUND ('+str(sig)+'):'+str(bnd))
            log.debug('CURRENT Cb:'+str(lp.Cb[bind[0]]))
            log.debug('OLD BOUND:'+str(oldbnd))
        if oldbnd is not None and\
            abs(to_double(oldbnd-bnd))/to_double(bnd) < threshold:
            log.debug('~~~~~~~~bound not varying much.  STOPPING')
            break
        if oldbnd is not None:
            log.debug('Bound change ratio: '+str(abs(to_double(oldbnd)-to_double(bnd))/to_double(bnd)))
        if bnd is not None:
            log.debug2('copying NP_PREC')
            oldbnd=NP_PREC.copy(bnd)
        log.stats('bound optimization (%d iterations)'%lp.iter)
    log.debug("status:"+str(lp.STATUS))
    log.debug("leading scalar (ope maxed): "+str(lp.Xb[0][2]))
    log.log_per_sig("leading scalar (ope maxed): "+str(lp.Xb[0][2]), sig)
    log.debug("leading scalar (feasible sol): "+str(lp.Xb[0][2]))
    t2 = datetime.datetime.now()
    tm=(t2-t1).seconds+((t2-t1).microseconds*1e-6)
    log.debug('done maximizing bound after %d iterations in %f sec.'%(lp.iter,
        tm))
    bind=[i for i in range(len(lp.Xb)) if lp.Xb[i][0]=='BOUND']
    if len(bind) > 0:
        bound=(sig, to_double(lp.curSol()[bind[0]]), lp.STATUS,
                to_double(lp.curCost())) 
        cT=(sig, (sig**2 * 4**dim * dim)/((dim-1)*TD(lp.curSol()[bind[0]])) ) 
        log.debug('BOUND: '+str(bound))
        log.debug('Cb:'+str(lp.Cb[bind[0]]))
    else:
        bound=(sig, 0)
        cT=(sig, 0)
        log.debug('***NO BOUND FOUND!!!')
    t3 = datetime.datetime.now()
    tm=(t3-t0).seconds+((t3-t0).microseconds*1e-6)
    log.debug2('about to leave findbounds')
    lp.__finalize__()
    #res=bound,cT,lp.Xb[0][2],lp.Xb,to_double(lp.curSol()),tm,cblen,lp.iter
    #return [to_double(r) for r in res]

def GFF_spectrum(sig, needed, lmax, ordering):
    """ ordering must be a function depending on l and delta, eg lambda l,delta: l+delta """
    def compare(tup):
        return ordering(tup[0],tup[1])
    
    spec0=[]
    l=0
    n=0
    while n<1000:# large number, way larger than the number of operators that we need
        spec0+=[(l,2*sig+l+n)]
        l+=2
        if l>lmax:
            l=0
            n+=2
            
    spec0.sort(key=compare)
    return [(x[0],pf(x[1])) for x in spec0[:needed]]


def replace_GFF(gff, initspec, ordering):
    def compare(tup):
        return ordering(tup[0],tup[1])
    gff.sort(key=compare)
    initspec.sort(key=compare)
    spec=initspec
    log.debug('Padding %d of %d operators using GFF spec.'%(len(gff) -len(initspec), 
             len(gff)))
    spec += gff[len(initspec):len(gff)]
    return spec


