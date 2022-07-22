#from __future__ import print_function
#import datetime
import config
import utils.logging as log

if config.precision:
    import solver.v1.mpfr_array.mpfr_array as NP_PREC
    from solver.v1.mpfr_array.mpfr_array import inv, inv_assign, lu_factor, lu_solve
    from solver.v1.mpfr_array.mpfr_array import to_double
    inf = NP_PREC.array(['inf']);
else:
    import numpy as NP_PREC
    from numpy.linalg import inv  
    from scipy.linalg import lu_solve, lu_factor
    def inv_assign(x,y): y[:] = inv(x)
    def to_double (x): return x
    inf = 'inf'

import numpy 
from scipy import interpolate
import scipy

from multiprocessing import Process, Queue, Value
import timeit
import sys

from utils.timing import t0, deltat
import contract_minimize
reload (contract_minimize)
from contract_minimize import contract_minimize_rel


#    import pdb; pdb.set_trace();

# function called by indiviual thread
def handleJob(jn, inQueue, outQueue, shared_cur_best, vecfunlist):
    vecfunlist_local = vecfunlist[:] # create a copy of vecfunlist
    while True:
        [d0, d1, l, rho0] = inQueue.get()
        rho=NP_PREC.array(rho0);# rho0 is a python double list for precicion=False
        #it is a python list of strings for precision = True 
        d, redcost = contract_minimize_rel(d0, d1, rho, vecfunlist_local[l],
                                            threshold=shared_cur_best.value,
                                            greed_factor=2.0,
                                            relative_diff_goal=0.01, 
                                            debug=False)
        shared_cur_best.value = min(shared_cur_best.value, redcost)
        outQueue.put( [redcost, l, d] )

class Hotstart_Data:
    def __init__(self, lp):
        self.Cb = NP_PREC.copy(lp.Cb)
        self.CbAux = NP_PREC.copy(lp.CbAux)
        self.Xb = [x for x in lp.Xb]
        self.toElim = lp.toElim
        self.AbT = NP_PREC.copy(lp.AbT)

#############################################
#    LP_Problem - main linear programming problem class
#############################################

STATUS_COST_MINIMIZED='Cost minimized'
STATUS_AUX_ELIMINATED='AUX variables eliminated'
STATUS_UNBOUNDED='Problem unbounded'
STATUS_STILL_RUNNING='LP still in progress'
STATUS_COSTFAILS='Too many cost increases'
STATUS_STALLED='LP stalled (cyclic solution)'

class LP_Problem:
    def __init__(self, spectrum, lp_table, parallel = True, pool_size = 8, useLU = True):      
        self.spectrum = spectrum
        self.lp_table = lp_table
        if self.spectrum.spacetimedim != self.lp_table.spacetimedim:
            sys.exit("LP_Problem: problem and table have different spacetimedim")
        self.CBlen = self.lp_table.CBlen
        T2 = NP_PREC.array(numpy.array(
                [interpolate.splev(lp_table.spacetimedim, self.lp_table.funtab[2][i])
                    for i in range(self.CBlen)]))  #here we are first interpolating, then multiplying
        self.b = - self.lp_table.unitCB
        #self.b = - T2#self.lp_table.unitCB
            
        self.parallel = parallel  
        self.useLU = useLU
        
        if self.parallel:
            self.findMRC = self.findMRC_rel_par #MRC stands for 'minimal reduced cost'
            #self.threadPool = []
            self.pool_size = pool_size
            #self.startThreads()
            self.shared_cur_best = Value ('d',0.0)
        else:
            self.findMRC = self.findMRC_rel
        self.reset()
        # used for bounds checking if we're doing that
        self.opedelta=-1 
        self.opel=-1 
        #self.opebound(2, 2)
#--------------------------------------------        
    def __finalize__(self):
        if self.parallel:
            endThreads()
            del self.shared_cur_best
#--------------------------------------------
    def psol(self):
        print 'cost:',self.curCost()
        print 'aux cost:',self.curAuxCost()
        print 'sol:',to_double(self.curSol())

    def opebound(self, delta, l):
        self.opedelta=delta
        self.opel=l
        for i in range(len(self.Xb)):
            if self.Xb[i][1]==self.opel and self.Xb[i][2]==self.opedelta:
                self.Xb[i]=('BOUND', self.opel, self.opedelta)
                self.Cb[i] = "-1"
        self.toElim +=1
        #print 'pre-ope-cost:',self.curCost()
        #print 'pre-ope-aux-cost:',self.curAuxCost()
        #print 'pre-ope-sol:',to_double(self.curSol())
        return

#        self.AIb = self.timesAbInv(self.b)
#        self.Aa = NP_PREC.array(numpy.array(
#            [interpolate.splev(delta, self.lp_table.funtab[l][i])
#                for i in range(self.CBlen)]))  #here we are first interpolating, then multiplying
#                    # at high precision; the procedure must be stable to such interchanges if the
#                    #method is to make sense
#        self.AIa = self.timesAbInv(self.Aa)           
#        
#        zero = NP_PREC.zeros((1,));
#        #inf = NP_PREC.array(['inf']);
#        #log.debug('xcr vals: ' +
#        #    str([x if x > zero else inf for x in self.AIb/self.AIa ]))
#        self.pivot, self.xcr = min(enumerate(
#            [x if x > zero else inf for x in self.AIb/self.AIa ]
#            ), key = lambda x: x[1])
#           
#        if self.xcr == inf:
#            self.STATUS = STATUS_UNBOUNDED
#            print 'UNBOUNDED NO SOLUTION!'
#        # perform basis exchange
#        self.AbT[self.pivot] = self.Aa
#        if self.Xb[self.pivot][0] == 'AUX': 
#            self.toElim -= 1
#        self.Xb[self.pivot] = ('BOUND', l, delta)
#        if to_double(self.b[self.pivot]) > 0:
#            self.Cb[self.pivot] = -NP_PREC.ones((1,))/self.b[self.pivot]
#        else:
#            self.Cb[self.pivot] = NP_PREC.ones((1,))/self.b[self.pivot]
#        self.CbAux[self.pivot] = "0"
#        self.opepivot=self.pivot
#        self.set_inverse()
#        #print 'BOUND pivot',self.pivot
#        #self.AbT[0] = NP_PREC.array(numpy.array(
#        #        [interpolate.splev(delta, self.lp_table.funtab[l][i])
#        #            for i in range(self.CBlen)]))  #here we are first interpolating, then multiplying
#        #self.toElim -= 1
#        #self.Xb[0] = ('OPE', l, delta)
#        #self.Cb[0] = "0"#- self.Cb[0]
#        #self.set_inverse()
#        #print 'post-ope-cost:',self.curCost()
#        #print 'post-ope-aux-cost:',self.curAuxCost()
#        #print 'post-ope-sol:',to_double(self.curSol())
#
#--------------------------------------------
    def startThreads(self):
        self.threadPool = []
        if len(self.threadPool) == 0:
            self.jobQueue = Queue()
            self.respQueue = Queue()
            #self.shared_cur_best = Value ('d',cur_best)
            self.threadPool = [Process(target=handleJob, 
                #args=(i, self.jobQueue, self.respQueue)) 
                args=(i, self.jobQueue, self.respQueue, self.shared_cur_best, self.lp_table.vecfunlist)) 
                for i in range(self.pool_size)]
            for p in self.threadPool:
                p.start()
#--------------------------------------------        
    def endThreads(self):
        for t in self.threadPool:
            t.terminate()
#--------------------------------------------                        
    def reset(self):      
        self.Cb = NP_PREC.ones((self.CBlen,)) / self.b
        #print 'initial CB',self.Cb
        # defined to only have non-zero values for AUX vars (relevant when doing
        # bounds)
        self.CbAux = NP_PREC.ones((self.CBlen,)) / self.b
       
        self.Xb = [("AUX", i) for i in range(self.CBlen)] 
        self.toElim = self.CBlen
        
        self.AbT = NP_PREC.empty((self.CBlen,self.CBlen))
        #self.AbInv = NP_PREC.empty((self.CBlen,self.CBlen))
        self.AbTInv = NP_PREC.empty((self.CBlen,self.CBlen))
        for i in range(self.CBlen): 
            for j in range(self.CBlen):
                if i==j:
                    self.AbT[i][j] = "1"
                    #self.AbInv[i][j] = "1"
                    self.AbTInv[i][j] = "1"
                else:
                    self.AbT[i][j] = "0"
                    #self.AbInv[i][j] = "0"
                    self.AbTInv[i][j] = "0"
                    
        if self.useLU:
            self.set_inverse()
        self.iter=0
        self.STATUS = STATUS_STILL_RUNNING
#--------------------------------------------
    def hotstart(self, spectrum, hotstart_data):
        self.spectrum = spectrum
        self.Cb = NP_PREC.copy(hotstart_data.Cb)
        self.CbAux = NP_PREC.copy(hotstart_data.CbAux)
        self.Xb[:] = hotstart_data.Xb
        self.toElim = hotstart_data.toElim
        self.AbT = NP_PREC.copy(hotstart_data.AbT)
        self.set_inverse()
        self.iter=0
        self.STATUS = STATUS_STILL_RUNNING      
#--------------------------------------------       
    def is_still_running(self):
        return (self.STATUS == STATUS_STILL_RUNNING)
#---------------------------------------------           
    def findMRC_rel(self, threshold = 0.0): 
        self.set_rho()
        
        self.RClist = []
        cur_best = threshold
        for p in self.spectrum.ilist:
            results = contract_minimize_rel(p.d0, p.d1, 
                                            self.rho, self.lp_table.vecfunlist[p.l],
                                            threshold = cur_best,
                                            greed_factor=2.0,
                                            relative_diff_goal=0.01,
                                            debug=True)
            d, redcost = results[0:2]
            
            if len(results)>2: # if debug=True
                bbprob = results[2]
            else:
                bbprob = None
                
            self.RClist.append([redcost, p.l, d, bbprob])
            cur_best = min(cur_best, redcost)
            
        self.MRC, self.l_MRC, self.d_MRC = min(self.RClist, key = lambda x: x[0] )[0:3]
#---------------------------------------------
    def findMRC_rel_par(self, threshold = 0.0):
        self.set_rho()
        self.rholist = self.rho.tolist()
        
        self.shared_cur_best.value = threshold
        for p in self.spectrum.ilist:
            self.jobQueue.put([p.d0, p.d1, p.l, self.rholist])
        
        self.RClist = []
        expected = len(self.spectrum.ilist)
        while len(self.RClist) < expected:
            self.RClist.append(self.respQueue.get())
        
        self.MRC, self.l_MRC, self.d_MRC = min(self.RClist, key = lambda x: x[0] )[0:3]
#--------------------------------------------
    def fewsteps(self,n=1,sort_after=True,status_after=True):
        count_step=0
        if self.parallel:
            self.startThreads()
        while count_step < n:
            if self.toElim == 0:
                self.STATUS = STATUS_AUX_ELIMINATED
                break
            
            #print self.curCost();
            #threshold = - 1.0e-04 * to_double(self.curCost())
            #threshold = - config.threshold_multiplier * to_double(self.curCost())
            threshold = - config.threshold_multiplier * to_double(self.curAuxCost())
            self.findMRC(threshold = threshold)
            
            if self.MRC >= threshold:
                self.STATUS = STATUS_COST_MINIMIZED
                break
            # find pivot    
            self.AIb = self.timesAbInv(self.b)
            self.Aa = NP_PREC.array(numpy.array(
                [interpolate.splev(self.d_MRC, self.lp_table.funtab[self.l_MRC][i])
                    for i in range(self.CBlen)]))  #here we are first interpolating, then multiplying
                        # at high precision; the procedure must be stable to such interchanges if the
                        #method is to make sense
            self.AIa = self.timesAbInv(self.Aa)           
            
            zero = NP_PREC.zeros((1,));
            #inf = NP_PREC.array(['inf']);
            #log.debug('xcr vals: ' +
            #    str([x if x > zero else inf for x in self.AIb/self.AIa ]))
            self.pivot, self.xcr = min(enumerate(
                [x if x > zero else inf for x in self.AIb/self.AIa ]
                ), key = lambda x: x[1])
               
            if self.xcr == inf:
                self.STATUS = STATUS_UNBOUNDED
                break
            # perform basis exchange
            self.AbT[self.pivot] = self.Aa
            if self.Xb[self.pivot][0] == 'AUX': 
                self.toElim -= 1
            if self.d_MRC == self.opedelta and self.l_MRC == self.opel:
                self.Xb[self.pivot] = ('BOUND', self.l_MRC, self.d_MRC)
                if to_double(self.b[self.pivot]) > 0:
                    self.Cb[self.pivot] = -NP_PREC.ones((1,))/self.b[self.pivot]
                else:
                    self.Cb[self.pivot] = NP_PREC.ones((1,))/self.b[self.pivot]
                #self.Cb[self.pivot] = -self.Cb[self.pivot]
                self.CbAux[self.pivot] = "0"
            else:
                self.Xb[self.pivot] = ('OPE', self.l_MRC, self.d_MRC)
                self.Cb[self.pivot] = "0"
                self.CbAux[self.pivot] = "0"
            self.set_inverse()
            self.iter += 1
            count_step += 1
        if 'threshold' in locals():
            print 'threshold',threshold, ' (aux cost:',to_double(self.curAuxCost()),')'
        if sort_after:
            self.sort_basis(False)
        if status_after:
            self.status()
        if self.parallel:
            self.endThreads()
#------------------------------------
    def set_inverse(self):
        if self.useLU:
            self.LUdata = lu_factor(self.AbT)   
            #inv_assign(self.AbT,self.AbTInv)  
        else:
            #self.AbTInv = inv(self.AbT) 
            inv_assign(self.AbT,self.AbTInv)
#------------------------------------
    def timesAbInv(self, vector):
        if self.useLU:
            return lu_solve(self.LUdata, vector, trans = 1)
            #print "x",max(x-NP_PREC.dot(vector, self.AbTInv))
            #return x
        else:
            return NP_PREC.dot(vector, self.AbTInv)          
#------------------------------------
    def set_rho(self):
        if self.useLU:
            self.rho = - lu_solve(self.LUdata, self.Cb, trans = 0)
            #print "rho",max(self.rho+NP_PREC.dot(self.AbTInv, self.Cb))
        else:
            self.rho = - NP_PREC.dot(self.AbTInv, self.Cb)
#--------------------------------------------
    def curSol(self): 
        return self.timesAbInv(self.b)
#---------------------------------------------
    def curCost(self):
        return NP_PREC.dot(self.Cb, self.curSol())
#---------------------------------------------
    def curAuxCost(self):
        """Only aux variables"""
        return NP_PREC.dot(self.CbAux, self.curSol())
#--------------------------------------------
    def status(self):
        print "Current cost:", self.curCost()
        print "AUX variables:", self.toElim
        print "Xb:", self.Xb
        print "Solution:", to_double(self.curSol())
        print "Iterations:", self.iter    
#--------------------------------------------
    def sort_basis(self,status_after=True):
        self.ordering = numpy.lexsort( [[x[2] if len(x)>2 else 0 for x in self.Xb],
                                [x[1] for x in self.Xb],
                                [0 if x[0]=='OPE' else 1 for x in self.Xb]])
        self.Xb = [ self.Xb[i] for i in self.ordering]
        self.Cb_temp = NP_PREC.empty_like(self.Cb)
        self.CbAux_temp = NP_PREC.empty_like(self.CbAux)
        self.AbT_temp = NP_PREC.empty_like(self.AbT)
        for i,ii in enumerate(self.ordering):
            self.Cb_temp[i] = self.Cb[ii]
            self.CbAux_temp[i] = self.CbAux[ii]
            self.AbT_temp[i] = self.AbT[ii]
        self.Cb = self.Cb_temp
        self.CbAux = self.CbAux_temp
        self.AbT = self.AbT_temp
        self.set_inverse()
        if status_after:
            self.status()



