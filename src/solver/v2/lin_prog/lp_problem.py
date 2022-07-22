#S.Rychkov, S. El-Showk May-July 2013
# Main linear programming class

# Edited by Carlos Cardona 2022


import solver.v2.mpfr_array.mpfr_array as mpfr_array
from solver.v2.mpfr_array.mpfr_array import inv_assign, lu_factor, lu_solve, to_double
import solver.v2.prec_float.prec_float as PF
from solver.v2.prec_float.prec_float import prec_float
import solver.v2.rho_repr.cb as cb
import solver.v2.lin_prog.spectrum as spectrum
from utils.memory import rss
import utils.logging as log
import utils.stats as stats



def pf(x):
        """ Converts a string or a float into a big float 
    input=x: String or Float
    output= prec_float : mpfr float with precision prec """
    return prec_float(x,prec=212)

# from multiprocessing import Process, Queue, Value, Array
from multiprocess import Process, Value, Array , Queue  # Carlos edit 
#import threading # Carlos edit 
import queue # Carlos edit

import timeit 
import datetime
import time
import sys
import numpy

from utils.timing import t0, deltat
import solver.v2.branch_bound.divide 
from importlib import reload            # Carlos add
reload(solver.v2.branch_bound.divide)
from solver.v2.branch_bound.divide import Divide
import struct

# function called by indiviual thread
def handleJob(jn, inQueue, outQueue, vecfunlist):
    """ Take the spectrum list as well as the conformal blocks in rho coordinates evaluated on such spectrum from inQueue
        and past it to the Divide algorithm to compute the minimun of the cost function for a given spin. Record each mininimun for each spin into outQueue"""
    # Now we need to transform vecfunlist to mpfr array. Figure also how to copy mpfr
    vecfunlist_local = vecfunlist[:] # create a local copy of vecfunlist
    #print("vecfunlist_str %s" %vecfunlist_str)
    #vecfunlist_local=[mpfr_array.array(item) for item in vecfunlist_str]
   

    bbprob = None
    fun=[None for l0 in range(100)]
    while True:
        [l, d0_str, d1_str, rho_str] = inQueue.get()
        d0=pf(d0_str)
        d1=pf(d1_str)
        #print("d1 is equal to : %s \n" %d1_str)
        rho=mpfr_array.array(rho_str)
        
        
        if fun[l]==None:
            fun[l] = cb.rhomult(rho, vecfunlist_local[l])
        else:
            cb.rhomult_update(rho,vecfunlist_local[l], fun[l])
        #funlist_c=[fun[l] if l0==l else l0 for l0 in range(0,100,2)]
        #print("funlist_c %s \n" %funlist_c)
            
        if bbprob==None:
            bbprob = Divide([fun[l] if l0==l else None for l0 in range(100)], [spectrum.SpecInt(l,d0,d1)])
            #bbprob = Divide([fun[l] if l0==l else None for l0 in range(0,100,2)], [spectrum.SpecInt(l,d0,d1)]) # Carlos edit
          
        else:
            bbprob.funlist = [fun[l] if l0==l else None for l0 in range(100)]
            # bbprob.funlist = [fun[l] if l0==l else None for l0 in range(0,100,2)] # Carlos edit
            bbprob.spectrum = [spectrum.SpecInt(l,d0,d1)]
            bbprob.reset()
   
        Imin = bbprob.findmin()
        #print("Imin : %s \n" %Imin )
        
        outQueue.put( [Imin.fmin.__repr__(), l, Imin.xmin.__repr__()] )


class Hotstart_Data:
    def __init__(self, lp):
        self.Cb = mpfr_array.copy(lp.Cb)
        self.Xb = [x for x in lp.Xb]
        self.toElim = lp.toElim
        self.AbT = mpfr_array.copy(lp.AbT)

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
        self.prec=self.lp_table.prec
        self.CBlen = self.lp_table.CBlen
        self.b = - self.lp_table.unitCB
            
        self.useLU = useLU
        self.parallel = parallel
        
        if self.parallel:
            self.findMRC = self.findMRC_parallel #MRC stands for 'minimal reduced cost'
            self.threadPool = []
            self.pool_size = pool_size
            self.startThreads()
        else:
            self.findMRC = self.findMRC_serial
        self.bbprob = None
        self.bbprob1 = None
        self.buffer0=" "*(2*self.prec)
        self.buffer1=" "*(2*self.prec)
        self.bufferrho=[" "*(2*self.prec) for i in range(self.CBlen)]
        # used for bounds checking if we're doing that
        self.opedelta=pf(-1) 
        self.opel=-1 
        self.reset()
#--------------------------------------------        
    def opebound_old(self, delta, l):
        self.opedelta=pf(delta)
        self.opel=l
        print ('Xb: ', self.Xb)
        for i in range(len(self.Xb)):
            if self.Xb[i][1]==self.opel and self.Xb[i][2]==self.opedelta:
                log.debug("found op in spectrum, maximising OPE for (%d, %s)."%(self.opel, 
                    self.opedelta))
                self.Xb[i]=('BOUND', self.opel, self.opedelta)
                self.Cb[i] = "-1"
        self.toElim +=1
        return
#--------------------------------------------        
    def opebound(self, delta, l):
        """ Look for bound operator in feasible basis. If find it, it gets tagged as bound OPE 
        if is not found it, it bringing it into the feasible basis"""
        self.opedelta=pf(delta)
        self.opel=l
        #log.debug('Xb: '+self.Xb.__repr__())
        
        for i in range(len(self.Xb)):
            if self.Xb[i][1]==self.opel and self.Xb[i][2]==self.opedelta:
                log.debug("found op in spectrum, maximising OPE for (%d, %s)."%(self.opel, 
                    self.opedelta))
                self.Xb[i]=('BOUND', self.opel, self.opedelta)
                self.Cb[i] = "-1"
                break
        else:
            log.debug("op (%d, %s) not in the spectrum, bringing it in"%(self.opel, 
                    self.opedelta))
            self.AIb = self.timesAbInv(self.b)
            self.Aa = self.lp_table.vecfunlist[self.opel].value(self.opedelta)
            
            self.AIa = self.timesAbInv(self.Aa)
            
            zero = mpfr_array.zeros((1,));
            inf = mpfr_array.array(['inf']);
            
            self.pivot, self.xcr = min(enumerate(
                [x if x > zero else inf for x in self.AIb/self.AIa ]
                ), key = lambda x: x[1])
            
            if self.xcr == inf:
                self.STATUS = STATUS_UNBOUNDED
                log.debug("ope coefficient unbounded")
                sys.exit()
            
            # perform basis exchange
            self.AbT[self.pivot] = self.Aa
            
            self.Xb[self.pivot] = ('BOUND', self.opel, self.opedelta)
            self.Cb[self.pivot] = "-1"
            self.set_inverse()
        
        self.iter=0
        self.invtime=0
        self.mrctime=0
        self.toElim =1 # there is nothing to eliminate in this case; this is just to make fewsteps happy
        self.STATUS = STATUS_STILL_RUNNING
    
#--------------------------------------------        

#--------------------------------------------        
    def __finalize__(self):
        """ signal end of the algorithm for cleaning purposes """
        if self.parallel:
            for t in self.threadPool:
                t.terminate()
#--------------------------------------------

    def startThreads(self):
        """ Distributes  handleJob function across cpu kernels """
        if len(self.threadPool) == 0:
            self.jobQueue = Queue()
            self.respQueue = Queue()
            self.threadPool = [Process(target=handleJob, 
                args=(i, self.jobQueue, self.respQueue,
                      self.lp_table.vecfunlist))     
                for i in range(self.pool_size)]
            for p in self.threadPool:
                p.start()        
############################################################ Carlos edit ################ 
#--------------------------------------------

    def reset(self):      
        self.Cb = mpfr_array.ones((self.CBlen,)) / self.b
       
        self.Xb = [("AUX", i) for i in range(self.CBlen)] 
        self.toElim = self.CBlen
        
        self.AbT = mpfr_array.empty((self.CBlen,self.CBlen))
        self.LUscratch = mpfr_array.empty((self.CBlen,self.CBlen))
        self.scratch = mpfr_array.empty((self.CBlen,self.CBlen))
        self.scratch1 = mpfr_array.empty((self.CBlen,self.CBlen))
        self.AbTInv = mpfr_array.empty((self.CBlen,self.CBlen))
        
        self.costs = []
        self.redcosts = []
        for i in range(self.CBlen): 
            for j in range(self.CBlen):
                if i==j:
                    self.AbT[i][j] = "1"
                    self.AbTInv[i][j] = "1"
                else:
                    self.AbT[i][j] = "0"
                    self.AbTInv[i][j] = "0"
                    
        if self.useLU:
            self.LUdata = mpfr_array.LUdcmp_data(self.LUscratch) # this LUscratch array cannot be used for
                        #anything else at this point
                         # 
            self.set_inverse()
        self.iter=0
        self.iter=0
        self.invtime=0
        self.mrctime=0
        self.STATUS = STATUS_STILL_RUNNING
#--------------------------------------------        
#################### Taken from cern_new_mrc_log ########################       
    def set_trial_spectrum(self, trial_spec):
        """ trial_spec must be an array of length CBlen with elements (l, delta (pf))"""
        self.reset() 
        #self.Xb = [("OPE",)+op for op in trial_spec]

        for i,op in enumerate(trial_spec):
            self.Xb[i] = ("OPE",)+op 
            self.AbT[i]=self.lp_table.vecfunlist[op[0]].value(op[1])
       
        # invert AbT (or setup LUdata)
        self.set_inverse()

        
        self.Cb = mpfr_array.zeros((self.CBlen,))
        
        # now compute the solution. If some coefficients are negative, flip the sign of
        # the corresponding vectors and set Cb to eliminate these coefficients
        
        #self.toElim=0
        #self.toElim=self.CBlen-len(trial_spec)
        sol = self.curSol()
        solstr=''
        for i in range(self.CBlen):
            if self.Xb[i][0] == 'AUX': 
                solstr += '("AUX", %d, %s),'%(self.Xb[i][1], str(sol[i]))
            else:
                solstr += '(%d, %s, %s),'%(self.Xb[i][1], str(self.Xb[i][2]),
                    str(sol[i]))
        log.debug2('Got solution: [%s]'%solstr)
        zero = mpfr_array.zeros((1,))
        for i in range(self.CBlen):
            if sol[i] < zero:
                self.AbT[i] = - self.AbT[i]
                #self.Xb[i] = ("AUX",self.Xb[i][1],self.Xb[i][2])
                self.Xb[i] = ("AUX", i)
                self.Cb[i] = "1"
                self.toElim+=1
        self.toElim=0
        for i in range(self.CBlen):
            if self.Xb[i][0]=='AUX':
                self.toElim += 1
        log.info('Still %d AUX vars to eliminate after setting trial spectrum.'%self.toElim)
        self.set_inverse()
        #self.status()         
#--------------------------------------------      

    def set_trial_spectrum_old(self, trial_spec):
        """ trial_spec must be an array of length CBlen with elements (l, delta (pf))"""
       
        self.Xb = [("OPE",)+op for op in trial_spec]
        
        for i,op in enumerate(trial_spec):
            self.AbT[i]=self.lp_table.vecfunlist[op[0]].value(op[1])
       
        # invert AbT (or setup LUdata)
        self.set_inverse()

        
        self.Cb = mpfr_array.zeros((self.CBlen,))
        
        # now compute the solution. If some coefficients are negative, flip the sign of
        # the corresponding vectors and set Cb to eliminate these coefficients
        
        self.toElim=0
        sol = self.curSol()
        zero = mpfr_array.zeros((1,))
        for i in range(self.CBlen):
            if sol[i] < zero:
                self.AbT[i] = - self.AbT[i]
                self.Xb[i] = ("AUX",self.Xb[i][1],self.Xb[i][2])
                self.Cb[i] = "1"
                self.toElim+=1
        
        self.set_inverse()
        #self.status()    
########################################################################        
        
#--------------------------------------------
    def hotstart(self, spectrum, hotstart_data):
        self.spectrum = spectrum
        self.Cb = mpfr_array.copy(hotstart_data.Cb)
        self.Xb[:] = hotstart_data.Xb
        self.toElim = hotstart_data.toElim
        self.AbT = mpfr_array.copy(hotstart_data.AbT)
        self.set_inverse()
        self.iter=0
        self.STATUS = STATUS_STILL_RUNNING      
#--------------------------------------------       
    def is_still_running(self):
        return (self.STATUS == STATUS_STILL_RUNNING)
#---------------------------------------------           
    def findMRC_serial(self):
        """ Take the list of cost functions for all spins in spectrum and pass it 
        to the minimization algorithm Divide, out put the global minimun."""
        if self.bbprob == None: 
            funlist = [cb.rhomult(self.rho,self.lp_table.vecfunlist[l]) if l%2 == 0 else None for l in range(self.lp_table.lmax+1)]
            self.bbprob = Divide(funlist, self.spectrum.ilist)
        
        else:
            for l in range(self.lp_table.lmax+1):
                if l%2==0:
                    cb.rhomult_update(self.rho,self.lp_table.vecfunlist[l], self.bbprob.funlist[l])
            self.bbprob.reset()
          
        #funlist1 = [fun1(f) if f != None else None for f in funlist]
        #intlist1 = [ spectrum.SpecInt(I.l, (-I.d1).exp(), (-I.d0).exp()) for I in self.spectrum.ilist]

        Imin = self.bbprob.findmin()
        self.MRC = Imin.fmin
        self.d_MRC = Imin.xmin
        self.l_MRC = Imin.l
        
#---------------------------------------------
    def findMRC_parallel(self):        
        """ Take the list of spectrum quantum numbers (i.e a spin l, 
         the dimension d0 at the unitarity bound for that given spin, 
         and the max dimension p.d1 )
         as well as the list of Conformal block values in rho coordinates for the given elements in the spectrum, 
         and put it into the jobQueue. Elements in jobQueue are handle by function handleJob, which fills
         the respQueue containing the minimun of the cost function for a given element of the spectrum (spin l).
         Finally, the absolute minimun is taking from respQueue and recorded in self.MRC, self.l_MRC, self.d_MRC"""
        
        self.rholist = self.rho.tolist1() #buf = self.bufferrho) #[x.__repr__() for x in self.rho ]
        for p in self.spectrum.ilist:
            self.jobQueue.put([p.l, p.d0.__repr__(),
                               #p.d0.bufrepr(self.buffer0),
                               p.d1.__repr__(),
                               #p.d1.bufrepr(self.buffer1),
                               self.rholist])        
            #print p, p.d0.bufrepr(self.buffer0), p.d1.bufrepr(self.buffer1)
        #a= raw_input("a")
        
        self.RClist = []
        expected = len(self.spectrum.ilist)
        while len(self.RClist) < expected:
            redcost_str, l, dstring = self.respQueue.get()
            self.RClist.append([pf(redcost_str),l,pf(dstring)])
        
        self.MRC, self.l_MRC, self.d_MRC = min(self.RClist, key = lambda x: x[0] )[0:3] 
               
#--------------------------------------------
    def fewsteps(self,n=1,sort_after=True,status_after=True):
        """ Perform n iterations of the linear programming main algorithm """
        #count_step=0
        count_step=1
        madds0,mmults0=stats.ndarray_adds,stats.ndarray_mults
        adds0,mults0=stats.pf_adds,stats.pf_mults
        # NS = Note on Simplex method.pdf
        while count_step < n:    
            if self.toElim == 0:
                self.STATUS = STATUS_AUX_ELIMINATED
                print("STATUS_AUX_ELIMINATED" )
                break
            
            self.set_rho() # create the an instance to compute cb.dot.AbInv as per eq 12 in NS
                     # name should be changed not to be confused with cross ratio rho
            t0=datetime.datetime.now()
            self.findMRC() # create an instance to compute the minimum cost Ca, eq 12 NS
            self.costs += [self.curCost()]
            self.redcosts += [self.MRC]
            t1=datetime.datetime.now()
            self.mrctime += (t1-t0).seconds+((t1-t0).microseconds*1e-6)
            if self.MRC >= - pf("1.0e-60"):#If minimun of the cost is positive, break
                print("self.MRC IS POSITIVE" )
                self.STATUS = STATUS_COST_MINIMIZED

                break
            # find pivot
            self.AIb = self.timesAbInv(self.b) # AbInv.dot.b (numerator eq 13 NS for every i)
 
            # Aa at denominator eq 13 NS. In other words, column vector Aa_i^*
            self.Aa = self.lp_table.vecfunlist[self.l_MRC].value(self.d_MRC) 

            # AbInv.dot.Aa (denominator eq 13 NS for every i)
            self.AIa = self.timesAbInv(self.Aa) 

            zero = mpfr_array.zeros((1,));
            inf = mpfr_array.array(['inf']);
            
            # pick minimum from eq 13 NS and its index 
            # self.pivot, self.xcr =(index, min_value)
            self.pivot, self.xcr = min(enumerate(
                [x if x > zero else inf for x in self.AIb/self.AIa ]
                ), key = lambda x: x[1])
               
            if self.xcr == inf:
                self.STATUS = STATUS_UNBOUNDED
                break
            # perform basis exchange
            self.AbT[self.pivot] = self.Aa # The row self.pivot the AbT is the 
                                           # is the columb vector Aa
            if self.Xb[self.pivot][0] == 'AUX': 
                self.toElim -= 1
            # OPE bound code
            # If OPE to be optimized enter the basic basis, record it.
            if self.d_MRC == self.opedelta and self.l_MRC == self.opel:
                log.debug("op just entered spectrum, maximising OPE for (%d, %s)."%(self.opel, 
                        self.opedelta))
                self.Xb[self.pivot] = ('BOUND', self.l_MRC, self.d_MRC)
                # for  sign(b )=pm  then set cb[out]=-mp/b
                if to_double(self.b[self.pivot]) > 0:
                    self.Cb[self.pivot] = -mpfr_array.ones((1,))/self.b[self.pivot]
                else:
                    self.Cb[self.pivot] = mpfr_array.ones((1,))/self.b[self.pivot]
            # for entering operator different than the to-be optimized, make cb[in]=0
            else:
                self.Xb[self.pivot] = ('OPE', self.l_MRC, self.d_MRC)
                self.Cb[self.pivot] = "0"
            self.set_inverse()
        
            self.iter += 1
            count_step += 1
    
        madds1,mmults1=stats.ndarray_adds,stats.ndarray_mults
        adds1,mults1=stats.pf_adds,stats.pf_mults
        log.debug("number of mpfr (add,mult) in %d iterations: (%d,%d)"%(count_step, 
            madds1-madds0, mmults1-mmults0))
        log.debug("number of pf (add,mult) in %d iterations: (%d,%d)"%(count_step,
            adds1-adds0, mults1-mults0))
        addavg=float(adds1+madds1-adds0-madds0)/count_step
        multavg=float(mults1+mmults1-mults0-mmults0)/count_step
        log.debug("avg (add,mult) per iteration: (%d,%d)."%(addavg,multavg))
        log.stats("fewsteps at %d iterations"%self.iter)
                
        if sort_after:
            self.sort_basis(False)
            #self.sort_basis(True) # Carlos edit
        if status_after:
            self.status()
        #return self.Xb
            
        #self.rotate_plane_test()
        #self.adjust_test()
#------------------------------------
    def adjust_test(self):
        # check if all operator except on the unitarity bound come in close pairs
        pairlist = []
        for i in range(self.CBlen-1):
            op = self.Xb[i]
            if (op[0] == 'AUX'):
                break # stop cycle when reaching AUX vectors
            
            lowerbound = [specint.d0 for specint in self.spectrum.ilist if specint.l == op[1]][0]
            
            if op[2] - lowerbound < pf("1.0e-10"):
                continue # skip operators at the lower end of the unitarity bound
            
            op1 = self.Xb[i+1]
            if (op1[0] == 'OPE' and op1[1] == op[1] and
                op1[2] - op[2] < pf("1.0e-2")):
                pairlist += [(i, op,op1)]
                i+=2
            else:
                flag = "not in pairs"
                break
        else:
            flag = "in pairs"
            print(pairlist)


#------------------------------------
    def set_inverse(self):
        """ Set inverse for AbT matrix """
        if self.useLU:
            mpfr_array.lu_factor_assign(self.AbT, self.LUdata)    
        else:
            inv_assign(self.AbT,self.AbTInv, self.scratch)
#------------------------------------
    def timesAbInv(self, vector):
        """ Multiplies vector times AbInv """
        if self.useLU:
            return lu_solve(self.LUdata, vector, trans = 1)
            #print "x",max(x-mpfr_array.dot(vector, self.AbTInv))
            #return x
        else:
            return mpfr_array.dot(vector, self.AbTInv)          
#------------------------------------
    def set_rho(self):
        """ Computes  AbTInv times Cb"""
        if self.useLU:
            self.rho = - lu_solve(self.LUdata, self.Cb, trans = 0)
            #print "rho",max(self.rho+mpfr_array.dot(self.AbTInv, self.Cb))
        else:
            self.rho = - mpfr_array.dot(self.AbTInv, self.Cb)
#--------------------------------------------
    def curSol(self): 
        """ Computes  AbInv times b"""
        return self.timesAbInv(self.b)
#---------------------------------------------
    def curCost(self):
        """ Current cost function """
        return mpfr_array.to_pf(mpfr_array.dot(self.Cb, self.curSol()))
#--------------------------------------------
    def status(self):
        """ Print status of the simplex algorithm """
        print ("Current cost: %s" %( self.curCost()  )) # Carlos edit
        print ("AUX variables:%s" %( self.toElim ))              # Carlos edit
        #print ("Xb: %s" %( self.Xb  ) )                    # Carlos edit
        print ("Solution: %s" %( to_double(self.curSol() ) ) )   # Carlos edit
        print ("Iterations: %s" %( self.iter  ))                # Carlos edit

#--------------------------------------------
    def get_status_old(self):
        status="Current cost:%s\nAUX variables:%s\nXb:\
                %s\nSolution:%s\nIterations:%s\n"%(str(self.curCost()),
            str(self.toElim), str(self.Xb), str(to_double(self.curSol())), 
            str(self.iter))
        return status
    
    def get_status(self):
        newline='\n'    
        status=f"Current cost: {self.curCost()}{newline}AUX variables: {self.toElim}{newline}Xb: {self.Xb}{newline}Solution: {to_double(self.curSol())}{newline} Iterations:{self.iter}"
        return status
#--------------------------------------------

#--------------------------------------------
    def sort_basis(self,status_after=True):
        """ Sort feasible basis first by spin then by dimension"""
        self.ordering = numpy.lexsort( [[PF.to_double(x[2]) if len(x)>2 else 0 for x in self.Xb],
                                [x[1] for x in self.Xb],
                                [0 if x[0]=='OPE' else 1 for x in self.Xb]])
        self.Xb = [ self.Xb[i] for i in self.ordering]
        self.Cb_temp = mpfr_array.empty_like(self.Cb)
        self.AbT_temp = mpfr_array.empty_like(self.AbT)
        for i,ii in enumerate(self.ordering):
            self.Cb_temp[i] = self.Cb[ii]
            self.AbT_temp[i] = self.AbT[ii]
        self.Cb = self.Cb_temp
        self.AbT = self.AbT_temp
        self.set_inverse()
        if status_after:
            self.status()

