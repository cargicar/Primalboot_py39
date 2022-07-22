#S.Rychkov, S. El-Showk May-July 2013
# Main linear programming class


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
    return prec_float(x,prec=212)

from multiprocessing import Process, Queue, Value, Array
#from multiprocess import Process, Queue, Value, Array # Carlos edit
import timeit 
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
# equivalent to def search pg 640 POOP book
def handleJob(jn, inQueue, outQueue, vecfunlist):
    vecfunlist_local = vecfunlist[:] # create a local copy of vecfunlist
    bbprob = None
    fun=[None for l0 in range(100)]
    ################################################## Carlos edit
    while True:
        if (var := inQueue.get()) is None:
            break
        [l, d0_str, d1_str, rho_str] = var
        
        # Ugly solution but seems to work
        # slicing the string from position 2 to -1 cut out the annyoing b' ' notation
        #print("rho_str %s\n" %(rho_str[0][2:-1]))
        d0=pf(d0_str[2:-1])
        d1=pf(d1_str[2:-1])
#         d0=struct.unpack('f', d0_str)
#         d1=struct.unpack('f', d1_str)
        rho_list=[pf(rho_item[2:-1]) for rho_item in rho_str] # Carlos edit
        rho=mpfr_array.array(rho_list) 
        ################################################## Carlos edit
        
        if fun[l]==None:
            fun[l] = cb.rhomult(rho, vecfunlist_local[l])
        else:
            cb.rhomult_update(rho,vecfunlist_local[l], fun[l])
            
            ####################################################### Carlos edit
        # funlist = [cb.rhomult(self.rho,self.lp_table.vecfunlist[l]) if l%2 == 0 else None for l in range(self.lp_table.lmax+1)]
        
        #funlist_c=[fun[l] if l0==l else None for l0 in range(100)]
        # Carlos edit here range has been readed manually from self.lp_table.lmax+1
        funlist_c =[fun[l] if l%2 == 0 else None for l in range(21)] 
        #funlist_c=list(filter((None).__ne__, funlist_c))
        
        if bbprob==None:
            #bbprob = Divide([fun[l] if l0==l else None for l0 in range(100)], [spectrum.SpecInt(l,d0,d1)])
            bbprob = Divide(funlist_c, [spectrum.SpecInt(l,d0,d1)]) # Carlos edit
        else:
            #bbprob.funlist = [fun[l] if l0==l else None for l0 in range(100)]
            bbprob.funlist = funlist_c
            bbprob.spectrum = [spectrum.SpecInt(l,d0,d1)]
            bbprob.reset()

            ####################################################### 
#             spec= [spectrum.SpecInt(l,d0,d1)]
#             Imin = bbprob.findmin()
#             print("spectrum.SpecInt(l,d0,d1) %s" %spectrum.SpecInt(l,d0,d1))
#             #print("Imin_MRC : %s, d_MRC : %s, l_MRC : %s" %(Imin.fmin, Imin.xmin, Imin.l))
         
        Imin = bbprob.findmin()
        # print("(Imin_MRC : %s, d_MRC : %s, l_MRC : %s) \n" %(Imin.fmin, Imin.xmin, Imin.l))
        
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
    def opebound(self, delta, l):
        self.opedelta=pf(delta)
        self.opel=l
        for i in range(len(self.Xb)):
            if self.Xb[i][1]==self.opel and self.Xb[i][2]==self.opedelta:
                log.debug("found op in spectrum, maximising OPE for (%d, %s)."%(self.opel, 
                    self.opedelta))
                self.Xb[i]=('BOUND', self.opel, self.opedelta)
                self.Cb[i] = "-1"
        self.toElim +=1
        return
############################################################ Carlos edit ################  
# 
#                         PARALLEL BLOCK
############################################################ Carlos edit ################ 
#--------------------------------------------    # Carlos edit    
#     def __finalize__(self):
#         if self.parallel:
#             for t in self.threadPool:
#                 t.terminate()
#--------------------------------------------   
############################################################ Carlos edit ################ 
    # Equivalent to def setup search in pg 642 POOP book
    def startThreads(self):
        if len(self.threadPool) == 0:
            # workers pool
            self.jobQueue = Queue()
            self.respQueue = Queue()
            self.threadPool = [Process(target=handleJob, 
                #args=(i, self.jobQueue, self.respQueue)) 
                args=(i, self.jobQueue, self.respQueue,
                      self.lp_table.vecfunlist)) 
                for i in range(self.pool_size)]
            for p in self.threadPool:
                p.start()          

#---------------------------------------------
############################################################ Carlos edit ################
    def teardown_jobs(self):
        self.rholist = self.rho.tolist1()#buf = self.bufferrho) #[x.__repr__() for x in self.rho ]
        # Signal process termination
        for p in self.spectrum.ilist:
            self.jobQueue.put(None)   
        # join threats  
        for p in self.threadPool:
                p.join()      
        
    # Equivalent to def search in pg 642 POOP book
    def findMRC_parallel(self):
        #print("self.lp_table.lmax+1 %s" %self.lp_table.lmax)
        self.rholist = self.rho.tolist1()#buf = self.bufferrho) #[x.__repr__() for x in self.rho ]
        
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
            
############################################################ Carlos edit ################          
#                  END PARALLEL BLOCK
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
        self.STATUS = STATUS_STILL_RUNNING
        
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
        
       
#--------------------------------------------
    def fewsteps(self,n=1,sort_after=True,status_after=True):
        #count_step=0
        count_step=1 # Carlos edit There was a division by zero at line 303
        madds0,mmults0=stats.ndarray_adds,stats.ndarray_mults
        adds0,mults0=stats.pf_adds,stats.pf_mults
        # NS = Note on Simplex method.pdf
        while count_step < n:    
            if self.toElim == 0:
                self.STATUS = STATUS_AUX_ELIMINATED
                break
            
            self.set_rho() # create the an instance to compute cb.dot.AbInv as per eq 12 in NS
                     # name should be changed not to be confused with cross ratio rho
            self.findMRC() # create an instance to compute the minimum cost Ca, eq 12 NS
            if self.parallel: # Carlos edit
                self.teardown_jobs()
                
            self.costs += [self.curCost()]
            self.redcosts += [self.MRC]
            
            if self.MRC >= - pf("1.0e-60"):#If minimun of the cost is positive, break
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
           # self.sort_basis(False)
            self.sort_basis(True) # Carlos edit
        if status_after:
            self.status()
        return self.Xb
            
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

#     def list_to_doubles(lis):              # Carlos edit
#         l_out=[to_double(x) for x in lis]
#         return l_out
#------------------------------------
    #def rotate_plane_test(self):
    #    # step 1 identify vectors where reduced cost function vanishes or has local negative minimum
    #    # list includes: 1) all vectors  at the unitarity bound (if they are in Xb)
    #    # 2) all vectors at deltamax for each spin
    #    # 3) all vectors where reduced cost function has a negative local minimum
    #    self.set_rho()
    #    
    #    if self.bbprob1 == None: 
    #        funlist = [cb.rhomult(self.rho,self.lp_table.vecfunlist[l]) if l%2 == 0 else None for l in range(self.lp_table.lmax+1)]
    #        self.bbprob1 = Divide(funlist, self.spectrum.ilist)
    #    else:
    #        for l in range(self.lp_table.lmax+1):
    #            if l%2==0:
    #                cb.rhomult_update(self.rho,self.lp_table.vecfunlist[l], self.bbprob1.funlist[l])
    #        self.bbprob1.reset()
    #        
    #    self.bbprob1.find_negative_local_minima()
    #    
    #    if self.CBlen - (len(self.bbprob1.ilist2) + len(self.spectrum.ilist)) >= 0:
    #        print "can try"
    #        self.rotW = mpfr_array.empty((self.CBlen,self.CBlen))
    #        self.rotWinv = mpfr_array.empty((self.CBlen,self.CBlen))
    #        
    #        for i in range(len(self.bbprob1.ilist2)):
    #            I = self.bbprob1.ilist2[i]
    #            for j in range(self.CBlen):
    #                self.rotW[i][j] = self.lp_table.vecfunlist[I.l].value(I.xmin)[j]
    #        
    #        for i in range(len(self.spectrum.ilist)):
    #            I = self.spectrum.ilist[i]
    #            for j in range(self.CBlen):
    #                self.rotW[len(self.bbprob1.ilist2)+i][j] = self.lp_table.vecfunlist[I.l].value(I.d1)[j]
    #            
    #        for i in range(self.CBlen - (len(self.bbprob1.ilist2) + len(self.spectrum.ilist))):
    #            I = self.spectrum.ilist[i]
    #            for j in range(self.CBlen): 
    #                self.rotW[(len(self.bbprob1.ilist2) + len(self.spectrum.ilist)) + i][j] = )
    #                    self.lp_table.vecfunlist[I.l].value(I.d0+pf("0.1"))[j])
    #        
    #        self.rotRHS0 = [None for i in range(self.CBlen)]
    #        
    #        for i in range(len(self.bbprob1.ilist2)):
    #            I = self.bbprob1.ilist2[i]
    #            if I.fmin < pf(0):
    #                self.rotRHS0[i] = -pf(2) * I.fmin
    #            else:
    #                self.rotRHS0[i] = pf(0)
    #        
    #        for i in range(len(self.spectrum.ilist),self.CBlen):
    #            self.rotRHS0[i] = pf(0)
    #            
    #        self.rotRHS = mpfr_array.array(self.rotRHS0)
    #        
    #        inv_assign(self.rotW, self.rotWinv, self.scratch1)   
    #        
    #        self.deltarho = mpfr_array.dot(self.rotRHS, self.rotWinv)
    #        
    #        print "cost=", mpfr_array.dot(self.rho,self.b)
    #        print "deltacost=", mpfr_array.dot(self.deltarho,self.b)
    #        
    #        print "done"
#------------------------------------
    def set_inverse(self):
        if self.useLU:
            mpfr_array.lu_factor_assign(self.AbT, self.LUdata)    
        else:
            inv_assign(self.AbT,self.AbTInv, self.scratch)
#------------------------------------
    def timesAbInv(self, vector):
        if self.useLU:
            return lu_solve(self.LUdata, vector, trans = 1)
            #print "x",max(x-mpfr_array.dot(vector, self.AbTInv))
            #return x
        else:
            return mpfr_array.dot(vector, self.AbTInv)          
#------------------------------------
    def set_rho(self):
        if self.useLU:
            self.rho = - lu_solve(self.LUdata, self.Cb, trans = 0)
            #print "rho",max(self.rho+mpfr_array.dot(self.AbTInv, self.Cb))
        else:
            self.rho = - mpfr_array.dot(self.AbTInv, self.Cb)
#--------------------------------------------
    def curSol(self): 
        return self.timesAbInv(self.b)
#---------------------------------------------
    def curCost(self):
        return mpfr_array.to_pf(mpfr_array.dot(self.Cb, self.curSol()))
#--------------------------------------------
    def status(self):
        print ("Current cost: %s" %( self.curCost()  )) # Carlos edit
        print ("AUX variables:%s" %( self.toElim ))              # Carlos edit
        #print ("Xb: %s" %( self.Xb  ) )                    # Carlos edit
        print ("Solution: %s" %( to_double(self.curSol() ) ) )   # Carlos edit
        print ("Iterations: %s" %( self.iter  ))                # Carlos edit
   

#--------------------------------------------
    def sort_basis(self,status_after=True):
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



