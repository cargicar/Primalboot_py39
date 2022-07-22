#S.Rychkov, S. El-Showk May-July 2013
# Main linear programming class


import solver.v2.mpfr_array.mpfr_array as mpfr_array
from solver.v2.mpfr_array.mpfr_array import inv_assign, lu_factor, lu_solve, to_double
import solver.v2.prec_float.prec_float as PF
from solver.v2.prec_float.prec_float import prec_float
import solver.v2.rho_repr.cb as cb
import spectrum
from utils.memory import rss

def pf(x):
    return prec_float(x,prec=212)

from multiprocessing import Process, Queue, Value, Array
import timeit
import time
import sys
import numpy

from utils.timing import t0, deltat
import solver.v2.branch_bound.divide 
reload(solver.v2.branch_bound.divide)
from solver.v2.branch_bound.divide import Divide

# function called by indiviual thread
def handleJob(jn, inQueue, outQueue, vecfunlist):
    vecfunlist_local = vecfunlist[:] # create a local copy of vecfunlist
    print jn
    bbprob = None
    fun=[None for l0 in range(100)]
    while True:
        [l, d0_str, d1_str, rho_str] = inQueue.get()
        d0=pf(d0_str)
        d1=pf(d1_str)
        rho=mpfr_array.array(rho_str)
        
        if fun[l]==None:
            fun[l] = cb.rhomult(rho, vecfunlist_local[l])
        else:
            cb.rhomult_update(rho,vecfunlist_local[l], fun[l])
            
        if bbprob==None:
            bbprob = Divide([fun[l] if l0==l else None for l0 in range(100)], [spectrum.SpecInt(l,d0,d1)])
        else:
            bbprob.funlist = [fun[l] if l0==l else None for l0 in range(100)]
            bbprob.spectrum = [spectrum.SpecInt(l,d0,d1)]
            bbprob.reset()
   
        Imin = bbprob.findmin()
        
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
        #self.mem = [0,0,0,0,0,0,0,0,0,0,0,0]
        self.buffer0=" "*(2*self.prec)
        self.buffer1=" "*(2*self.prec)
        self.bufferrho=[" "*(2*self.prec) for i in range(self.CBlen)]
        
        self.reset()
#--------------------------------------------        
    def opebound(self, delta, l):
        self.opedelta=delta
        self.opel=l
        for i in range(len(self.Xb)):
            if self.Xb[i][1]==self.opel and self.Xb[i][2]==self.opedelta:
                self.Xb[i]=('BOUND', self.opel, self.opedelta)
                self.Cb[i] = "-1"
        self.toElim +=1
        return
#--------------------------------------------        
    def __finalize__(self):
        if self.parallel:
            for t in self.threadPool:
                t.terminate()
#--------------------------------------------
    def startThreads(self):
        if len(self.threadPool) == 0:
            self.jobQueue = Queue()
            self.respQueue = Queue()
            self.threadPool = [Process(target=handleJob, 
                #args=(i, self.jobQueue, self.respQueue)) 
                args=(i, self.jobQueue, self.respQueue,
                      self.lp_table.vecfunlist)) 
                for i in range(self.pool_size)]
            for p in self.threadPool:
                p.start()
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
            #mem=rss()
            for l in range(self.lp_table.lmax+1):
                if l%2==0:
                    cb.rhomult_update(self.rho,self.lp_table.vecfunlist[l], self.bbprob.funlist[l])
            self.bbprob.reset()
            #m1=rss()
            #self.mem[7] += m1-mem
             
        #funlist1 = [fun1(f) if f != None else None for f in funlist]
        #intlist1 = [ spectrum.SpecInt(I.l, (-I.d1).exp(), (-I.d0).exp()) for I in self.spectrum.ilist]
        
        
        #mem=rss()
        Imin = self.bbprob.findmin()
        #m1=rss()
        #self.mem[8] += m1-mem

        self.MRC = Imin.fmin
        self.d_MRC = Imin.xmin
        self.l_MRC = Imin.l
        
#---------------------------------------------
    def findMRC_parallel(self):
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

        
#--------------------------------------------
    def fewsteps(self,n=1,sort_after=True,status_after=True):
        count_step=0
        while count_step < n:    
            if self.toElim == 0:
                self.STATUS = STATUS_AUX_ELIMINATED
                break
            #mem0=rss()
            self.set_rho()
            self.findMRC()
            #mem1=rss()
            #self.mem[0] += mem1-mem0
            #
            self.costs += [self.curCost()]
            self.redcosts += [self.MRC]
            
            if self.MRC >= - pf("1.0e-60"):
                self.STATUS = STATUS_COST_MINIMIZED
                break
            # find pivot
            self.AIb = self.timesAbInv(self.b)
            
            #mem2=rss()
            #self.mem[1] += mem2-mem1
            #
            self.Aa = self.lp_table.vecfunlist[self.l_MRC].value(self.d_MRC)
            
            #mem3=rss()
            #self.mem[2] += mem3-mem2
            #
            self.AIa = self.timesAbInv(self.Aa)
            
            #mem4=rss()
            #self.mem[3] += mem4-mem3
            #
            zero = mpfr_array.zeros((1,));
            inf = mpfr_array.array(['inf']);
            
            #mem5=rss()
            #self.mem[4] += mem5-mem4
            #
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
            self.Xb[self.pivot] = ('OPE', self.l_MRC, self.d_MRC)
            self.Cb[self.pivot] = "0"
            self.set_inverse()
            
            #mem6=rss()
            #self.mem[5] += mem6-mem5
            #
            self.iter += 1
            count_step += 1
                
        if sort_after:
            self.sort_basis(False)
        if status_after:
            self.status()
            
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
            print pairlist

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
        print "Current cost:", self.curCost()
        print "AUX variables:", self.toElim
        print "Xb:", self.Xb
        print "Solution:", to_double(self.curSol())
        print "Iterations:", self.iter    
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



