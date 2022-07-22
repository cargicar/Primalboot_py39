#cython: profile=True

###
# S. Rychkov June-July 2013
# Contains two classes: CB which encodes vector of conformal block components
# (contracted with sigma-tensor) for a given spin
# CB_Component - encodes a single CB component for a given spin, or a CB contracted with
# rho-vector

from libc.stdlib cimport malloc, free
from libc.stdio cimport stdout, sprintf, fprintf
import io # Carlos add

from c_mpfr cimport *
from c_mpfi cimport * 
cimport c_polyval 

import solver.v2.mpfr_array.mpfr_array as mpfr_array 
cimport solver.v2.mpfr_array.mpfr_array as mpfr_array 

import solver.v2.prec_float.prec_float as PF
cimport solver.v2.prec_float.prec_float as PF

import numbers
import numpy as np
#from matplotlib.pyplot import *

#import cbdata as cbdata

import solver.v2.rho_repr.cbdata as cbdata # Carlos edit 


def pf(x):
    return PF.prec_float(x,prec=212)


#
## this function multiples CB by rho-vector and returns CB_Component

# Carlos comment: In the future the rho in rhomult should be renamed (noit to get confused with the radial cross ratio)
def rhomult( rho, cb):
    return CB_Component( cb.prec, cb.delta0, mpfr_array.dot(<mpfr_array.ndarray>rho, <mpfr_array.ndarray>cb.Ps),
                        cb.poles,
                        mpfr_array.dot(<mpfr_array.ndarray>rho, <mpfr_array.ndarray>cb.polecoeffs))

## this function multiples CB by rho-vector and updates the existing CB_Component
# (for reducing memory leaks)
def rhomult_update( rho, cb, cbcomp):
    (<mpfr_array.ndarray>rho)._vm_dot_assign(<mpfr_array.ndarray>cb.Ps,
                                                                   <mpfr_array.ndarray>cbcomp.P)
    (<mpfr_array.ndarray>rho)._vm_dot_assign(<mpfr_array.ndarray>cb.polecoeffs,
                                                                 <mpfr_array.ndarray>cbcomp.polecoeffs)
    return cbcomp
 

###########################################################
# each conformal block component has the form
# rho^x (p(x-delta0) +q(x)). x= operator dimension, rho = 3-sqrt(8)
# p(x) is a polynomial of degree 2 nmax + mmax
# q(x) is a rational function of the form sum c(i)/(x-x(i))
#(simple poles only - this condition won't work in d=2,4)

cdef class CB:
    cdef public mpfr_array.ndarray Ps,poles,polecoeffs,maxlist
    cdef public list pstring  # Carlos edit
    cdef public int prec
    cdef public int CBlen
    cdef public PF.prec_float eps
    cdef public PF.prec_float rho
    cdef public PF.prec_float delta0
    
    def __init__(self, source):
        #if isinstance(source, 'file'):
        if isinstance(source,io.IOBase): # Carlos edit
            self.from_file(source)
        elif isinstance(source, cbdata.CBLdata): # Carlos edit: Comment out this two lines to run test_solve.py. 
            self.from_cbldata(source)
    
    # CBLdata should _already_ have its sigma set!
    def from_cbldata(self, cbldata):
        self.prec = cbldata.prec
        self.CBlen = cbldata.CBlen
        
        #point around which the polynomials are expanded
        #chosen at the unitarity bound
        self.delta0 = cbldata.delta0
        
        ncoeffsP = cbldata.ncoeffsP
        self.Ps = cbldata.sigPs
        #read coefficients of all polynomials p(x-delta0) (one polynomial per each of CBlen components)
        #and store in a two-dimensional mpfr_array
        #npoles = cbldata.npoles
        self.poles = cbldata.goodpoles
        # read the poles x(i) and store in mpfr_array. Poles are the same for each component
        
        self.polecoeffs = cbldata.goodpolecoeffs
        # read the coefficients of the poles c(i) (different for each of CBlen components)
        #and store in mpfr_array
        
        self.rho = PF.prec_float(3,self.prec)- PF.sqrt(PF.prec_float(8,self.prec))
        #initialize rho=3-sqrt(8) for further use

    def from_file(self, FILE):
        #initialize from file
        self.prec = int(FILE.readline()) #precision
        self.CBlen = int(FILE.readline()) #vector length
        
        self.delta0 = PF.prec_float(FILE.readline(),self.prec)
        #point around which the polynomials are expanded
        #chosen at the unitarity bound
        
        ncoeffsP = int(FILE.readline()) #number of nonzero polynomial coefficients, 2nmax+mmax
        
        Pstring = [ [FILE.readline() for i in range(ncoeffsP)] for n in range(self.CBlen)]
        self.pstring=Pstring                # Carlos edit
        self.Ps = mpfr_array.array(Pstring)
        
        #read coefficients of all polynomials p(x-delta0) (one polynomial per each of CBlen components)
        #and store in a two-dimensional mpfr_array
        
        
        npoles = int(FILE.readline()) # number of poles x(i)
        
        poles_string = [ FILE.readline() for i in range(npoles)] 
        self.poles = mpfr_array.array(poles_string)
        # read the poles x(i) and store in mpfr_array. Poles are the same for each component
        
        polecoeffstring = [ [FILE.readline() for i in range(npoles)] for n in range(self.CBlen)]
        self.polecoeffs = mpfr_array.array(polecoeffstring)
        # read the coefficients of the poles c(i) (different for each of CBlen components)
        #and store in mpfr_array
        
        self.rho = PF.prec_float(3,self.prec)- PF.sqrt(PF.prec_float(8,self.prec))
        #initialize rho=3-sqrt(8) for further use
        
    def value(self, PF.prec_float opdim):
        #this function computes the value of CB for a given delta;
        #it's used when updating linear programming basis matrix at the end of each step
        cdef PF.prec_float pvalue, polevalue, fac
        cdef int n 
        
        ret_list = [ "0" for n in range(self.CBlen)]
        pvalue = PF.prec_float(0, prec = self.prec)
        polevalue = PF.prec_float(0, prec = self.prec)
        fac = self.rho**opdim
        for n in range(self.CBlen):
            c_polyval.p0(self.delta0.data, (<mpfr_array.ndarray>self.Ps[n]).data, (<mpfr_array.ndarray>self.Ps[n]).size, opdim.data, pvalue.data)
            c_polyval.r0(self.poles.data, (<mpfr_array.ndarray>self.polecoeffs[n]).data, self.poles.size, opdim.data, polevalue.data)
            ret_list[n] = fac * (pvalue+polevalue)
        
        return mpfr_array.array(ret_list)
     
###########################################################
cdef class CB_Component: #inidividual CB component or CB contracted with a rho-vector
    cdef public mpfr_array.ndarray P,poles,polecoeffs
    cdef public int prec
    cdef public PF.prec_float delta0, rho, logrho
    cdef PF.prec_interval rhoint, logrhoint
    
    def __init__(self, prec, delta0, Parray, polearray, polecoeffarray):
        self.prec = prec
        self.delta0 = delta0
        self.P = Parray 
        self.poles = polearray
        self.polecoeffs = polecoeffarray
        #these variables have the same meaning as for CB class
        
        self.rho = PF.prec_float(3,self.prec) - PF.sqrt(PF.prec_float(8,self.prec))
        self.logrho = self.rho.log()
        self.rhoint= PF.prec_interval(self.rho, self.rho, prec=self.prec)
        self.logrhoint= PF.prec_interval(self.logrho, self.logrho, prec=self.prec)
        #just initialize rho, logrho once and for all for future uses
            
    cpdef PF.prec_float value(self, PF.prec_float opdim):
        # computes the value of the component for a given Delta=opdim
        cdef PF.prec_float pvalue, polevalue
        pvalue = PF.prec_float(0, prec = self.prec)
        polevalue = PF.prec_float(0, prec = self.prec)
        c_polyval.p0(self.delta0.data, self.P.data, self.P.size, opdim.data, pvalue.data)
        c_polyval.r0(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, polevalue.data)
        return (self.rho**opdim) * (pvalue + polevalue)
    
    cpdef PF.prec_float valueder(self, PF.prec_float opdim):
        # computes the value of the derivative for a given Delta=opdim
        cdef PF.prec_float p0,p1,r0,r1
        p0 = PF.prec_float(0, prec = self.prec)
        r0 = PF.prec_float(0, prec = self.prec)
        p1 = PF.prec_float(0, prec = self.prec)
        r1 = PF.prec_float(0, prec = self.prec)
        c_polyval.p0(self.delta0.data, self.P.data, self.P.size, opdim.data, p0.data)
        c_polyval.r0(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r0.data)
        c_polyval.p1(self.delta0.data, self.P.data, self.P.size, opdim.data, p1.data)
        c_polyval.r1(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r1.data)
        return (self.rho**opdim) * (p1 + r1 + (p0 + r0) * self.logrho)
    
    cpdef PF.prec_float valueder2(self, PF.prec_float opdim):
        # computes the value of the second derivative for a given Delta=opdim
        cdef PF.prec_float p0,p1,r0,r1,p2,r2
        p0 = PF.prec_float(0, prec = self.prec)
        r0 = PF.prec_float(0, prec = self.prec)
        p1 = PF.prec_float(0, prec = self.prec)
        r1 = PF.prec_float(0, prec = self.prec)
        p2 = PF.prec_float(0, prec = self.prec)
        r2 = PF.prec_float(0, prec = self.prec)
        c_polyval.p0(self.delta0.data, self.P.data, self.P.size, opdim.data, p0.data)
        c_polyval.r0(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r0.data)
        c_polyval.p1(self.delta0.data, self.P.data, self.P.size, opdim.data, p1.data)
        c_polyval.r1(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r1.data)
        c_polyval.p2(self.delta0.data, self.P.data, self.P.size, opdim.data, p2.data)
        c_polyval.r2(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r2.data)
        return (self.rho**opdim) * (p2 + r2
                                    + PF.prec_float(2, prec = self.prec) * self.logrho * (p1 + r1)
                                    + (p0 + r0) * self.logrho * self.logrho)
    
    cpdef PF.prec_float valueder3(self, PF.prec_float opdim):
        # computes the value of the third derivative for a given Delta=opdim
        cdef PF.prec_float p0,p1,r0,r1,p2,r2,p3,r3
        p0 = PF.prec_float(0, prec = self.prec)
        r0 = PF.prec_float(0, prec = self.prec)
        p1 = PF.prec_float(0, prec = self.prec)
        r1 = PF.prec_float(0, prec = self.prec)
        p2 = PF.prec_float(0, prec = self.prec)
        r2 = PF.prec_float(0, prec = self.prec)
        p3 = PF.prec_float(0, prec = self.prec)
        r3 = PF.prec_float(0, prec = self.prec)
        c_polyval.p0(self.delta0.data, self.P.data, self.P.size, opdim.data, p0.data)
        c_polyval.r0(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r0.data)
        c_polyval.p1(self.delta0.data, self.P.data, self.P.size, opdim.data, p1.data)
        c_polyval.r1(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r1.data)
        c_polyval.p2(self.delta0.data, self.P.data, self.P.size, opdim.data, p2.data)
        c_polyval.r2(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r2.data)
        c_polyval.p3(self.delta0.data, self.P.data, self.P.size, opdim.data, p3.data)
        c_polyval.r3(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r3.data)
        return (self.rho**opdim) * (p3 + r3
                                    + PF.prec_float(3, prec = self.prec) * self.logrho * (p2 + r2)
                                    + PF.prec_float(3, prec = self.prec) * self.logrho * self.logrho * (p1 + r1)
                                    + (p0 + r0) * self.logrho * self.logrho * self.logrho)
    
    def valueder123(self, PF.prec_float opdim):
        # computes the value of the 1st, 2nd, 3rd derivative in one go
        cdef PF.prec_float p0,p1,r0,r1,p2,r2,p3,r3
        p0 = PF.prec_float(0, prec = self.prec)
        r0 = PF.prec_float(0, prec = self.prec)
        p1 = PF.prec_float(0, prec = self.prec)
        r1 = PF.prec_float(0, prec = self.prec)
        p2 = PF.prec_float(0, prec = self.prec)
        r2 = PF.prec_float(0, prec = self.prec)
        p3 = PF.prec_float(0, prec = self.prec)
        r3 = PF.prec_float(0, prec = self.prec)
        c_polyval.p0(self.delta0.data, self.P.data, self.P.size, opdim.data, p0.data)
        c_polyval.r0123(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, r0.data, r1.data, r2.data,r3.data)
        c_polyval.p1(self.delta0.data, self.P.data, self.P.size, opdim.data, p1.data)
        c_polyval.p2(self.delta0.data, self.P.data, self.P.size, opdim.data, p2.data)
        c_polyval.p3(self.delta0.data, self.P.data, self.P.size, opdim.data, p3.data)
        return (
            (self.rho**opdim) * (p1 + r1 + (p0 + r0) * self.logrho),
            (self.rho**opdim) * (p2 + r2
                                    + PF.prec_float(2, prec = self.prec) * self.logrho * (p1 + r1)
                                    + (p0 + r0) * self.logrho * self.logrho),
            (self.rho**opdim) * (p3 + r3
                                    + PF.prec_float(3, prec = self.prec) * self.logrho * (p2 + r2)
                                    + PF.prec_float(3, prec = self.prec) * self.logrho * self.logrho * (p1 + r1)
                                    + (p0 + r0) * self.logrho * self.logrho * self.logrho)
            )
    
    def valueder123fast(self, PF.prec_float opdim):
        # computes the value of the 1st, 2nd, 3rd derivative in one go
        #using compiled C function
        cdef PF.prec_float d1,d2,d3
        
        d1 = PF.prec_float(0, prec = self.prec)
        d2 = PF.prec_float(0, prec = self.prec)
        d3 = PF.prec_float(0, prec = self.prec)
        
        c_polyval.valueder123(self.rho.data, self.logrho.data,
                               self.delta0.data, self.P.data, self.P.size,
                               self.poles.data, self.polecoeffs.data, self.poles.size,
                               opdim.data,
                               d1.data, d2.data, d3.data
                               )
        return (d1,d2,d3)
        
    def valueder123fast1(self, PF.prec_float opdim, PF.prec_float d1, PF.prec_float d2, PF.prec_float d3):
        # computes the value of the 1st, 2nd, 3rd derivative in one go
        #using compiled C function
        
        c_polyval.valueder123(self.rho.data, self.logrho.data,
                               self.delta0.data, self.P.data, self.P.size,
                               self.poles.data, self.polecoeffs.data, self.poles.size,
                               opdim.data,
                               d1.data, d2.data, d3.data
                               )
    
    cpdef PF.prec_interval bounds(self, PF.prec_interval opdim):
        # computes the bounds of the component for a given Delta=opint
        # NOT CURRENTLY USED
        cdef PF.prec_interval pvalue, polevalue
        pvalue = PF.prec_interval(0, 0, prec = self.prec)
        polevalue = PF.prec_interval(0, 0, prec = self.prec)
        c_polyval.i_p0(self.delta0.data, self.P.data, self.P.size, opdim.data, pvalue.data)
        c_polyval.i_r0(self.poles.data, self.polecoeffs.data, self.poles.size, opdim.data, polevalue.data)
        if pvalue.length() > polevalue.length():
            self.cnt[0]+=1 #check which error dominates; one sees that it's always the first one
        else:
            self.cnt[1]+=1
        return (self.rhoint**opdim) * (pvalue + polevalue)
    
    cpdef PF.prec_interval boundsder(self, PF.prec_interval opdim):
        # computes the bounds of the component derivative for a given Delta=opint
        # NOT CURRENTLY USED
        cdef PF.prec_interval p0,p1, q0,q1
        p0 = PF.prec_interval(0, 0, prec = self.prec)  
        p1 = PF.prec_interval(0, 0, prec = self.prec)
        q0 = PF.prec_interval(0, 0, prec = self.prec)
        q1 = PF.prec_interval(0, 0, prec = self.prec)
        c_polyval.i_p01(self.delta0.data, self.P.data, self.P.size, opdim.data, p0.data, p1.data)
        #c_polyval.i_p0(self.delta0.data, self.P.data, self.P.size, opdim.data, p0.data)
        #c_polyval.i_p1(self.delta0.data, self.P.data, self.P.size, opdim.data, p1.data)
        c_polyval.i_r01(self.poles.data, self.polecoeffs.data,
                                  self.poles.size, opdim.data, q0.data , q1.data)
        
        #c_polyval.i_r1together(self.poles.data, self.polecoeffs.data,
                                  #self.poles.size, opdim.data, q0.data, q1.data)
        return (self.rhoint**opdim) * (p1 + q1 + (p0 + q0) * self.logrhoint)

    def exam(self):
        # output parameters of the component for further printing; for debugging purposes 
        return [self.delta0,
                mpfr_array.to_double(self.P),
         mpfr_array.to_double(self.poles),
         mpfr_array.to_double(self.polecoeffs)]
    
#    def plot(self, x0, x1, plotpoints = 100):
#        # plot the component, for debugging purposes
#        if isinstance(x0,str):
#            x0pf = PF.prec_float(x0,prec=212)
#            x1pf = PF.prec_float(x1,prec=212)
#        x = [x0pf + pf(k)*(x1pf-x0pf)/pf(plotpoints)   for k in range(plotpoints)]
#        y = [self.value(xi) for xi in x]
#        xnumpy = np.array([PF.to_double(xi) for xi in x])
#        ynumpy = np.array([PF.to_double(xi) for xi in y])
#        plot( xnumpy,ynumpy)
#        show()
        
                


