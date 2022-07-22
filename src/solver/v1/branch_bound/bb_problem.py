import numpy as np
from scipy import interpolate
import operator
import sys
import Taylor_methods
reload(Taylor_methods)
from Taylor_methods import cubic,candList

import utils.logging as log

COMP_TOLERANCE = 1.e-12 # comparison tolerance

np.seterr(all='raise')

################## BBinterval class definition
class BB_Interval:
    """will use as a structure with fields
    x0,x1,i0,i1,Lb,Ub,min_attained,x_min"""
    def __init__(self,x0,x1):
        self.x0 = x0 
        self.x1 = x1 # beginning and end of the interval
        self.i0= None
        self.i1= None # (xs[i0],xs[i1]) is the interval covering (x0,x1)
        self.Lb = None 
        self.Ub = None # lower and upper bound on the minimum
        self.x_min = None # the point where Lb is attained
    
    def __repr__(self): # used when printing
        return ("<(%s,%s) Lb=%s,Ub=%s,x_min=%s,i0=%s,i1=%s>" %
        (self.x0,self.x1,self.Lb,self.Ub,self.x_min,self.i0,self.i1))


###################### BB class definition

class BB_Problem:
    """ Branch and bound function, relevant for our applications. 
    it takes as arguments
    1) cubic interpolated function
    2) endpoints of the intreval on which minimization has to be performed
    3) method keyword
    """    
#----------------------------------------------------------
    def __init__(self, f_int, x0, x1, method = "local"):

        if (f_int[2]!=3):
            sys.exit("BBproblem.__init__:interpolation order must be 3")
        self.f_int = f_int
        self.x0, self.x1 = x0,x1
 
        xs = np.unique(f_int[0])  #make a deduplicated list of knots

        if (xs[0]>self.x0 + COMP_TOLERANCE or xs[-1]<self.x1 - COMP_TOLERANCE):
            errstr='BBprob: (%f, %f) vs (%f,%f)'%(xs[0], xs[-1],
                    self.x0+COMP_TOLERANCE, self.x1-COMP_TOLERANCE)
            print errstr
            log.critical(errstr)
            #self.x0 = max(xs[0], self.x0 + COMP_TOLERANCE)
            #self.x1 = min(xs[-1], self.x1 - COMP_TOLERANCE)
            #log.critical('BBprob: manually setting requested range within bounds: (%f,%f)'%
            #        (self.x0, self.x1))
            #sys.exit("BBproblem.__init__:minimization interval\
            # extends out of the range of the interpolation function")
            raise ValueError(errstr)
        # determine the smallest set of interpolation nodes covering (x0,x1):
        
        self.xlist = np.array(xs[max(np.searchsorted(xs,self.x0) - 1,0):
                        np.searchsorted(xs,self.x1) + 1])
        # shift x-  values to the middle
        # of the intervals to avoid ambiguities when computing 3rd derivative :
    
        self.der3 = interpolate.splev((self.xlist[0:-1] + self.xlist[1:])/2.,
                                    self.f_int, 3)
               
        self.config_minmax_est(method)
        
        I0 = BB_Interval(x0,x1) # first interval of BB hierarchy
        I0.i0, I0.i1 = 0,len(self.xlist)-1 # the range of xlist covering I0
        self.minmax_est(I0) # set minmax_est variables
        
        self.ilist = [ I0 ]
        self.iL = 0
        self.iter = 0  
#----------------------------------------------------------
    def config_minmax_est (self,method):
        self.METHOD = method
        if method == "global":
            self.minder3 = min( self.der3 )
            self.maxder3 = max( self.der3 )
            self.minmax_est = self.minmax_est_global
        elif method == "semilocal":
            self.minder3 = min( self.der3 )
            self.maxder3 = max( self.der3 )
            self.minmax_est = self.minmax_est_semilocal
        elif method == "local":
            self.minmax_est = self.minmax_est_local
        else:
            sys.exit("config_minmax_est called for unknown method")
#----------------------------------------------------------
    def bb_step(self):
        I = self.ilist[self.iL] # interval to subdivide
        t0, t1 = I.x0, I.x1
        Ileft = BB_Interval(t0, (t0+t1)/2.)
        Iright = BB_Interval((t0+t1)/2., t1) # initialize the subintervals
        if self.METHOD in ["local", "semilocal"]: 
            deltai = np.searchsorted(self.xlist[I.i0:I.i1+1],(t0+t1)/2.)
            Ileft.i0, Ileft.i1 = I.i0, I.i0+deltai
            Iright.i0, Iright.i1 = I.i0+deltai-1, I.i1 
        # for each half, find the covering union of interpolating subintervals
        self.minmax_est(Ileft)
        self.minmax_est(Iright) #set the minmax variables
        self.ilist = self.ilist[:self.iL] + [ Ileft,Iright ] +\
                                                        self.ilist[self.iL+1:]
        self.set_bb_vars() 
        self.iter += 1
#----------------------------------------------------------       
    def set_bb_vars(self):
        """sets the values of the internal variables:
        Lb,Ub the minimal lower and upper bound; 
        iL, iU - indices of the intervals on which these bounds are realized"""
        self.iL = min( enumerate([x.Lb for x in self.ilist]), key=lambda x: x[1])[0]
#----------------------------------------------------------           
    def findmin(self, diff_goal=1.e-14  ,steps = 10000, threshold = float('inf')):
        for iter in range(steps): 
            I = self.ilist[self.iL]
            if I.Ub - I.Lb > diff_goal and I.Lb < threshold: 
            # only minima below threshold are interesting
                self.bb_step()
            else:
                break
#----------------------------------------------------------           
    def findmin_relative_negative(self, relative_diff_goal=0.01, steps = 10000, 
                    threshold=0., greed_factor=2.):
        """find an estimate for the minimum terminating the search if
        - cannot improve at least greed_factor over current threshold
        assumed <=0 (in particular if Lb is positive)
        - vertical relative precision relative_diff_goal is reached
        """
        if not threshold <= 0.:
            errstr="threshold must be <=0:  "+str(threshold)+", FAILING..."
            log.critical(errstr)
            raise ValueError(errstr)
            #sys.exit("threshold must be <=0:  "+str(threshold))
        if not greed_factor > 1.:
            errstr="greed_factor must be >1 (currently %f), failing"%greed_factor
            log.critical(errstr)
            raise ValueError(errstr)
            #sys.exit("greed_factor must be >1")
        if not 0. < relative_diff_goal < 1. - 1./greed_factor:
            errstr="Must have 0. < relative_diff_goal < 1. - 1./greed_factor\
                    (currently %f)"%relative_diff_goal
            log.critical(errstr)
            raise ValueError(errstr)
            #sys.exit("Must have 0. < relative_diff_goal < 1. - 1./greed_factor")
            
        for iter in range(steps): 
            I = self.ilist[self.iL]
            if I.Lb < threshold*greed_factor and\
                I.Ub > (1. - relative_diff_goal)*I.Lb: 
                # i.e. 1) only minima significantly below current best are interesting
                # and 2) terminate when [Lb,Ub] becomes short compared to its bottom
                self.bb_step()
            else:
                break
#----------------------------------------------------------    
    def status(self,print_ilist=False):
        print "=======BB problem status"
        print "x0, x1=", self.x0, self.x1
        print "METHOD:", self.METHOD
        print len(self.ilist), "intervals"
        print "iL, Lb, Ub =", self.iL, self.ilist[self.iL].Lb, self.ilist[self.iL].Ub
        print "iter=", self.iter
        print "accuracy:", self.ilist[self.iL].Ub - self.ilist[self.iL].Lb
        if print_ilist:
            print "interval list:" 
            print self.ilist
############# Three minmax_est functions defined here  
    def minmax_est_global (self, I):
            """estimate min & max using global min and max of der3"""
            a = [interpolate.splev(I.x0 , self.f_int, k) for k in range(3)]
            a[2]=a[2]/2.
            I.Lb, I.Ub, I.x_min =\
            Taylor_minmax_est(a, [self.minder3/6., self.maxder3/6.], I.x0, I.x1)
#----------------------------------------------------------          
    def minmax_est_semilocal ( self, I ):
            """estimate min & max using global min and max of der3,
            except if the interval is entirely inside one interpolation
            subinterval"""
            if (self.xlist[I.i0] > I.x0 + COMP_TOLERANCE 
            or self.xlist[I.i1] < I.x1 - COMP_TOLERANCE):
                errstr="minmax_est_semilocal: error in BB interval\
                        specifications: xlist entries (%f,%f) vs x0,x1(%f,%f)"\
                        %(self.xlist[I.i0], self.xlist[I.i1], I.x0 +
                                COMP_TOLERANCE, I.x1 - COMP_TOLERANCE)
                log.critical(errstr)
                raise ValueError(errstr)
                #sys.exit("minmax_est_semilocal: error in BB interval\
                # specifications")
            a = [interpolate.splev(I.x0 , self.f_int, der=k) for k in range(3)]
            a[2]=a[2]/2.
            if I.i1 == I.i0+1: #search interval inside one interpolation subinterval
                minder3 = maxder3 = self.der3[I.i0]
            else:
                minder3 = self.minder3
                maxder3 = self.maxder3
            I.Lb, I.Ub, I.x_min =\
            Taylor_minmax_est(a, [minder3/6., maxder3/6.], I.x0, I.x1)
#----------------------------------------------------------  
    def minmax_est_local ( self, I ):
            """estimate min & max using min and max of der3
            over the union of interpolation subintervals covering I"""
            if (self.xlist[I.i0] > I.x0 + COMP_TOLERANCE 
                or self.xlist[I.i1] < I.x1 - COMP_TOLERANCE):
                errstr="minmax_est_local: error in BB interval\
                        specifications: xlist entries (%f,%f) vs x0,x1(%f,%f)"\
                        %(self.xlist[I.i0], self.xlist[I.i1], I.x0 +
                                COMP_TOLERANCE, I.x1 - COMP_TOLERANCE)
                log.critical(errstr)
                raise ValueError(errstr)
                #sys.exit("minmax_est_local: error in BB interval specifications")
            a = [interpolate.splev((I.x0 + I.x1)/2., self.f_int, k) for k in range(3)]
            a[2]=a[2]/2.
            
            der3need = self.der3[I.i0:I.i1]
            I.minder3 = np.amin( der3need )
            I.maxder3 = np.amax( der3need )
            #I.minder3 = min( der3need )
            #I.maxder3 = max( der3need )
            
            I.Lb, I.Ub, I.x_min =\
            Taylor_minmax_est(a, [I.minder3/6., I.maxder3/6.], I.x0, I.x1)
            
############# End of minmax_est function definitions            


def Taylor_minmax_est(a, minmax_a3, x0, x1):
    """Returns lower and upper bounds for min and max of function f on [x0,x1],
    given 
    1) the array of its TAYLOR COEFFS a = [a0,a1,a2] at x=xm=(x0+x1)/2 and 
    2) lower ad upper bound for the third TAYLOR COEFF a3: 
    minmax_a3=[mina3,maxa3]"""
    a0, a1, a2 = a
    min3, max3 = minmax_a3 
    xm = (x0 + x1)/2.
    cand_list_max = candList(a1, a2, max3, x0, x1)
    cand_list_min = candList(a1, a2, min3, x0, x1)
    # _cubic(a,max3,x) is an upper bound for x>=xm, lower bound x<xm,
    # vice versa for _cubic(a,min3,x) 
    x_min_est, Ub = min([[x, cubic(a0, a1, a2, max3, x - xm)] for x in cand_list_max if x >= xm]+
                         [[x, cubic(a0, a1, a2, min3, x - xm)] for x in cand_list_min if x < xm], 
                         key = lambda v: v[1])
    Lb           = min([cubic(a0, a1, a2, max3, x - xm) for x in cand_list_max if x < xm]+
                         [cubic(a0, a1, a2, min3, x - xm) for x in cand_list_min if x >= xm])
                         
 
    # NB:1) x_min_est must be computed as the point where Ub (not Lb!) estimate is minimal
    # 2) If Ub=Lb then x_min_est is the real minimum, this can happen if min3=max3 
    #(so that full function is a cubic polynomial) but also if max3-min3 is small and x_min_est=x0
    return [ Lb, Ub, x_min_est]    
 
############################  end BB class

