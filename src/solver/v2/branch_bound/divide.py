import operator
import sys
import solver.v2.prec_float.prec_float as PF
import numpy
#from matplotlib.pyplot import *
import resource
from utils.memory import rss

def pf(x):
    return PF.prec_float(x,prec=212)

COMP_TOLERANCE = pf(1.e-12) # comparison tolerance

number_BB_Interval = 0
################## BBinterval class definition
class BB_Interval:
    def __init__(self,l,x0,x1):
        global number_BB_Interval
        self.l = l #spin
        self.x0 = x0 
        self.x1 = x1 # beginning and end of the interval
        self.type = "unknown"
        self.f0L = None
        self.f1L = None
        self.f2L = None
        self.f0R = None
        self.f1R = None
        self.f2R = None
        self.f1c = None
        self.f2c = None
        self.xmin = None
        self.fmin = None
        self.leftborder = False
        self.rightborder = False
        number_BB_Interval +=1
    def __repr__(self): # used when printing
        return ("<l=%s (%s,%s) %s %s %s %s %s>" %
        (self.l, self.x0,self.x1,self.type, self.leftborder, self.rightborder, self.xmin, self.fmin))
    
    def __del__(self):
        global number_BB_Interval
        number_BB_Interval -=1
###################### BB class definition


class Divide:
    """ takes at initialization:
    1) list of function in rho_representation
    2) list of intervals on which the functions are defined"""
#----------------------------------------------------------
    def __init__(self, funlist, spectrum):

        self.funlist = funlist
        self.spectrum = spectrum
        
        self.ilist = [ BB_Interval(specint.l,specint.d0,specint.d1) for specint in self.spectrum]
        # print("ilist %s \n" %self.ilist )
        
        for x in self.ilist:
            x.leftborder = True
            x.rightborder = True
            
        self.iter = 0
        self.curint =0
        
    #----------------------------------------------------------
    def reset(self):
        
        for specint in self.ilist:
            del specint
            
        self.ilist = [ BB_Interval(specint.l,specint.d0,specint.d1) for specint in self.spectrum]
        
        for x in self.ilist:
            x.leftborder = True
            x.rightborder = True
            
        self.iter = 0
        self.curint =0
    
    def good_interval(self,n):
        I = self.ilist[n] # interval to check  n=interval label
                                    #initially each interval is label by spin, but 
            # later on interval will be subdivided.
                            
        
        if I.f1L == None: # set left 0,1derivatives
            I.f1L = self.funlist[I.l].valueder(I.x0)
        
        if I.f1R == None: # set right 0,1 derivatives
            I.f1R = self.funlist[I.l].valueder(I.x1)
            
        I.c = (I.x0+I.x1)/pf(2)
        dx = (I.x1 - I.x0)/pf(2)
        #(I.f1c, I.f2c, I.f3c) = self.funlist[I.l].valueder123fast(I.c)
        I.f1c = pf(0)
        I.f2c = pf(0)
        I.f3c = pf(0)
        self.funlist[I.l].valueder123fast1(I.c,I.f1c,I.f2c,I.f3c)
   
        if (
            (I.f1c + dx * I.f2c + dx*dx/pf(2) * I.f3c - I.f1R).abs() < pf("0.1") * I.f1R.abs()
            and 
            (I.f1c - dx * I.f2c + dx*dx/pf(2) * I.f3c - I.f1L).abs() < pf("0.1") * I.f1L.abs()
           ):
             #the above two conditions suggest that the first derivative is well approximated
                # by a quadratic polynomial P(x)
                # If f1L and f1R have different sign, we conclude it's a good interval.
                #
                # in the case that f1L and f1R have the same sign (say positive),
                # we want to check additionally that the minimum of P(x) is
                # larger than 0.5 times the min(f1L,f1R)
                # This should be a safeguard that the function f1 does not cross the zero
                # within the interval.
                # This extra check need to be performed only if f3c has the same sign as f1L and f1R
                # (otherwise is automatic)
            if I.f1L > pf(0) and I.f1R > pf(0) and I.f3c < pf(0):
                dxmin = -I.f2c/I.f3c # where P(x) has minimum
                if dxmin.abs() < dx: #minimum within the interval
                    Pmin = I.f1c + dxmin * I.f2c/pf(2)
                    if Pmin > pf("0.5") * min(I.f1L,I.f1R):
                        return True # good interval
                    else:
                        return False # continue divisions 
            if I.f1L < pf(0) and I.f1R < pf(0) and I.f3c > pf(0):
                dxmin = -I.f2c/I.f3c # where P(x) has minimum
                if dxmin.abs() < dx:
                    Pmin = I.f1c + dxmin * I.f2c/pf(2)
                    if Pmin < pf("0.5") * max(I.f1L,I.f1R):
                        return True # good interval
                    else:
                        return False # continue divisions
            return True
        else:
            return False
        
    def good_interval1(self,n):
        I = self.ilist[n] # interval to check
        
        if I.f1L == None: # set left 0,1derivatives
            I.f1L = self.funlist[I.l].valueder(I.x0)
        
        if I.f1R == None: # set right 0,1 derivatives
            I.f1R = self.funlist[I.l].valueder(I.x1)
            
        I.c = (I.x0+I.x1)/pf(2)
        I.f1c = self.funlist[I.l].valueder(I.c)
        I.f2c = self.funlist[I.l].valueder2(I.c)
        # Carlos commment: Notice that this function compares only up to linear order
        if (
            (I.f1c + (I.x1 - I.x0)/pf(2) * I.f2c - I.f1R).abs() < pf("0.1") * I.f1R.abs()
            and 
            (I.f1c - (I.x1 - I.x0)/pf(2) * I.f2c - I.f1L).abs() < pf("0.1") * I.f1L.abs()
           ):
            return True
        else:
            return False
            
    def declare_interval(self,n):
        I = self.ilist[n]
        if I.f1L< pf(0):
            if I.f1R < pf(0):
                I.type = "decr"
            else:
                I.type = "min"
        else:
            if I.f1R< pf(0):
                I.type = "max"
            else:
                I.type = "incr"
            
#----------------------------------------------------------
    def divide_step(self): 
        if self.good_interval(self.curint):
            self.declare_interval(self.curint)
            self.curint += 1 
        else:
            I = self.ilist[self.curint]
            Ileft = BB_Interval(I.l, I.x0, (I.x0+I.x1)/pf(2))
            Iright = BB_Interval(I.l, (I.x0+I.x1)/pf(2), I.x1) # initialize the subintervals
            Ileft.f1L = I.f1L
            Ileft.f1R = I.f1c
            Ileft.leftborder = I.leftborder
            
            Iright.f1R = I.f1R
            Iright.f1L = I.f1c
            Iright.rightborder = I.rightborder
            
            self.ilist = self.ilist[:self.curint] + [ Ileft,Iright ] +\
                                                        self.ilist[self.curint+1:]

#----------------------------------------------------------           
    def findmin(self):
        
        while self.curint < len(self.ilist):
            # This loop perform all the subdivitions and label the final intervals
            self.divide_step()
            
        # at this point we have to select the intervals which can have minimum
         # left-rightborder bolean variable declares that the border of the given interval is 
        # the min/max value of dimension for a given spin (i.e d_min, d_max)
        
        self.ilist1 = [x for x in self.ilist if (x.type == "min"
                                                # or ( x.type == "decr" and x.rightborder == True) # Carlos edit
                                                 or ( x.type == "incr" and x.leftborder == True) )]

        
        if self.ilist1 == []: # this can occur in parallel version 
            I = BB_Interval(0,pf(0),pf(0)) # create fictitious interval
            I.fmin = pf(1)
            return I
        

        for I in self.ilist1:
            if I.type == "min":
                self.Newton(I)
            else: # x.type == "incr"
                I.xmin = I.x0
                I.fmin = self.funlist[I.l].value(I.x0)
                
        return min(self.ilist1, key = lambda x: x.fmin)
    
#----------------------------------------------------------           
    def find_negative_local_minima(self):
        while self.curint < len(self.ilist):
            self.divide_step()
            
        # at this point we have to select the intervals which can have minimum
        
        self.ilist1 = [x for x in self.ilist if (x.type == "min"
                                                 # or ( x.type == "decr" and x.rightborder == True)
                                                 or ( x.type == "incr" and x.leftborder == True) )]
        self.ilist2=[]
        
        for I in self.ilist1:
            if I.type == "min":
                self.Newton(I)
                if I.fmin < pf("1.0e-10"):
                    self.ilist2 += [I]
            elif I.type == "incr":
                I.xmin = I.x0
                I.fmin = self.funlist[I.l].value(I.x0)
                if I.fmin < pf("1.0e-10"):
                    self.ilist2 += [I]
            else:
                #x.type="decr"
                pass
                #I.xmin = I.x1
                #I.fmin = pf("0")
                #self.ilist2 += [I]
                #        
#---------------------------------
    def Newton(self,I):
        x = I.c
        f1 = I.f1c
        f2 = I.f2c
        #dx = - f1/f2
        #f3 is known : use it in the first step:
        dx = (PF.sqrt(f2*f2 - pf(2)*f1*I.f3c) - f2)/I.f3c # out of two roots this one has a finite limit for f3->0
        while dx.abs() > pf("1.0e-16"):
            xold = x
            x = x + dx
            if x < I.x0  or x > I.x1:
                #print I, I.f1L, I.f1R, f1,f2, I.f3c
                #
                #xarray = [I.x0 + pf(k)*(I.x1-I.x0)/pf(20)   for k in range(21)]
                #
                ##yarray = [self.funlist[I.l].value(xi) for xi in xarray]
                #
                #y1array = [self.funlist[I.l].valueder(xi) for xi in xarray]
                #y1appr_array = [I.f1c + I.f2c*(xi-I.c)+ I.f3c * (xi-I.c)*(xi-I.c)/pf(2) for xi in xarray]
                #
                #xnumpy = numpy.array([PF.to_double(xi) for xi in xarray])
                ##ynumpy = numpy.array([PF.to_double(xi) for xi in yarray])
                #y1numpy = numpy.array([PF.to_double(xi) for xi in y1array])
                #y1appr_numpy = numpy.array([PF.to_double(xi) for xi in y1appr_array])
                #
                #plot( xnumpy,y1numpy)
                #plot( xnumpy,y1appr_numpy)
                #show()
                #
                #raise ValueError("Newton method brings out of the interval")
                
                #take a bisection step in the direction of zero of f1
                print ("Newton method brings out of the interval")
                if f1 < pf(0):
                    x = (xold + I.x1)/pf(2)
                else:
                    x = (xold + I.x0)/pf(2)
            f1 = self.funlist[I.l].valueder(x)
            f2 = self.funlist[I.l].valueder2(x)
            dx = - f1/f2
        I.xmin = x
        I.fmin = self.funlist[I.l].value(x)
        
               
    def plot_funlist(self, fun,d0,d1,label, title): 
        import numpy as np

        xarray = [d0 + pf(k)*(d1-d0)/pf(100)   for k in range(101)]
        yarray = [fun.value(xi) for xi in xarray]
        xnp = np.array([PF.to_double(xi) for xi in xarray])
        ynp = np.array([PF.to_double(xi) for xi in yarray])

        plt.subplots()
        plt.plot(xnp, ynp )
        plt.title(title)
        plt.savefig(label)
        plt.close()
        #plt.show()

