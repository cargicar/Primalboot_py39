# cython: profile=False
# internal methods for Taylor_minmax_est

def candList(double a1, double a2, double a3, double x0, double x1):
    """Private method for Taylor_minmax_est
    Returns the set of possible extrema points for the cubic polynomial 
    const + a1(x-xm)+a2(x-xm)^2+a3(x-xm)^3 on the interval [x0,x1] where xm=(x0+x1)/2
    This set includes endpoints of the interval and the zeros of the derivative, 
    if they lie within the interval"""
    xm = (x0 + x1)/2.
    t = [x0, x1]
    disc = a2 * a2 - 3 * a3 * a1
    if a3 != 0:
        if disc >= 0:
            SQRT = pow(disc,0.5) #np.sqrt(disc)
            t += [x for x in 
                [xm + (-a2 - SQRT)/(3*a3),\
                xm + (-a2 + SQRT)/(3*a3)]
                if x0 < x < x1]
    elif a2 != 0:
        x = xm - a1/(2*a2)
        if  x0 < x < x1:
            t += [x]  
    return t
 
#----------------------------------------------------------  
#def cubic(a0, a1, a2, a3, x):
#    return a0 + x*(a1 + x*(a2 + x*a3))
#----------------------------------------------------------  
#cimport cython

#@cython.profile(False)
def cubic(double a0, double a1, double a2, double a3, double x):
    return a0 + x*(a1 + x*(a2 + x*a3))

    
     
