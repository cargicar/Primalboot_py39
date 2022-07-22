#from lp_problem import precision
from config import precision

if precision:
    import solver.v1.mpfr_array.mpfr_array as NP_PREC
    from solver.v1.mpfr_array.mpfr_array import to_double
else:
    import numpy as NP_PREC
    def to_double (x): return x
from solver.v1.branch_bound.bb_problem import BB_Problem, COMP_TOLERANCE

from scipy import interpolate

#import branch_bound.bb_problem as bb_problem
#reload(bb_problem)
###############################
#
#  contract_minimize function
#  takes:
#  d0, d1, rho, and a list of interpolated functions (one per each der.index)
#  functions must be cubic and have the same knots
#
##############################
#np.seterr(all='raise')

def contract_minimize_rel(d0, d1, rho, vecfun,
                        threshold=0.,greed_factor=2., relative_diff_goal=0.01 , 
                        debug=False):
    
    contract_list = to_double (NP_PREC.dot(vecfun[1],rho)) # list of doubles holding contraction
    
    tck = (vecfun[0], contract_list , 3)   

    if d1-d0 < COMP_TOLERANCE: # interval consisting of one point
        return [d0, interpolate.splev (d0, tck), None ]
    
    bbprob = BB_Problem(tck, d0, d1, "local")
   
    bbprob.findmin_relative_negative(threshold=threshold, greed_factor=greed_factor, 
                                relative_diff_goal=relative_diff_goal)
        
    I = bbprob.ilist[bbprob.iL]
    if debug:
        return [I.x_min, interpolate.splev (I.x_min, tck), bbprob] #send bbprob back for debugging
    else:
        return [I.x_min, interpolate.splev (I.x_min, tck)] 
