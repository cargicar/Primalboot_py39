#from lin_prog.lp_problem import precision
from config import precision

if precision:
    import solver.v1.mpfr_array.mpfr_array as NP_PREC
else:
    import numpy as NP_PREC
    
import scipy




############################################################################
#  
# LP_TABLE
# -------------
#
# The lp_cbtable class represents a CB table convolved with a particular
# value of sigma.  So a new instance is needed for each value of sigma.
#
# Right now initialized by passing sigma_table object
#
# vectab[l][k] - vector of convolved CB ders for spin l and dimension dlist[l][k] 
# funtab[l][i] - interpolation function for spin l and derivative component i 
#               (as a function of dimension)  
# vecfunlist[l] - interpolation data for spin l = (x_list, [[ y1 array], [y2 array],...]
#                where yi includes all derivatives (better for later dot blas multiplications)
#               
#############################################################################
# vecfunlist and unitCB are potentially at higher precision



class LP_Table:

    def __init__(self,sigma_table):
        self.dlist = sigma_table.dlist
        self.lmax = len(self.dlist)-1
        self.vectab = sigma_table.table
        self.CBlen = sigma_table.CBlen
        self.funtab = [[scipy.interpolate.splrep(self.dlist[l],
                        [self.vectab[l][d_index][i] for d_index in range(len(self.dlist[l]))],k=3)
                            for i in range(self.CBlen)]
                            if l%2 == 0 else None for l in range(self.lmax+1)]
        self.vecfunlist = [[   self.funtab[l][0][0], 
                               NP_PREC.array([[self.funtab[l][der][1][t_index] for der in range(self.CBlen)]
                                            for t_index in range(len(self.funtab[l][0][0])) ]) ]
                               if l%2 == 0 else None for l in range(self.lmax+1)]
        self.unitCB = NP_PREC.array(sigma_table.unitCB) 
        self.spacetimedim = 2*sigma_table.eps + 2
        
# ------------- END LP_TABLE

