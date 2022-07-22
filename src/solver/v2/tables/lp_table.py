############################################################################
#  
# LP_TABLE
# -------------
#
# The lp_cbtable class represents a CB table convolved with a particular
# value of sigma.  So a new instance is needed for each value of sigma.
#
# This is a new version using polynomials in rho-repr

import solver.v2.rho_repr.cb as cb
import solver.v2.rho_repr.cbdata as cbdata
import solver.v2.mpfr_array.mpfr_array as mpfr_array
import solver.v2.prec_float.prec_float as PF
import utils.logging as log
import version 

def pf(x):
    return PF.prec_float(x,prec=212)


class LP_Table:

    # ds should be a pf
    def __init__(self,filename, ds=None): #initialized by reading from file
        #self.spacetimedim=3
        with open(filename,"r") as FILE:
            # file has ds hard-coded (old file format)
            if ds is None:
                self.eps = PF.prec_float(FILE.readline(), prec=212) #spacetime dimension (d-2)/2
                self.lmax = int(FILE.readline()) # maximal spin
                self.vecfunlist = [ cb.CB(FILE)
                        if l%2 == 0 else None for l in range(self.lmax+1)]
                # total length of vectors to be fed into linear programming
                self.CBlen = self.vecfunlist[0].CBlen 
                self.prec = self.vecfunlist[0].prec # precision
                self.unitCB = mpfr_array.array([ FILE.readline() for i in range(self.CBlen)], 
                        prec=self.prec) # unit CB (contracted)
            else:
                log.stats('about read in cb file')
                cbd=cbdata.CBdata(FILE)
                log.stats('read in cb file')
                    
                cbd.set_ds(ds)
                log.stats('set ds')
                self.lmax=cbd.lmax
                self.eps=cbd.eps
                self.vecfunlist = [ cb.CB(cbd.CBs[l])
                        if l%2 == 0 else None for l in range(self.lmax+1)]
                log.stats('generated vecfunlist')
                # total length of vectors to be fed into linear programming
                self.CBlen = self.vecfunlist[0].CBlen 
                self.prec = self.vecfunlist[0].prec # precision
                self.unitCB = cbd.unitCB()
            self.spacetimedim= pf(2) * self.eps + pf(2)
            
            
        
# ------------- END LP_TABLE

