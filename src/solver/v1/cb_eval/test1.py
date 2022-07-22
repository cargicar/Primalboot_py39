import cb_eval
import cb_eval_mpfr 
import mpfr.prec_float as PF
reload(cb_eval_mpfr)
import numpy as np
from timing import t0, deltat

eps=1
nmax=10
mmax=1
#def flatten (List):
#    [item for sublist in List for item in sublist]

deltamax=40


t0()
res = [[cb_eval_mpfr.cb_ders(eps,l,opdim,10,1,prec=212) for opdim in np.arange(l+2*eps,deltamax,0.05)] for l in range(0,2,2)]
print deltat()

