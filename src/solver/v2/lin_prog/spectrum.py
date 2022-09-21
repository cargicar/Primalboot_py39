#########
#
#   here define a class which specifies the problem to solve
#   1) spacetime dimension
#   2) spectrum whose consistency we have to check
#   for the moment it's just an array of intervals 
#
##########
from solver.v2.prec_float.prec_float import prec_float

def pf(x):
 #  """ Converts a string or a float into a big float input=x: String or Float output= prec_float : mpfr float with precision prec """
       return prec_float(x,prec=212)

class SpecInt:
    def __init__(self,l, deltamin, deltamax):
        if deltamin > deltamax:
            sys.exit("specint.__init__: deltamin > deltamax")
        self.l = l
        self.d0 = deltamin
        self.d1 = deltamax
    def __repr__(self):
        return "<l=%s (%s %s)>" % (self.l, self.d0, self.d1)

class Spectrum:
    def __init__(self, spacetimedim = pf(2), deltamax=pf(42), lmax=40, de=pf(1)):    
        self.spacetimedim = spacetimedim
        self.ilist = [ SpecInt(l, pf(l)+spacetimedim-pf(2) if l>0 else de , deltamax)
                    for l in range(0,lmax+1,2) ] 
    

