import config
import utils.logging as log

#########
#
#   here define a class which specifies the problem to solve
#   1) spacetime dimension
#   2) spectrum whose consistency we have to check
#   for the moment it's just an array of intervals 
#
##########

class SpecInt:
    def __init__(self,l, deltamin, deltamax):
        if deltamin > deltamax:
            errstr="specint.__init__: deltamin > deltamax (%d: %f > %f)"\
                    %(l, deltamin, deltamax)
            log.critical(errstr)
            raise ValueError(errstr)
            #sys.exit("specint.__init__: deltamin > deltamax")
        self.l = l
        self.d0 = deltamin
        self.d1 = deltamax
    def __repr__(self):
        return "<l=%s (%s %s)>" % (self.l, self.d0, self.d1)

class Spectrum:
    def __init__(self, spacetimedim =2, deltamax=42., lmax=40., de=1.0):    
        self.spacetimedim = spacetimedim
        self.lmax=lmax
        self.deltamax=deltamax
        self.ilist = [ SpecInt(l, l+spacetimedim-2. if l>0 else de , deltamax)
                    for l in range(0,lmax+1,2) ] 
    
    def setgaps(self,gaps):
        self.ilist = [ SpecInt(l, l+self.spacetimedim-2.+config.fudge if l not in gaps.keys() else
            gaps[l], self.deltamax)
                    for l in range(0,self.lmax+1,2) ]

    def add_interval(self, l, d0, d1):
        self.ilist += [SpecInt(l, d0, d1)]

    def __repr__(self):
        rpr=''
        for il in self.ilist:
            rpr+= str(il)+','
        return rpr
