import tables.tables as cb
reload(cb)
from lin_prog.spectrum import Spectrum

print "Problem spec called."
threeD=True

def specfunc2d(sig, eps):
    spectrum = Spectrum(spacetimedim = 2,
                        deltamax = 41.5,
                        lmax = 40,
                        de = eps   )
    return spectrum

def specfunc3d(sig, eps):
    spectrum = Spectrum(spacetimedim = 3,
                        deltamax = 41.5,
                        lmax = 40,
                        de = eps   )
    return spectrum

if threeD:
    tab = cb.CB_Table(FILE = "../data/CBeps0.5n8m1.txt")
    specfunc=specfunc3d
else:
    tab = cb.CB_Table(FILE = "../data/CBeps0.0n8m1.txt")
    specfunc=specfunc2d


