import setpath
import config

# which version of the code to use
version=2


if version == 1:
    from solver.v1.lin_prog.spectrum import Spectrum
    import solver.v1.lin_prog.lp_problem as lp_problem
    #from solver.v1.mpfr_array.mpfr_array import to_double1
    import solver.v1.mpfr_array.mpfr_array as NP_PREC
    import solver.v1.mpfr_array as mpfr_array
    import solver.v1.prec_float.prec_float as prec_float
else:
    from solver.v2.lin_prog.spectrum import Spectrum
    import solver.v2.lin_prog.lp_problem as lp_problem
    #from solver.v2.mpfr_array.mpfr_array import to_double2
    import solver.v2.mpfr_array.mpfr_array as NP_PREC
    import solver.v2.mpfr_array as mpfr_array
    import solver.v2.prec_float.prec_float as prec_float

def to_double(obj):
    if isinstance(obj, NP_PREC.ndarray):
        return NP_PREC.to_double(obj)
    elif isinstance(obj, prec_float.prec_float):
        return prec_float.to_double(obj)
    else:
        return obj

if config.precision is False:
    TD=np.array
else:
    TD=to_double


def pfwrap(val):
    if version==1:
        return val
    else:
        return prec_float.prec_float(val,prec=212)
