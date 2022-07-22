from c_mpfr cimport *
from c_mpfi cimport *

cdef class prec_float:
    cdef mpfr_t data
    cdef public int prec

cdef class prec_interval:
    cdef mpfi_t data
    cdef public int prec
   