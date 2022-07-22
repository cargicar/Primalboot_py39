#cython: profile=False
#cdef extern from "stdio.h":
#    cdef struct FILE:
#        pass
#    cdef FILE* stdout
#    int printf(const char *format, ...)

from c_mpfr cimport *
from c_mpfi cimport *
    
cdef extern from "polyval.h":
    void p0( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t)
    void p1( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t)
    void p2( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t)
    void p3( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t)
    void i_p0( mpfr_t, mpfr_t *, int, mpfi_t, mpfi_t)
    void i_p1( mpfr_t, mpfr_t *, int, mpfi_t, mpfi_t)
    void i_p01(mpfr_t delta0, mpfr_t *p, int n, mpfi_t inter, mpfi_t ret0, mpfi_t ret1)
    
    void r0( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t)
    void r1( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t)
    void r2( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t)
    void r3( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t)
    void r0123( mpfr_t *, mpfr_t *, int , mpfr_t, mpfr_t,mpfr_t,mpfr_t,mpfr_t)
    void i_r0( mpfr_t *, mpfr_t *,int, mpfi_t, mpfi_t)
    void i_r1( mpfr_t *, mpfr_t *,int, mpfi_t, mpfi_t)
    void i_r01( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfi_t x, mpfi_t r0, mpfi_t r1)
    void i_r1together( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfr_t logrho, mpfi_t x, mpfi_t ret)
    
    void valueder123(mpfr_t rho, mpfr_t logrho, mpfr_t delta0, mpfr_t *p, int n,
                 mpfr_t *poles, mpfr_t *polecoeffs, int npoles, mpfr_t x, mpfr_t d1, mpfr_t d2, mpfr_t d3)



