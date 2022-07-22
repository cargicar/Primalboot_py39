from libc.stdio cimport stdout, FILE, printf

from c_mpfr cimport *

cdef extern from "mpfi.h":
    ctypedef struct mpfi_t:
        pass
    void mpfi_init2(mpfi_t, mpfr_prec_t)
    void mpfi_interv_si(mpfi_t, long int, long int)
    void mpfi_interv_d(mpfi_t, double, double)
    void mpfi_interv_fr(mpfi_t, mpfr_t, mpfr_t)
    void mpfi_neg(mpfi_t, mpfi_t)
    void mpfi_add(mpfi_t, mpfi_t, mpfi_t)
    void mpfi_sub(mpfi_t, mpfi_t, mpfi_t)
    void mpfi_mul(mpfi_t, mpfi_t, mpfi_t)
    void mpfi_div(mpfi_t, mpfi_t, mpfi_t)
    int mpfi_exp(mpfi_t, mpfi_t)
    int mpfi_log(mpfi_t, mpfi_t)
    int mpfi_get_left(mpfr_t, mpfi_t)
    int mpfi_get_right(mpfr_t, mpfi_t)
    size_t mpfr_out_str(FILE*, int, size_t, mpfr_t, mpfr_rnd_t)
    void mpfi_clear (mpfi_t)
    int mpfi_set_str (mpfi_t rop, char *s, int base)

