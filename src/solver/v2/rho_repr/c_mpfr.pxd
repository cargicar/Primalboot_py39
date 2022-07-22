from libc.stdio cimport stdout, FILE, printf, sprintf

cdef extern from "mpfr.h":
    ctypedef struct mpfr_t:
        pass
    ctypedef short mpfr_prec_t
    ctypedef enum mpfr_rnd_t:
        MPFR_RNDN=0,   #round to nearest, with ties to even */
        MPFR_RNDZ,     #round toward zero */
        MPFR_RNDU,     #round toward +Inf */
        MPFR_RNDD,     #round toward -Inf */
        MPFR_RNDA,     #round away from zero */
        MPFR_RNDF,     #faithful rounding (not implemented yet) */
        MPFR_RNDNA=-1  #round to nearest, with ties away from zero (mpfr_round) */

    void mpfr_init2(mpfr_t, mpfr_prec_t)
    void mpfr_set_d(mpfr_t, double, mpfr_rnd_t)
    double mpfr_get_d(mpfr_t,mpfr_rnd_t)
    void mpfr_set(mpfr_t, mpfr_t, mpfr_rnd_t)
    void mpfr_div(mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t)
    void mpfr_mul(mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t)
    void mpfr_add(mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t)
    void mpfr_sub(mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t)
    int mpfr_neg (mpfr_t, mpfr_t, mpfr_rnd_t)
    size_t mpfr_out_str(FILE*, int, size_t, mpfr_t, mpfr_rnd_t)
    int mpfr_init_set_str (mpfr_t x, char *s, int base, mpfr_rnd_t rnd)
    int mpfr_set_str (mpfr_t x, char *s, int base, mpfr_rnd_t rnd)
    int mpfr_set_si (mpfr_t, long int, mpfr_rnd_t)
    void mpfr_set_zero (mpfr_t, int) # int=1 to set to +0
    int mpfr_snprintf(char *buf, size_t n, char *template, ...)
    int mpfr_sprintf(char *buf, char *template, ...)
    int mpfr_cmp_si (mpfr_t, long int)
    int mpfr_cmp (mpfr_t, mpfr_t)
    void mpfr_clear (mpfr_t)
    int mpfr_pow_si (mpfr_t, mpfr_t, long int, mpfr_rnd_t)
    int mpfr_pow (mpfr_t, mpfr_t, mpfr_t, mpfr_rnd_t)
    int mpfr_log (mpfr_t, mpfr_t, mpfr_rnd_t)
    int mpfr_sqrt (mpfr_t, mpfr_t, mpfr_rnd_t)
    mpfr_prec_t mpfr_get_prec (mpfr_t)
    int mpfr_log(mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd)
    int mpfr_ui_div (mpfr_t rop, unsigned long int op1, mpfr_t op2, mpfr_rnd_t rnd)
    int mpfr_mul_si (mpfr_t rop, mpfr_t op1, long int op2, mpfr_rnd_t rnd)
   