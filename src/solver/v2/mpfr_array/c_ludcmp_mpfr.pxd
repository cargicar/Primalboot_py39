from c_mpfr cimport *

cdef extern from "ludcmp_mpfr.h":
    int LUdcmp_in_place( mpfr_t *, int , int *)
    int LUsolve_in_place(mpfr_t *, int , int* , mpfr_t *)
    int LUinverse(mpfr_t *, int , int*, mpfr_t *)
    int inverse(mpfr_t *, int , mpfr_t *, mpfr_t *)
    int LUdcmp( mpfr_t *, int , mpfr_t *, int *)
    int ULsolve_in_place(mpfr_t *, int , int* , mpfr_t *)