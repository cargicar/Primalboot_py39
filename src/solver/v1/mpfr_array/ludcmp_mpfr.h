/* LU decomposition code based on NR3 but turned into mpfr_t and C
See also NR2 where a different loop order is used  */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

#define PREC 212
#define RND MPFR_RNDN
#define ABS(X, Y)   mpfr_abs(X, Y, RND)
#define SET(X, Y)   mpfr_set(X, Y, RND)
#define GTR(X,Y)    (mpfr_cmp(X,Y) > 0)  
#define ISZERO(X)   mpfr_zero_p(X) 
#define INVERSE(X, Y)  mpfr_si_div(X , (long int) 1, Y, RND)
#define MUL(X,Y,Z)    mpfr_mul(X,Y,Z,RND)
#define DIV(X,Y,Z)    mpfr_div(X,Y,Z,RND)
#define SUB(X,Y,Z)    mpfr_sub(X,Y,Z,RND)
#define ADD(X,Y,Z)    mpfr_add(X,Y,Z,RND)

 
int LUdcmp( mpfr_t *, int , mpfr_t *, int *); 
int LUdcmp_in_place( mpfr_t *, int , int *); 
int LUsolve_in_place(mpfr_t *, int , int* , mpfr_t *);
int ULsolve_in_place(mpfr_t *, int , int* , mpfr_t *);

int LUinverse(mpfr_t *, int , int*, mpfr_t *);
int inverse(mpfr_t *, int , mpfr_t *);

