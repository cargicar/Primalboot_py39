/* polynomial evaluation   */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>


// MPFR functions
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

// MPFI functions
#define iSET(X, Y, Z)  mpfi_set(X, Y, Z)
#define iSETfr(X, Y)  mpfi_set_fr(X, Y)
#define iMUL(X,Y,Z)    mpfi_mul(X,Y,Z)
#define iDIV(X,Y,Z)    mpfi_div(X,Y,Z)
#define iSUB(X,Y,Z)    mpfi_sub(X,Y,Z)
#define iADD(X,Y,Z)    mpfi_add(X,Y,Z)
#define iADDfr(X,Y,Z)    mpfi_add_fr(X,Y,Z)

 
void p0( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t);
void p1( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t);
void p2( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t);
void p3( mpfr_t, mpfr_t *, int, mpfr_t, mpfr_t);
void i_p0( mpfr_t, mpfr_t *, int, mpfi_t, mpfi_t);
void i_p1( mpfr_t, mpfr_t *, int, mpfi_t, mpfi_t);
void i_p01(mpfr_t, mpfr_t *, int , mpfi_t , mpfi_t, mpfi_t );

void r0( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t);
void r1( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t);
void r2( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t);
void r3( mpfr_t *, mpfr_t *,int, mpfr_t, mpfr_t);
void r0123( mpfr_t *, mpfr_t *, int , mpfr_t, mpfr_t,mpfr_t,mpfr_t,mpfr_t);
void i_r0( mpfr_t *, mpfr_t *,int, mpfi_t, mpfi_t);
void i_r1( mpfr_t *, mpfr_t *,int, mpfi_t, mpfi_t);
void i_r01( mpfr_t *, mpfr_t *, int , mpfi_t , mpfi_t , mpfi_t );
void i_r1together( mpfr_t *, mpfr_t *, int n, mpfr_t , mpfi_t , mpfi_t);
void r0123( mpfr_t *, mpfr_t *, int , mpfr_t, mpfr_t,mpfr_t,mpfr_t,mpfr_t);
void valueder123(mpfr_t, mpfr_t, mpfr_t, mpfr_t *, int, \
                 mpfr_t *, mpfr_t *, int , mpfr_t , mpfr_t , mpfr_t , mpfr_t );

