/* polynomial evaluation   */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>

#include "polyval.h"

/* compute the value of the polynomial (p0), its first (p1), second (p2) and third (p3)
derivative using Horner's method */
/* the polynomial is assumed stored as an mpfr_t array starting from zeroth term.
delta0 is the expansion point */
void p0( mpfr_t delta0, mpfr_t *p, int n, mpfr_t x, mpfr_t ret)
/* delta0 - polynomial expansion point
 * p polynomial coefficient array
 * n length of p
 * x - value of dimension at which p is needed
 * ret variable in which result is returned */{
    int i;
    mpfr_t temp;
    mpfr_init2(temp,mpfr_get_prec(x));
    SUB(temp,x,delta0);
    SET(ret, p[n-1]);
    for(i=n-2;i>=0;i--) {
        MUL(ret, ret, temp);
        ADD(ret, ret, p[i]);
    }
    mpfr_clears(temp,(mpfr_ptr)0);
}

void p1( mpfr_t delta0, mpfr_t *p, int n, mpfr_t x, mpfr_t ret) {
    int i;
    mpfr_t diff, temp;
    mpfr_init2(temp,mpfr_get_prec(x));
    mpfr_init2(diff,mpfr_get_prec(x));
    SUB(diff,x,delta0);
    mpfr_mul_ui(ret, p[n-1], n-1,MPFR_RNDN);
    for(i=n-2;i>=1;i--) {
        MUL(ret, ret, diff);
        mpfr_mul_ui(temp, p[i], i,MPFR_RNDN);
        ADD(ret, ret, temp);
    }
    mpfr_clears(temp,diff,(mpfr_ptr)0);
}

void p2( mpfr_t delta0, mpfr_t *p, int n, mpfr_t x, mpfr_t ret) {
    int i;
    mpfr_t diff, temp;
    mpfr_init2(temp,mpfr_get_prec(x));
    mpfr_init2(diff,mpfr_get_prec(x));
    SUB(diff,x,delta0);
    mpfr_mul_ui(ret, p[n-1], (n-1)*(n-2),MPFR_RNDN);
    for(i=n-2;i>=2;i--) {
        MUL(ret, ret, diff);
        mpfr_mul_ui(temp, p[i], i*(i-1),MPFR_RNDN);
        ADD(ret, ret, temp);
    }
    mpfr_clears(temp,diff,(mpfr_ptr)0);
}

void p3( mpfr_t delta0, mpfr_t *p, int n, mpfr_t x, mpfr_t ret) {
    int i;
    mpfr_t diff, temp;
    mpfr_init2(temp,mpfr_get_prec(x));
    mpfr_init2(diff,mpfr_get_prec(x));
    SUB(diff,x,delta0);
    mpfr_mul_ui(ret, p[n-1], (n-1)*(n-2)*(n-3),MPFR_RNDN);
    for(i=n-2;i>=3;i--) {
        MUL(ret, ret, diff);
        mpfr_mul_ui(temp, p[i], i*(i-1)*(i-2),MPFR_RNDN);
        ADD(ret, ret, temp);
    }
    mpfr_clears(temp,diff,(mpfr_ptr)0);
}

/* compute the value of the rational function (r0), its first (r1), second (r2) and third (r3)
derivative */
/* poles - positions of simple poles
 * c - residues of the poles
 * n number of poles
 * x - value of dimension at which r is needed
 * ret variable in which result is returned */
void r0( mpfr_t *poles, mpfr_t *c, int n, mpfr_t x, mpfr_t ret ) {
    int i;
    mpfr_t tmp;
    mpfr_init2(tmp,mpfr_get_prec(x));
    
    mpfr_set_zero(ret,1);
    for(i=0;i<n;i++) {
        SUB(tmp,x,poles[i]);
        DIV(tmp,c[i],tmp);
        ADD(ret, ret, tmp);
    }
    mpfr_clears(tmp,(mpfr_ptr)0);
}

void r1( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfr_t x, mpfr_t r1) {
    int i;
    mpfr_t diff, tmp;
    mpfr_init2(diff,mpfr_get_prec(x));
    mpfr_init2(tmp,mpfr_get_prec(x));
    
    mpfr_set_zero(r1,1);
    for(i=0;i<n;i++) {
        SUB(diff,x,poles[i]);
        
        DIV(tmp,polecoeffs[i],diff);
        
        DIV(tmp,tmp,diff);
        SUB(r1, r1, tmp);
    }
    mpfr_clears(tmp,diff,(mpfr_ptr)0);
}

void r2( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfr_t x, mpfr_t r2) {
    int i;
    mpfr_t diff, tmp;
    mpfr_init2(diff,mpfr_get_prec(x));
    mpfr_init2(tmp,mpfr_get_prec(x));
    
    mpfr_set_zero(r2,1);
    for(i=0;i<n;i++) {
        SUB(diff,x,poles[i]);
        
        DIV(tmp,polecoeffs[i],diff);
        DIV(tmp,tmp,diff);
        DIV(tmp,tmp,diff);
        
        mpfr_mul_si(tmp,tmp,2,MPFR_RNDN);
        ADD(r2, r2, tmp);
    }
    mpfr_clears(tmp,diff,(mpfr_ptr)0);
}

void r3( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfr_t x, mpfr_t r2) {
    int i;
    mpfr_t diff, tmp;
    mpfr_init2(diff,mpfr_get_prec(x));
    mpfr_init2(tmp,mpfr_get_prec(x));
    
    mpfr_set_zero(r2,1);
    for(i=0;i<n;i++) {
        SUB(diff,x,poles[i]);
        
        DIV(tmp,polecoeffs[i],diff);
        DIV(tmp,tmp,diff);
        DIV(tmp,tmp,diff);
        DIV(tmp,tmp,diff);
        
        mpfr_mul_si(tmp,tmp,6,MPFR_RNDN);
        SUB(r2, r2, tmp);
    }
    mpfr_clears(tmp,diff,(mpfr_ptr)0);
}

/* same as r0-r3 but for all derivatives up to order 3 simultaneously */
void r0123( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfr_t x, mpfr_t r0,mpfr_t r1,mpfr_t r2,mpfr_t r3) {
    int i;
    mpfr_t diff, tmp;
    mpfr_init2(diff,mpfr_get_prec(x));
    mpfr_init2(tmp,mpfr_get_prec(x));
    
    mpfr_set_zero(r0,1);
    mpfr_set_zero(r1,1);
    mpfr_set_zero(r2,1);
    mpfr_set_zero(r3,1);
    
    for(i=0;i<n;i++) {
        SUB(diff,x,poles[i]);
        
        DIV(tmp,polecoeffs[i],diff);
        ADD(r0, r0, tmp);
        
        DIV(tmp,tmp,diff);
        
        SUB(r1, r1, tmp);
        
        DIV(tmp,tmp,diff);
        mpfr_mul_si(tmp,tmp,2,MPFR_RNDN);
        
        ADD(r2, r2, tmp);
        
        DIV(tmp,tmp,diff);
        mpfr_mul_si(tmp,tmp,3,MPFR_RNDN);
        
        SUB(r3, r3, tmp);
    }
    mpfr_clears(tmp,diff,(mpfr_ptr)0);
}

/* computes derivatoves of order 1,2,3 of rho^x (p(x)+r(x)) */
/* rho, logrho - passed from outside
 * delta0,p,n - same parameters as in p0
 * poles,polecoeffs, npoles - same params as in r0
 * d1,d2,d3 - where the answer is returned */
void valueder123(mpfr_t rho, mpfr_t logrho, \
                 mpfr_t delta0, mpfr_t *p, int n, \
                 mpfr_t *poles, mpfr_t *polecoeffs, int npoles, \
                 mpfr_t x, \
                 mpfr_t d1, mpfr_t d2, mpfr_t d3) {
        mpfr_t p0val,p1val,p2val,p3val,r0val,r1val,r2val,r3val, pr0, pr1, pr2, pr3, tmp, fac;
        
        mpfr_inits2(mpfr_get_prec(x),p0val,p1val,p2val,p3val,r0val,r1val,r2val,r3val, \
                    pr0, pr1, pr2, pr3, tmp, fac,(mpfr_ptr)0);
        p0(delta0, p, n, x, p0val);
        p1(delta0, p, n, x, p1val);
        p2(delta0, p, n, x, p2val);
        p3(delta0, p, n, x, p3val);
        r0123(poles, polecoeffs, npoles, x, r0val, r1val, r2val, r3val);
        
        SET(fac,rho);
        mpfr_pow(fac,fac,x,MPFR_RNDN);
        
        SET(pr0,r0val);
        ADD(pr0,pr0,p0val);
        
        SET(pr1,r1val);
        ADD(pr1,pr1,p1val);
        
        SET(pr2,r2val);
        ADD(pr2,pr2,p2val);
        
        SET(pr3,r3val);
        ADD(pr3,pr3,p3val);
        
        SET(d1,pr0);
        MUL(d1,d1,logrho);
        ADD(d1,d1,pr1);
        MUL(d1,d1,fac);
        
        SET(d2,pr0);
        MUL(d2,d2,logrho);
        mpfr_mul_si(tmp,pr1,2,MPFR_RNDN);
        ADD(d2,d2,tmp);
        MUL(d2,d2,logrho);
        ADD(d2,d2,pr2);
        MUL(d2,d2,fac);
        
        SET(d3,pr0);
        MUL(d3,d3,logrho);
        mpfr_mul_si(tmp,pr1,3,MPFR_RNDN);
        ADD(d3,d3,tmp);
        MUL(d3,d3,logrho);
        mpfr_mul_si(tmp,pr2,3,MPFR_RNDN);
        ADD(d3,d3,tmp);
        MUL(d3,d3,logrho);
        ADD(d3,d3,pr3);
        MUL(d3,d3,fac);
        
        mpfr_clears(p0val,p1val,p2val,p3val,r0val,r1val,r2val,r3val, pr0, pr1, pr2, pr3, tmp, fac,(mpfr_ptr)0);
                
        //return (
        //    (self.rho**opdim) * (p1 + r1 + (p0 + r0) * self.logrho),
        //    (self.rho**opdim) * (p2 + r2
        //                            + PF.prec_float(2, prec = self.prec) * self.logrho * (p1 + r1)
        //                            + (p0 + r0) * self.logrho * self.logrho),
        //    (self.rho**opdim) * (p3 + r3
        //                            + PF.prec_float(3, prec = self.prec) * self.logrho * (p2 + r2)
        //                            + PF.prec_float(3, prec = self.prec) * self.logrho * self.logrho * (p1 + r1)
        //                            + (p0 + r0) * self.logrho * self.logrho * self.logrho)
        //    )
}













/* INTERVAL FUNCTIONS FOLLOW (CURRENTLY 16.07.2013 NO LONGER USED)
/* compute interval bounds on the value of the polynomial using Horner's method */
/* the polynomial is assumed stored as an mpfr_t array starting from zeroth term */
/* NB: no care for memory leak below */
void i_p0(mpfr_t delta0, mpfr_t *p, int n, mpfi_t inter, mpfi_t ret) {
    int i;
    mpfi_t diff;
    mpfi_init2(diff,mpfr_get_prec(p[0]));
    mpfi_sub_fr(diff,inter,delta0);
    iSETfr(ret, p[n-1]);
    for(i=n-2;i>=0;i--) {
        iMUL(ret, ret, diff);
        iADDfr(ret, ret, p[i]);
    }
}

/* compute the interval value of the first derivative */
void i_p1(mpfr_t delta0, mpfr_t *p, int n, mpfi_t inter, mpfi_t ret) {
    int i;
    mpfr_t temp;
    mpfi_t diff;
    mpfi_init2(diff,mpfr_get_prec(p[0]));
    mpfi_sub_fr(diff,inter,delta0);
    mpfr_init2(temp,mpfr_get_prec(p[n-1]));
    mpfr_mul_ui(temp,p[n-1],n-1,MPFR_RNDN);
    iSETfr(ret, temp);
    for(i=n-2;i>=1;i--) {
        iMUL(ret, ret, diff);
        mpfr_mul_ui(temp,p[i],i,MPFR_RNDN);
        iADDfr(ret, ret, temp);
    }
    mpfr_clear(temp);
}

void i_p01(mpfr_t delta0, mpfr_t *p, int n, mpfi_t inter, mpfi_t ret0, mpfi_t ret1) {
    int i;
    mpfr_t temp;
    mpfi_t diff;
    mpfi_init2(diff,mpfr_get_prec(p[0]));
    mpfi_sub_fr(diff,inter,delta0);
    
    iSETfr(ret0, p[n-1]);
    
    mpfr_init2(temp,mpfr_get_prec(p[n-1]));
    mpfr_mul_ui(temp,p[n-1],n-1,MPFR_RNDN);
    iSETfr(ret1, temp);
    
    for(i=n-2;i>=1;i--) {
        iMUL(ret0, ret0, diff);
        iADDfr(ret0, ret0, p[i]);
    
        iMUL(ret1, ret1, diff);
        mpfr_mul_ui(temp,p[i],i,MPFR_RNDN);
        iADDfr(ret1, ret1, temp);
    }
    
    iMUL(ret0, ret0, diff);
    iADDfr(ret0, ret0, p[0]);
    
    mpfr_clear(temp);
}


/* intervalizes rational function:*/
void i_r0( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfi_t x, mpfi_t ret) {
    int i;
    mpfi_t tmp;
    mpfi_init2(tmp,mpfr_get_prec(polecoeffs[0]));
    mpfi_set_str (tmp, "0", 10);
    for(i=0;i<n;i++) {
        mpfi_sub_fr(tmp,x,poles[i]);
        mpfi_fr_div(tmp,polecoeffs[i],tmp);
        iADD(ret, ret, tmp);
    }
}

void i_r1( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfi_t x, mpfi_t ret) {
    int i;
    mpfi_t tmp;
    mpfi_init2(tmp,mpfr_get_prec(polecoeffs[0]));
    mpfi_set_str (tmp, "0", 10);
    for(i=0;i<n;i++) {
        mpfi_sub_fr(tmp,x,poles[i]);
        iMUL(tmp,tmp,tmp);
        mpfi_fr_div(tmp,polecoeffs[i],tmp);
        iSUB(ret, ret, tmp);
    }
}

void i_r01( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfi_t x, mpfi_t r0, mpfi_t r1) {
    int i;
    mpfi_t diff, tmp;
    mpfi_init2(diff,mpfr_get_prec(polecoeffs[0]));
    mpfi_init2(tmp,mpfr_get_prec(polecoeffs[0]));
    
    mpfi_set_str (r0, "0", 10);
    mpfi_set_str (r1, "0", 10);
   
    for(i=0;i<n;i++) {
        mpfi_sub_fr(diff,x,poles[i]);
        
        mpfi_fr_div(tmp,polecoeffs[i],diff);
        iADD(r0, r0, tmp);
        
        iDIV(tmp,tmp,diff);
        iSUB(r1, r1, tmp);
    }
}

/* this function computes logrho* r0+r1 */
void i_r1together( mpfr_t *poles, mpfr_t *polecoeffs, int n, mpfr_t logrho, mpfi_t x, mpfi_t ret) {
    int i;
    mpfi_t invdiff, tmp;
    mpfi_init2(invdiff,mpfr_get_prec(polecoeffs[0]));
    mpfi_init2(tmp,mpfr_get_prec(polecoeffs[0]));
    
    mpfi_set_str (ret, "0", 10);
    
    for(i=0;i<n;i++) {
        mpfi_sub_fr(invdiff,x,poles[i]);
        mpfi_ui_div(invdiff,1,invdiff); /* now invdiff contains 1/(x-pole[i]) */
        
        mpfi_fr_sub(tmp, logrho, invdiff);
        mpfi_mul_fr(tmp,tmp,polecoeffs[i]);
        iMUL(tmp,invdiff,tmp);
        iADD(ret, ret, tmp);
    }
}