/* LU decomposition code based on NR3 but turned into mpfr_t and C
See also NR2 where a different loop order is used  */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

#include "ludcmp_mpfr.h"
/* --------------------------------------------------
test block follows
----------------------------------------------------*/

 void print_matrix(mpfr_t* a, int n_rows, int n) 
{
    int i,j;
    printf("[");
    for (i=0; i<n_rows; i++) {
        printf("[");
        for (j=0; j<n; j++) {
            mpfr_out_str((FILE*)0,10,0,a[i*n+j],RND);
            printf(" ");
        }
        printf("]\n");
     }
     printf("]\n");
}
 

double randDouble()
{
   return ((double)rand()/(double)RAND_MAX);
}

int main(void)
{
    int i,j,k;
    const int n= 3;
    mpfr_t *a1, *a2, *a3, *ainv, *b;
    int* indx;
    mpfr_t sum,temp,temp1;
    
    printf("Working at %i bits of precision\n",PREC);
    
    mpfr_inits2(PREC, sum,temp,temp1,(mpfr_ptr)0);
    
    a1=(mpfr_t*)malloc(n*n*sizeof(mpfr_t));
    a2=(mpfr_t*)malloc(n*n*sizeof(mpfr_t));
    a3=(mpfr_t*)malloc(n*n*sizeof(mpfr_t));
    ainv=(mpfr_t*)malloc(n*n*sizeof(mpfr_t));
    b=(mpfr_t*)malloc(n*sizeof(mpfr_t));
    indx = (int*)malloc(n*sizeof(int));
    
    /* initialize matrices a1=a2 */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            mpfr_inits2(PREC,a1[i*n+j],a2[i*n+j],ainv[i*n+j],a3[i*n+j],(mpfr_ptr)0);
            mpfr_set_d(a1[i*n+j],randDouble(),RND);
            SET(a2[i*n+j],a1[i*n+j]);
        }
        mpfr_init2(b[i],PREC); /* b not used below */
        mpfr_set_d(b[i],randDouble(),RND);
    }
    printf ("Here's a %i x %i matrix:\n",n,n);
    print_matrix(a1,n,n);
    
    
    LUdcmp_in_place(a1,n,indx); /* a1 now holds the LU decomp, a2 is the original matrix */
    LUinverse(a1,n,indx,ainv);
    printf ("Here's the inverse computed via LU decomposition:\n",n,n);
    print_matrix(ainv,n,n);
    
    inverse(a2,n,a3);
    printf ("Here's the same inverse computed via a wrapper function:\n");
    print_matrix(a3,n,n);
    
    /* compute a3=a2.ainv -Identity */
     for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (i==j) 
                mpfr_set_si(sum,(long int)-1,RND);
            else 
                mpfr_set_zero(sum,1);
            for (k=0; k<n; k++) {
                MUL(temp,a2[i*n+k],ainv[k*n+j]);
                ADD(temp1,sum,temp);
                SET(sum,temp1);
            }
            SET(a3[i*n+j],sum);
        }
     }
     
     printf ("Here's A.A^{-1) - Id:\n");
     print_matrix(a3,n,n);
}


