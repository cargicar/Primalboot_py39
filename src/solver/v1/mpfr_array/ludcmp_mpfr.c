/* LU decomposition code based on NR3 but turned into mpfr_t and C
See also NR2 where a different loop order is used  */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>

#include "ludcmp_mpfr.h"

void pr(mpfr_t a)
{
    mpfr_out_str((FILE*)0,10,0,a,RND);
}

int LUdcmp_in_place( mpfr_t *lu, int n, int *indx) {
/* performs LU decomposition in place, destroying the initial matrix */
	int i,imax,j,k;
	mpfr_t big,temp,temp1,temp2;
	mpfr_t* vv;
	
	mpfr_inits2(PREC, big,temp,temp1,temp2,(mpfr_ptr)0);
	
	vv = (mpfr_t *) malloc(n*sizeof(mpfr_t)); 
	if (vv==NULL) {
	        fprintf(stderr, "Memory allocation error in LUdcmp");
            return(0);
	}
	

	for (i=0;i<n;i++) {
		mpfr_set_zero(big, 1);
		
		for (j=0;j<n;j++) {
		    ABS(temp, lu[i*n+j]);
			if ( GTR(temp,big) ) SET(big,temp);
		}
		if ( ISZERO(big)) { 
		    fprintf(stderr, "Singular matrix in LUdcmp (check 1)");
            return(0);
		}
		mpfr_init2(vv[i],PREC);
		INVERSE(vv[i], big);
	}
	for (k=0;k<n;k++) {
		mpfr_set_zero(big, 1);
		for (i=k;i<n;i++) {
		    ABS(temp1,lu[i*n+k]);
			MUL(temp,vv[i],temp1);
			if ( GTR(temp,big) ) {
				SET(big,temp);
				imax=i;
			}
		}
		if (k != imax) {
			for (j=0;j<n;j++) {
				SET(temp,lu[imax*n+j]);
				SET(lu[imax*n+j],lu[k*n+j]);
				SET(lu[k*n+j],temp);
			}
			SET(vv[imax],vv[k]);
		}
		indx[k]=imax;
		if (ISZERO(lu[k*n+k])) {
		    fprintf(stderr, "Singular matrix in LUdcmp (check 2)");
            return(0);
		}
		for (i=k+1;i<n;i++) {
			SET(temp,lu[i*n+k]); 
			DIV(lu[i*n+k],temp,lu[k*n+k]);
			SET(temp,lu[i*n+k]);
			for (j=k+1;j<n;j++) {
			    MUL(temp1,temp,lu[k*n+j]);
			   /* SUB(temp2,lu[i*n+j], temp1);
				SET(lu[i*n+j],temp2); */
				SUB(lu[i*n+j],lu[i*n+j], temp1);
			}
		}
	}
	/* clear memory */
	for (i=0;i<n;i++) 
		mpfr_clear(vv[i]);
	free(vv);
	mpfr_clears(big,temp,temp1,temp2,(mpfr_ptr)0);
	return(1);
}


int LUdcmp( mpfr_t *a, int n,  mpfr_t *lu, int *indx) 
/* returns lu decomposition in lu */
/* a is not destroyed */
{
    int i,j;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            SET(lu[i*n+j],a[i*n+j]);
        }
    }
     
    if(LUdcmp_in_place(lu, n, indx))
        return 1;
    else return 0;    
}

int LUsolve_in_place(mpfr_t *lu, int n, int* indx, mpfr_t *x)
/* x is the RHS and also the solution is returned in x*/
{
	int i,ii=0,ip,j;
	mpfr_t sum,temp,temp1;
	mpfr_inits2 (PREC, sum, temp,temp1, (mpfr_ptr)0);
	for (i=0;i<n;i++) {
		ip=indx[i];
		SET(sum, x[ip]);
		SET(x[ip],x[i]);
		if (ii != 0)
			for (j=ii-1;j<i;j++) {
			    MUL(temp,lu[i*n+j],x[j]);
			    SUB(temp1, sum, temp);
			    SET(sum, temp1);
			 }
		else if (! ISZERO(sum))
			ii=i+1;
		SET(x[i],sum);
	}
	for (i=n-1;i>=0;i--) {
		SET(sum,x[i]);
		for (j=i+1;j<n;j++) {
		     MUL(temp,lu[i*n+j],x[j]);
			 SUB(temp1, sum, temp);
			 SET(sum, temp1);
		}
		DIV(x[i],sum,lu[i*n+i]);
	}
	mpfr_clears (sum, temp,temp1, (mpfr_ptr)0);
	return 1;
}


int ULsolve_in_place(mpfr_t *lu, int n, int* indx, mpfr_t *x)
/* solves the system A^T y = x where A = PLU with LU decomposition stored in lu
and P permutation matrix stored in vector indx (row i is permuted with row indx[i])
x is the RHS and also the solution is returned in x*/
{
	int i,j,temp;
	mpfr_t sum,tmp;
	mpfr_t* xcopy;
	mpfr_inits2 (PREC, sum, tmp, (mpfr_ptr)0);
	int* piv;
	for (i=0;i<n;i++) {
		SET(sum,x[i]);
		for (j=0;j<i;j++) 
		    {   
		        MUL(tmp, lu[j*n+i], x[j]);
		        SUB(sum,sum,tmp);
		    }
		DIV(x[i],sum,lu[i*n+i]);
	}
	for (i=n-1;i>=0;i--) {
		SET(sum,x[i]);
		for (j=i+1;j<n;j++) 
		 {   
		        MUL(tmp, lu[j*n+i], x[j]);
		        SUB(sum,sum,tmp);
		    }
		SET(x[i],sum);
	}
	/* now it remains to unscramble row permutation which is a bit subtle since it's stored as a sequence
	pairwise permutations */
	/* so we apply a sequence of these permutations to the vector piv initially set to (0,1,2,3...n-1)
	to find out where each row ends up */ 
	piv = (int*) malloc(n*sizeof(int));
	for (i=0;i<n;i++) 
	    piv[i]=i;
    for (i=0;i<n;i++) {
        temp=piv[i];
        piv[i]=piv[indx[i]];
        piv[indx[i]]=temp;
        }
	/* now unscramble*/
	xcopy = (mpfr_t *) malloc(n*sizeof(mpfr_t)); 
	for (i=0;i<n;i++) {
	    mpfr_inits2 (PREC, xcopy[i], (mpfr_ptr)0);
	    SET(xcopy[i], x[i]); 
	    }
    	for (i=0;i<n;i++) { 
	    SET(x[piv[i]],xcopy[i]);
	    }
	    
	/* free memory */
	for (i=0;i<n;i++) {
	    mpfr_clear (xcopy[i]);
	    }
	
	mpfr_clears (sum, tmp, (mpfr_ptr)0);
	free(xcopy); 
	free(piv);
	return 1;
}


int LUinverse(mpfr_t *lu, int n, int* indx, mpfr_t *y)
/* y is the inverse matrix */
{
    int i,j;
    mpfr_t * col;
    
    col = (mpfr_t *) malloc(n*sizeof(mpfr_t)); 
	if (col==NULL) {
	        fprintf(stderr, "Memory allocation error in LUinverse");
            return(0);
	}
	for(i=0;i<n;i++) 
	    mpfr_init2(col[i],PREC);
    for(j=0;j<n;j++) {
        for(i=0;i<n;i++) 
           mpfr_set_zero(col[i],1);
        mpfr_set_si(col[j],(long int)1,RND);
        LUsolve_in_place(lu,n,indx,col); 
        for(i=0;i<n;i++) 
            SET(y[i*n+j],col[i]);
    }
    for(i=0;i<n;i++) 
	    mpfr_clear(col[i]);
    free(col);
    return(1);
}   

int inverse(mpfr_t *a, int n, mpfr_t *y)
/* returns y as the inverse matrix */
/* a is not destroyed */
{
    int i,j;
    mpfr_t *lu;
    int *indx;
   
    lu=(mpfr_t*)malloc(n*n*sizeof(mpfr_t));
    indx = (int*)malloc(n*sizeof(int));
    if (lu==NULL || indx == NULL) {
	        fprintf(stderr, "Memory allocation error in inverse");
            return(0);
	}
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            mpfr_init2(lu[i*n+j],PREC);
            SET(lu[i*n+j],a[i*n+j]);
        }
    }
     
    if( LUdcmp_in_place(lu, n, indx) && LUinverse(lu, n, indx, y) ){
        for(i=0;i<n;i++) {
            for(j=0;j<n;j++) {
                mpfr_clear(lu[i*n+j]);
            }
        }
        free(lu);
        free(indx);
        return 1;
    }
    else return 0;    
}


