#include <stdio.h>
#include "mpfr_funcs.h"
#include <gmp.h>
#include <mpfr.h>

void testmpfr()
{
  unsigned int i;
  mpfr_t s, t, u;

  mpfr_init2 (t, 200);
  mpfr_set_d (t, 1, MPFR_RNDD);
  mpfr_init2 (s, 200);
  mpfr_set_d (s, 1, MPFR_RNDD);
  mpfr_init2 (u, 200);
  for (i = 1; i <= 100; i++)
    {
      mpfr_mul_ui (t, t, i, MPFR_RNDU);
      mpfr_set_d (u, 1.0, MPFR_RNDD);
      mpfr_div (u, u, t, MPFR_RNDD);
      mpfr_add (s, s, u, MPFR_RNDD);
    }
  printf ("\n MPFR sum is ");
  mpfr_out_str (stdout, 10, 0, s, MPFR_RNDD);
  putchar ('\n');
  mpfr_clear (s);
  mpfr_clear (t);
  mpfr_clear (u);
  mpfr_t x;
  mpfr_init2(x, 200);
  mpfr_set_str(x, "1.234234223", 10, MPFR_RNDD);
  printf ("\n set MPFR is ");
  mpfr_out_str (stdout, 10, 0, x, MPFR_RNDD);
  printf("\n");
  //return 0;
}


/* ------------------ OLD ------------------------ */


double* printarray(int len, double* arr) {
	int i;
	for (i=0; i < len; i++) {
		printf("%.13f ", arr[i]);
		arr[i] +=1.;
	}
	printf("\n");
	return arr;
}


