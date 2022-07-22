#!/bin/sh

rm ./solver/v1/mpfr_array/mpfr_array.so
rm ./solver/v1/prec_float/prec_float.so
rm ./solver/v1/branch_bound/Taylor_methods.so
rm ./solver/v2/mpfr_array/mpfr_array.so
rm ./solver/v2/prec_float/prec_float.so
rm ./solver/v2/rho_repr/cb.so
rm ./solver/v2/rho_repr/cbdata.so

rm -r ./solver/v1/branch_bound/build
rm -r ./solver/v1/mpfr_array/build/
rm -r ./solver/v1/prec_float/build/
rm -r ./solver/v2/rho_repr/build/
rm -r ./build/
