
import sys, os, stat
# import commands
#import subprocess # Carlos Edit
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# we'd better have Cython installed, or it's a no-go
try:
    from Cython.Distutils import build_ext
except:
    print ("You don't seem to have Cython installed. Please get a")
    print ("copy from www.cython.org and install it")
    sys.exit(1)


    
extensions = [
    Extension(
    "solver.v2.rho_repr.cb",                 
    ['solver/v2/rho_repr/cb.pyx', 'solver/v2/rho_repr/polyval.c'],    
    include_dirs=['/sw/include','.'],  
    library_dirs =['/sw/lib'],     
    libraries=["mpfr", "mpfi"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3']
    )
    ,
    Extension(
    "solver.v2.rho_repr.cbdata",                 
    ['solver/v2/rho_repr/cbdata.pyx', 'solver/v2/rho_repr/polyval.c'],    
    include_dirs=['/sw/include','.'],  
    library_dirs =['/sw/lib'],     
    libraries=["mpfr", "mpfi"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3']
    )
    ,
    Extension(
    "solver.v2.prec_float.prec_float",                 
    ['solver/v2/prec_float/prec_float.pyx'],    
    include_dirs=['/sw/include','.'],  
    library_dirs =['/sw/lib'],     
    libraries=["mpfr", "mpfi"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3']
    )
    ,
    Extension(
    "solver.v2.mpfr_array.mpfr_array",                 
    ['solver/v2/mpfr_array/mpfr_array.pyx','solver/v2/mpfr_array/ludcmp_mpfr.c'],    
    include_dirs=['/sw/include','.'],  
    library_dirs =['/sw/lib'],     
    libraries=["mpfr", "mpfi"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3']
    ),
    Extension(
    "solver.v1.prec_float.prec_float",                 
    ['solver/v1/prec_float/prec_float.pyx'],    
    include_dirs=['/sw/include','.'],  
    library_dirs =['/sw/lib'],     
    libraries=["mpfr", "mpfi"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3']
    )
    ,
    Extension(
    "solver.v1.mpfr_array.mpfr_array",                 
    ['solver/v1/mpfr_array/mpfr_array.pyx','solver/v1/mpfr_array/ludcmp_mpfr.c'],    
    include_dirs=['/sw/include','.'],  
    library_dirs =['/sw/lib'],     
    libraries=["mpfr", "mpfi"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3']
    ),
    Extension(
        "solver.v1.branch_bound.Taylor_methods", 
        ["solver/v1/branch_bound/Taylor_methods.pyx"],
        extra_compile_args=['-O3']
    )
    ]

# finally, we can pass all this to distutils
setup(
  name="bootstrap",
  packages=["solver.v2.rho_repr", "solver.v2.prec_float",
      "solver.v2.mpfr_array", "solver.v1.prec_float", "solver.v1.mpfr_array",
      "solver.v1.branch_bound.Taylor_methods"],
  ext_modules=extensions,
  cmdclass = {'build_ext': build_ext},
)
