import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#sourcefiles = ['mpfr_wrap.pyx', 'mpfr_funcs.c']

ext1 = Extension(
    "mpfr_wrap",                 
    ['mpfr_wrap.pyx', 'mpfr_funcs.c'],    
    include_dirs=['/sw/include',numpy.get_include()],   
    library_dirs =['/sw/lib'],   
    libraries=["mpfr"],             
    #extra_link_args=[...],       # if needed
    #cmdclass = {'build_ext': build_ext}
    )

ext2 = Extension(
    "mpfr_array",                 
    ['mpfr_array.pyx'],    
    include_dirs=['/sw/include',numpy.get_include()],  
    library_dirs =['/sw/lib'],       
    libraries=["mpfr"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3'],
    #cmdclass = {'build_ext': build_ext}
    )

setup(
    cmdclass = {'build_ext': build_ext},
    include_dirs = ['/sw/include',numpy.get_include(), "qd"],
    #library_dirs =['sw/lib'],
    #ext_modules = [Extension("example", sourcefiles)],
    ext_modules = [ext1, ext2],
)


