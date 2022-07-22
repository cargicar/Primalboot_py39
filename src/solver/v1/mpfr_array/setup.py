import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

#sourcefiles = ['mpfr_wrap.pyx', 'mpfr_funcs.c']

ext1 = Extension(
    "mpfr_array_test",                 
    ['mpfr_array_test.pyx', 'mpfr_funcs.c'],    
    #include_dirs=[numpy.get_include()],         
    libraries=["mpfr"],             
    include_dirs=['/home/selshowk/Utils/local/include','/sw/include',numpy.get_include()],   
    library_dirs =['/home/selshowk/Utils/local/lib','/sw/lib'],   
    #extra_link_args=[...],       # if needed
    cmdclass = {'build_ext': build_ext}
    )

ext2 = Extension(
    "mpfr_array",                 
    ['mpfr_array.pyx', 'ludcmp_mpfr.c'],    
    #include_dirs=[numpy.get_include()],       
    include_dirs=['/home/selshowk/Utils/local/include','/sw/include',numpy.get_include()],   
    library_dirs =['/home/selshowk/Utils/local/lib','/sw/lib'],   
    libraries=["mpfr"],             
    #extra_link_args=[...],       # if needed
    extra_compile_args=['-O3'],
    cmdclass = {'build_ext': build_ext}
    )

setup(
    cmdclass = {'build_ext': build_ext},
    #include_dirs = [numpy.get_include(), "qd"],
    include_dirs = ['/sw/include',numpy.get_include(), "qd"],
    #ext_modules = [Extension("example", sourcefiles)],
    ext_modules = [ext1, ext2],
)


