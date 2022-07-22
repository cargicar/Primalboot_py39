import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


ext2 = Extension(
    "prec_float",                 
    ['prec_float.pyx'],    
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
    ext_modules = [ext2],
)


