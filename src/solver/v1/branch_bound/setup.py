from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

Ext_modules = [Extension(
            "Taylor_methods", 
            ["Taylor_methods.pyx"],
            extra_compile_args=['-O3'])]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = Ext_modules
)

#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
#from Cython.Build import cythonize

#ext_modules = [Extension("bb_problem", ["profile_lp_c.pyx"])]

#setup(
#    cmdclass = {'build_ext': build_ext},
#    ext_modules = cythonize(ext_modules)
#)


