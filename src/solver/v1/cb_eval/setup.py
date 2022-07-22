from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("cb_eval_mpfr", ["cb_eval_mpfr.pyx"])]

setup(
  name = 'cb_eval_mpfr',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
