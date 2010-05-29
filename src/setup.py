from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("transform", ["transform.pyx"])]

setup(
  name = 'Joint Transformations',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

