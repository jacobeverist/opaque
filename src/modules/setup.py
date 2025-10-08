from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = [Extension("transform", ["transform.pyx"]),
			Extension("stability", ["stability.pyx", "ValueStability.cpp"], language="c++"),
			Extension("reference", ["reference.pyx", "RefNode.cpp"], language="c++"),
			Extension("servo", ["servo.pyx"]),
			Extension("func", ["func.pyx"]),
			Extension("icp", ["icp.pyx"])]

setup(
  name = 'Joint Transformations',
  ext_modules = cythonize(ext_modules, compiler_directives={'language_level': "3"})
)

