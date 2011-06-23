from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("transform", ["transform.pyx"]),
			Extension("stability", ["stability.pyx", "ValueStability.cpp"], language="c++"),
			Extension("reference", ["reference.pyx", "RefNode.cpp"], language="c++"),
			Extension("servo", ["servo.pyx"]),
			Extension("icp", ["icp.pyx"])]
#ext_modules = [Extension("transform", ["transform.pyx"]), Extension("stability", ["stability.pyx"], language="c++")]

setup(
  name = 'Joint Transformations',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

