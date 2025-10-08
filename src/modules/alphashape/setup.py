from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os


if os.name == 'posix':
	ext_modules = [Extension("alphamod", ["alphamod.pyx", "alpha2.cpp"],
				language="c++",
		)]

elif os.name == 'nt':


	ext_modules = [Extension("alphamod", ["alphamod.pyx", "alpha2.cpp"],
	                        language="c++",
				extra_compile_args=['/EHsc'],
				include_dirs=['CGAL-4.2/include', 'boost_1_53_0', 'CGAL-4.2/auxiliary/gmp/include'],
				library_dirs=['CGAL-4.2/lib', '../../binaries/boost_1_53_0/lib']
				)]


#extra_link_args=['-framework', 'CoreFoundation']


setup(
  name = 'Alpha Shape 2D',
  ext_modules = cythonize(ext_modules, compiler_directives={'language_level': "3"})
)

