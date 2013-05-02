from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("alphamod", ["alphamod.pyx", "alpha2.cpp"],
			language="c++",
			extra_compile_args=['/EHsc'],
			include_dirs=['CGAL-4.1/include', 'boost_1_53_0', 'CGAL-4.1/auxiliary/gmp/include'],
			library_dirs=['CGAL-4.1/lib', 'boost_1_53_0/stage/lib']
			)]
#extra_link_args=['-framework', 'CoreFoundation']

#ext_modules = [Extension("toromod", ["toromod.pyx", "toro.cpp", "posegraph2.cpp", "treeoptimizer2.cpp", "toro.cpp"], language="c++")]
#ext_modules = [Extension("transform", ["transform.pyx"]), Extension("stability", ["stability.pyx"], language="c++")]

setup(
  name = 'Alpha Shape 2D',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

