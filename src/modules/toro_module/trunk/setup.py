from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules = [Extension("toromod", ["toromod.pyx", "toro.cpp", "posegraph2.cpp", "treeoptimizer2.cpp"],
			language="c++",
			extra_compile_args=['/EHsc']
			)]
#extra_link_args=['-framework', 'CoreFoundation']

#ext_modules = [Extension("toromod", ["toromod.pyx", "toro.cpp", "posegraph2.cpp", "treeoptimizer2.cpp", "toro.cpp"], language="c++")]
#ext_modules = [Extension("transform", ["transform.pyx"]), Extension("stability", ["stability.pyx"], language="c++")]

setup(
  name = 'TORO Optimizer',
  ext_modules = cythonize(ext_modules, compiler_directives={'language_level': "3"})
)

