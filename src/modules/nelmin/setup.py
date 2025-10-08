from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# if GPU
#ext_modules = [Extension("nelminICP", ["nelminICP.pyx"],
#		library_dirs = ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\\v4.2\lib\Win32'],
#		libraries = ['cuda', 'cudart','python26'],
#		extra_link_args=['/NODEFAULTLIB:libcmt'],
#		extra_objects=["runNelmin.lib"])]

# else

ext_modules = [Extension("nelminICP",
		["nelminICP.pyx", "runHostNel.c"],
language='c++',
		)] 




setup(
  name = 'nelmin optimization',
  ext_modules = cythonize(ext_modules, compiler_directives={'language_level': "3"})
)

