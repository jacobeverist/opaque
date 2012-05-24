from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# if GPU
"""
ext_modules = [Extension("nelminICP", ["nelminICP.pyx"],
		library_dirs = ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\\v4.2\lib\Win32'],
		libraries = ['cuda', 'cudart','python26'],
		extra_link_args=['/NODEFAULTLIB:libcmt'],
		extra_objects=["runNelmin.lib"])]
"""

# else

ext_modules = [Extension("nelminICP", ["nelminICP.pyx", "runHostNel.c"])] 





setup(
  name = 'nelmin optimization',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

