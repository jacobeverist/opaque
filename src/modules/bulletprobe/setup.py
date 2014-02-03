from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

if os.name == 'posix':

	BULLET_DIR = '/usr/local/include/bullet'

	ext_modules = [Extension("bulletsnakeprobe",
			["bulletsnakeprobe.pyx", "BulletSnake.cpp"],
			language="c++",
			#extra_compile_args=['-fPIC'],
			include_dirs=[BULLET_DIR],
			libraries=["BulletDynamics", "LinearMath", "BulletCollision"],
								)]

elif os.name == 'nt':

	BULLET_DIR = '.\\bullet-2.81'

	ext_modules = [Extension("bulletsnakeprobe",
			["bulletsnakeprobe.pyx", "BulletSnake.cpp"],
			language="c++",
			include_dirs=[BULLET_DIR + '\\src'],
			#libraries=["BulletDynamics_Debug", "LinearMath_Debug", "BulletCollision_Debug"],
			libraries=["BulletDynamics", "linearmath", "BulletCollision"],
			#libraries=["BulletDynamics_vs2008_debug", "LinearMath_vs2008_debug", "BulletCollision_vs2008_debug"],
			#libraries=["BulletDynamics_vs2008", "LinearMath_vs2008", "BulletCollision_vs2008"],
			#library_dirs=[BULLET_DIR + '\\lib'],
			#library_dirs=[BULLET_DIR + '\\lib\\Release'],
			library_dirs=[BULLET_DIR + '\\lib\\Debug'],
								)]

setup(
  name = 'Bullet Snake Probe',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

