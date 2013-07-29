from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

if os.name == 'posix':

	OGRE_INC = '/usr/local/include/OGRE'
	OIS_INC = '/usr/include/ois'

	ext_modules = [Extension("ogreprobe",
			["ogreprobe.pyx", "BaseApplication.cpp", "ProbeApp.cpp"],
			language="c++",
			include_dirs=[OGRE_INC, OIS_INC],
			libraries=["OgreMain","OIS"],
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
  name = 'Ogre3D Snake Probe',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

