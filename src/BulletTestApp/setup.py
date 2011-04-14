from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

BULLET_DIR = 'C:\\Users\\everist\\Desktop\\bullet-2.77'


#C:\Users\everist\Desktop\bullet-2.77\build\lib\Debug

ext_modules = [Extension("bulletsnakeprobe",
		["bulletsnakeprobe.pyx", "BulletSnake.cpp"],
		language="c++",
		include_dirs=[BULLET_DIR + '\\src'],
		libraries=["BulletDynamics", "linearmath", "BulletCollision"],
		library_dirs=[BULLET_DIR + '\\lib\\Debug'],
								)]

setup(
  name = 'Bullet Snake Probe',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

