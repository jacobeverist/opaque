from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

#-I$PHYSX_DIR\SDKs\Foundation\include -I$PHYSX_DIR\SDKs\Physics\include -I$PHYSX_DIR\SDKs\Cooking\include -I$PHYSX_DIR\SDKs\PhysXLoader\include -I$PHYSX_DIR\SDKs\NxCharacter\include

#userinc=distutils.sysconfig.get_python_inc(prefix=os.path.expanduser('~'))
#lib_dirs=['/usr/lib','/usr/local/lib',os.path.expanduser('~/lib')]

#include_dirs=['include']
#print userinc
#print lib_dirs

PHYSX_DIR = os.environ.get("PHYSX_DIR")


ext_modules = [Extension("nxsnakeprobe",
								["nxsnakeprobe.pyx", "NxSnake.cpp", "Stream.cpp"],
								language="c++",
								include_dirs=[PHYSX_DIR + '\SDKs\Foundation\include', PHYSX_DIR + '\SDKs\Physics\include', PHYSX_DIR + '\SDKs\Cooking\include', PHYSX_DIR + '\SDKs\PhysXLoader\include', PHYSX_DIR + '\SDKs\NxCharacter\include', PHYSX_DIR + '\TrainingPrograms\Programs\Shared_Source', PHYSX_DIR + '\TrainingPrograms\Programs\Shared_Source\User_Reports', PHYSX_DIR + '\TrainingPrograms\Programs\Shared_Source\Timer', PHYSX_DIR + '\TrainingPrograms\Programs\Shared_Source\win', PHYSX_DIR + '\Graphics\include\win32'],
								libraries=["PhysXLoader", "PhysXCooking", "NxCharacter"],
								library_dirs=[PHYSX_DIR + '\SDKs\lib\Win32', PHYSX_DIR + '\Graphics\lib\win32\glut'],
								)]

#ext_modules = [Extension("transform", ["transform.pyx"]), Extension("stability", ["stability.pyx"], language="c++")]

#print ds.get_config_vars()
setup(
  name = 'PhysX Snake Probe',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

