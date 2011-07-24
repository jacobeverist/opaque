from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

OPENCV_DIR = 'C:\\OpenCV2.2'

ext_modules = [Extension("medialaxis",
		["medialaxis.pyx", "skeletonizeCV.cpp"],
		language="c++",
		include_dirs=[OPENCV_DIR + '\\include'],
		libraries=["opencv_core220", "opencv_highgui220"],
		library_dirs=[OPENCV_DIR + '\\lib'],
								)]

setup(
  name = 'Medial Axis',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

