from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

if os.name == 'posix':

	ext_modules = [Extension("medialaxis",
			["medialaxis.pyx", "skeletonizeCV.cpp"],
			language="c++",
			libraries=["opencv_core", "opencv_highgui", "opencv_imgproc"],
	)]

elif os.name == 'nt':

	#OPENCV_DIR = '.'
	#OPENCV_DIR = 'C:\\OpenCV2.2'
	OPENCV_DIR = 'C:\\opencv-2.4.6.0'
	#OPENCV_DIR = '/opt/local'

	ext_modules = [Extension("medialaxis",
			["medialaxis.pyx", "skeletonizeCV.cpp"],
			language="c++",
			#include_dirs=[OPENCV_DIR + '/include'],
			include_dirs=[OPENCV_DIR + '/build/include'],
			libraries=["opencv_core246", "opencv_highgui246"],
			#libraries=["opencv_core220", "opencv_highgui220"],
			#library_dirs=[OPENCV_DIR + '/lib'],
			library_dirs=[OPENCV_DIR + '/build/x86/vc9/lib'],
								)]

setup(
  name = 'Medial Axis',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

