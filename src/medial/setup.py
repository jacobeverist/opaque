from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import os
from Cython.Distutils import build_ext

OPENCV_DIR = 'C:\\OpenCV2.2'
#OPENCV_DIR = '/opt/local'

#ext_modules = [Extension("medialaxis",
#		["medialaxis.pyx", "skeletonizeCV.cpp"],
#		language="c++",
#		include_dirs=[OPENCV_DIR + '/include'],
#		libraries=["opencv_core.2.3.1", "opencv_highgui.2.3.1", "opencv_imgproc.2.3.1"],
#		library_dirs=[OPENCV_DIR + '/lib'],

ext_modules = [Extension("medialaxis",
		["medialaxis.pyx", "skeletonizeCV.cpp"],
		language="c++",
		include_dirs=[OPENCV_DIR + '/include'],
		libraries=["opencv_core220", "opencv_highgui220"],
		library_dirs=[OPENCV_DIR + '/lib'],
								)]

setup(
  name = 'Medial Axis',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

