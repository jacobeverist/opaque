# The name of our project is "HELLO". CMakeLists files in this project can 
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and 
# to the root binary directory of the project as ${HELLO_BINARY_DIR}. 
cmake_minimum_required (VERSION 2.6) 
project (SKELETON) 


IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	include_directories (C:/OpenCV2.2/include)
	link_directories (C:/OpenCV2.2/lib)
	add_executable (cvSkel skeletonizeCV.c) 
	target_link_libraries (cvSkel opencv_core220 opencv_highgui220 opencv_imgproc220) 
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	include_directories (/opt/local/include)
	link_directories (/opt/local/lib)
	add_executable (cvSkel skeletonizeCV.c) 
	target_link_libraries (cvSkel opencv_core.2.3.1 opencv_highgui.2.3.1 opencv_imgproc.2.3.1) 
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")



