# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( alpha2 )

set(CGAL_DIR C:/Users/everist/workspace/opaque/src/alpha2/CGAL-4.1)
#LINK_DIRECTORIES(C:/Users/everist/workspace/opaque/src/alpha2/boost_1_53_0/stage/lib)

set(Boost_DIR c:/Users/everist/workspace/opaque/src/alpha2/boost_1_53_0)
set(BOOST_LIBRARYDIR c:/Users/everist/workspace/opaque/src/alpha2/boost_1_53_0/stage/lib)
set(BOOST_INCLUDEDIR c:/Users/everist/workspace/opaque/src/alpha2/boost_1_53_0)
#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_MULTITHREADED ON)
#set(WITH_GMP OFF)
#set(WITH_MPFR OFF

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
 
if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package


# Creating entries for target: alpha2
# ############################

add_executable( alpha2  alpha2.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS alpha2 )

# Link the executable to CGAL and third-party libraries
target_link_libraries(alpha2   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

