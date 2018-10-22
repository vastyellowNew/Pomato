cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project( PMATE )

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(XAVIER32 FALSE)
  set(XAVIER64 TRUE)
else (CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(XAVIER32 TRUE)
  set(XAVIER64 FALSE)
endif (CMAKE_SIZEOF_VOID_P MATCHES "8")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release Public"
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/xavier/cmake" )

option( BUILD_OPTION_OPENMP "use openmp parallelization" OFF)

set( LIBMODE SHARED )
set( ALGO_LIBS "")

#######################################################################

# find_package( OpenGL REQUIRED )
find_package( Boost   1.46.1  )
find_package( Teem   REQUIRED )
# find_package( GLUT   REQUIRED )
find_package( Eigen  REQUIRED )
find_package( VTK )
# find_package( GLUI )

set( XAVIER_DIR ${CMAKE_SOURCE_DIR}/xavier )
set( NVIS_DIR ${CMAKE_SOURCE_DIR}/nvis )

# build configuration

include_directories( ${CMAKE_SOURCE_DIR} )
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${NVIS_DIR} )
include_directories( ${Eigen_INCLUDE_DIR} )
# include_directories( ${GLUT_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR}/GL)
# include_directories( ${GLUI_INCLUDE_DIR} )
# include_directories( ${OPENGL_INCLUDE_DIR} )
include_directories( ${Teem_INCLUDE_DIR} )
include_directories( ${XAVIER_DIR} )

link_directories( ${CMAKE_BINARY_DIR}/lib )
link_directories( ${Boost_LIBRARY_DIRS} )

set( CMAKE_CXX_FLAGS_DEBUG   "-g " )
set( CMAKE_CXX_FLAGS_RELEASE "-O2 -funroll-loops -fstrict-aliasing -DNODEBUG " )
set( CMAKE_CXX_FLAGS_PUBLIC "-O2 -funroll-loops -fstrict-aliasing -DNDEBUG -DNODEBUG " )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-result -Wno-deprecated-declarations" )

if(BUILD_OPTION_OPENMP)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

if( BUILD_PMATE_PROFILE )
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -pg")
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -pg")
endif()

#######################################################################
# set build product destinations

set( LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib )
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )

add_subdirectory( xavier )

######################################################################
# determine targets to build and define their linking
######################################################################

IF( NOT PMATE_INSTALL_PRESET_IS_SET )
  SET( CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Install path prefix, prepended onto install directories." FORCE )
  SET( PMATE_INSTALL_PRESET_IS_SET 1 CACHE INTERNAL "Default install path written" FORCE )
ENDIF()

install( DIRECTORY bin DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY scripts DESTINATION . USE_SOURCE_PERMISSIONS)
install( DIRECTORY utils/GLSLShaderCode DESTINATION utils USE_SOURCE_PERMISSIONS)

######################################################################
# install all binaries
######################################################################
foreach( target ${PMATE_TARGETS} )
    install( TARGETS ${target}
             RUNTIME DESTINATION bin
    )
endforeach()
