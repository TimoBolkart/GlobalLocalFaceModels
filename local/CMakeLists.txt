
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(LocalWaveletPCA)

INCLUDE(ConfigureANN.cmake)
INCLUDE(ConfigureCLAPACK.cmake)
INCLUDE(ConfigureGLEW.cmake)

FIND_PACKAGE(OpenCV)
IF(OpenCV_FOUND)
  INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
  LINK_DIRECTORIES(${OpenCV_LIBRARY_DIRS})
ELSE(OpenCV_FOUND)
  Message("OpenCV not found. Please install OpenCV or set OpenCV_ROOT_DIR properly. Some applications and/or libraries need it to compile.")
ENDIF(OpenCV_FOUND)

FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
  LINK_DIRECTORIES(${GLUT_LIBRARY_DIR})
ELSE(GLUT_FOUND)
  Message("Warning. GLUT not found. Some applications and libraries will not be compiled.")
ENDIF(GLUT_FOUND)

IF(ANN_FOUND AND CLAPACK_FOUND AND OpenCV_FOUND AND GLUT_FOUND AND GLEW_FOUND)
  ADD_SUBDIRECTORY(LocalWaveletPCA)
ELSE(ANN_FOUND AND CLAPACK_FOUND AND OpenCV_FOUND AND GLUT_FOUND AND GLEW_FOUND)
  Message("PROBLEM: One of the required libraries not found. StereoFace will not be compiled.")
ENDIF(ANN_FOUND AND CLAPACK_FOUND AND OpenCV_FOUND AND GLUT_FOUND AND GLEW_FOUND)