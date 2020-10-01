SET( LocalWaveletPCAFiles
	BSGWFactory.cpp
	FaceData.cpp
	GlobalImageMatch.cpp
	NearestNeighborAssistant.cpp
	Profile.cpp
	StereoFace.cpp
	WaveletShapeSampler.cpp
	BSGWFactory.h
	BSplineGridWavelet.h
	FaceData.h
	GlobalImageMatch.h
	NearestNeighborAssistant.h
	NearestNeighborEnergy.h
	Profile.h
	StereoFace.h
	WaveletShapeFitter.h
	WaveletShapeSampler.h
)

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CLAPACK_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GLEW_INCLUDE_DIR})

ADD_LIBRARY(LocalWaveletPCALib ${LocalWaveletPCAFiles})
TARGET_LINK_LIBRARIES(LocalWaveletPCALib ${ANN_LIBRARIES} ${CLAPACK_LIBRARIES} ${OpenCV_LIBS} ${GLUT_LIBRARY} ${GLEW_LIBRARY})

ADD_EXECUTABLE(LocalWavletPCA LocalWavletPCA.cpp)
TARGET_LINK_LIBRARIES(LocalWavletPCA LocalWaveletPCALib)



	