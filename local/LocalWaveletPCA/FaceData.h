///////////////////////////////////////////////////////////////////////////////
//
//	FaceData.h
//
//	Header file for the CFaceData class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __FACEDATA_H__
#define __FACEDATA_H__


#include "StereoFace.h"

#include <vector>

class CFaceData
{
public:

	static int*								m_pIndices;
	static int								m_nIndices;

	static double*							m_pAvgVertices;
	static int								m_nVertices;

	static double*							m_pAvgNormals;

	static double							m_centroidX, m_centroidY, m_centroidZ;

	const static int						m_iTemplateAlignNoseTip = 0;
	const static int						m_iTemplateAlignNoseBridge = 1;
	const static int						m_iTemplateAlignOutsideLeftEyeCorner = 2;
	const static int						m_iTemplateAlignOutsideRightEyeCorner = 3;
	const static int						m_nTemplateAlignMinIndices = 4;

	static bool init(const char* cstrFaceFileName);

	static void saveFaceOFFIndices(const char* cstrFileName, int nVertices, float* pVertices, int nIndices, unsigned int* pIndices, int nPolygonVertices, unsigned int indexIgnoreValue = (unsigned int)(-1));

	static void drawFaceWithTexCoords(double* pVertices, double* pTexCoord, double* pNormals = NULL);

	static void empericalPlaneFromLandmarkIndices(double* pNormal, int* pLandmarkIndices, int numIndices, double* pmatTransform);

	static void stereographicProject(double* pVertices, int nVertices, double* pmatTransform, double* pVertStereo, double* pXmin, double* pYmin, double* pXmax, double* pYmax);

	static bool readOFFMeshAsPatch(char *szFilename, std::vector<SStereoPatch> &patches, float g_nInDataScale);
};


#endif //__FACEDATA_H__


