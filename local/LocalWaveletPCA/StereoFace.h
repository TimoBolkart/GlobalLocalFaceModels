///////////////////////////////////////////////////////////////////////////////
//
//	StereoFace.h
//
//	Main header file for the StereoFace program
//	contains common declarations, definitions and includes
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __STEREOFACE_H__
#define __STEREOFACE_H__


#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include <string>
#include <vector>
#include <iostream>

#ifdef WIN32
#define WINDOWS_LEAN_AND_MEAN
#define NOMINMAX
#include <tchar.h>
#include "windows.h"
#else
#endif //WIN32

#include "glew.h"
#include "glut.h"

#include "cv.h"
#include "highgui.h"

//includes for random variables
#include "GaussianSample.h"

//preprocessor constants
#define STEREOFACE_CAPTURE					0
#define STEREOFACE_DEFAULT_STRING_LENGTH	2048

#define STEREOFACE_USECUDA					0

#define MAX_INFO_LOG_SIZE					4096

#pragma warning (disable : 4244)


//global constants
const double					g_epsilon = 1e-10;

//vector types, etc
#include "abutil.h"

//function prototypes
GLchar *LoadShaderText(const char *fileName);

//framebuffer functions
void fbInit(GLuint & nfb);
void fbTerm(GLuint & nfb);
void fbBind(GLuint & nfb);
void fbUnbind();
int fbMaxAttachments();
void fbAttachTexture(GLuint & nfb, GLenum attachment, GLenum texTarget, GLuint texId, int mipLevel = 0, int zSlice = 0);
void fbDetach(GLenum attachment, GLuint & nfb);
void fbDetachAll(GLuint & nfb);

//shader functions
GLuint loadVertexShader(const char* szFilename);
GLuint loadFragmentShader(const char* szFilename);
void programError(const char* szMessage, GLuint program);
void shaderError(const char* szMessage, GLuint shader);

//gl error checking
GLenum checkGLError(const char* szLastCommand = NULL);

//misc math functions
//void svdOrthogonalize(double* pmatR);

//inline helper functions
void cross_prod_3(double a[3], double b[3], double c[3]);
void normalize(double a[3]);
void mul_3x3_3x3(double* pOut, double* pA, bool bTransA, double* pB, bool bTransB);
void mul_4x4_4x4(double* pOut, double* pA, bool bTransA, double* pB, bool bTransB);

//structures
struct SStereoPatch
{
	float px, py, pz;					//position
	float nx, ny, nz;					//normal
	float tx, ty, tz;					//tangent
	float bx, by, bz;					//binormal
	float u, v;							//texture coordinates
	float r, g, b;						//static color information
	bool bLandmark;
};

struct SPatchDepth
{
	int iPatch;							//patch index
	float z;							//patch depth from viewing position
};

struct SBundlerKeypoint
{
	double x, y, z;
	unsigned char r, g, b;
	int nViews;
	int* pViews;
	double* puv;

	SBundlerKeypoint()
	{
		x = y = z = 0.0;
		r = g = b = 0;
		nViews = 0;
		pViews = NULL;
		puv = NULL;
	}
	~SBundlerKeypoint()
	{
		if (pViews != NULL)
			delete [] pViews;
		if (puv != NULL)
			delete [] puv;
	}

	void initViewList(int numViews)
	{
		nViews = numViews;
		pViews = new int[nViews];
		_ASSERT(pViews != NULL);
		puv = new double[nViews * 2];
		_ASSERT(puv != NULL);
	}
};


//functions for reading point cloud data
void readPointCloudFile(std::string strFile, std::vector<SStereoPatch>& patches, float dataScale, bool& bNormals);


#endif //__STEREOFACE_H__


