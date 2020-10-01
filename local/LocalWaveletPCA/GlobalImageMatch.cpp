///////////////////////////////////////////////////////////////////////////////
//
//	GlobalImageMatch.cpp
//
//	Source file for the CGlobalImageMatch class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#include "GlobalImageMatch.h"


///////////////////////////////////////////////////////////////////////////////
//member functions
///////////////////////////////////////////////////////////////////////////////

//void CGlobalImageMatch::initTextures()
//{
//	int i;
//	uchar* pData;
//	int step;
//
//	m_pntexImages = new GLuint[m_nImages];
//	_ASSERT(m_pntexImages != NULL);
//	glGenTextures(m_nImages, m_pntexImages);
//
////	m_pntexDisparity = new GLuint[m_nImages];
////	_ASSERT(m_pntexDisparity != NULL);
////	glGenTextures(m_nImages, m_pntexDisparity);
//	m_pntexDepth = new GLuint[m_nImages];
//	_ASSERT(m_pntexDepth != NULL);
//	glGenTextures(m_nImages, m_pntexDepth);
//
//	for (i = 0; i < m_nImages; i++)
//	{
//		cvGetRawData(m_ppImages[i], &pData, &step);
//
//		glBindTexture(GL_TEXTURE_2D, m_pntexImages[i]);
//		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
//		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_ppImages[i]->width, m_ppImages[i]->height, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, pData);
//
////		glBindTexture(GL_TEXTURE_2D, m_pntexDisparity[i]);
//		glBindTexture(GL_TEXTURE_2D, m_pntexDepth[i]);
//		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
////		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, m_ppImages[i]->width, m_ppImages[i]->height, 0, GL_RGBA, GL_FLOAT, NULL);
//		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, m_ppImages[i]->width, m_ppImages[i]->height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
//	}
//
//	glGenTextures(1, &m_ntexDepth);
//	glBindTexture(GL_TEXTURE_2D, m_ntexDepth);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32, m_ppImages[0]->width, m_ppImages[0]->height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
//
//	glGenTextures(1, &m_ntexPredicted);
//	glBindTexture(GL_TEXTURE_2D, m_ntexPredicted);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_ppImages[0]->width, m_ppImages[0]->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
//
//	glGenTextures(1, &m_ntexDissimilarity);
//	glBindTexture(GL_TEXTURE_2D, m_ntexDissimilarity);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, m_ppImages[0]->width, m_ppImages[0]->height, 0, GL_RGBA, GL_FLOAT, NULL);
//
//	m_pDissimilarityBuffer = new float[m_ppImages[0]->height * m_ppImages[0]->width * 4];
//	_ASSERT(m_pDissimilarityBuffer != NULL);
//}
//
//void CGlobalImageMatch::clearTextures()
//{
//	if (m_pntexImages != NULL)
//	{
//		glDeleteTextures(m_nImages, m_pntexImages);
//		delete [] m_pntexImages;
//		m_pntexImages = NULL;
//	}
////	if (m_pntexDisparity != NULL)
////	{
////		glDeleteTextures(m_nImages, m_pntexDisparity);
////		delete [] m_pntexDisparity;
////		m_pntexDisparity = NULL;
////	}
//	if (m_pntexDepth != NULL)
//	{
//		glDeleteTextures(m_nImages, m_pntexDepth);
//		delete [] m_pntexDepth;
//		m_pntexDepth = NULL;
//	}
//
//	if (m_ntexPredicted != 0)
//	{
//		glDeleteTextures(1, &m_ntexPredicted);
//		m_ntexPredicted = 0;
//	}
//	if (m_ntexDissimilarity != 0)
//	{
//		glDeleteTextures(1, &m_ntexDissimilarity);
//		m_ntexDissimilarity = 0;
//	}
//
//	if (m_ntexDepth != 0)
//	{
//		glDeleteTextures(1, &m_ntexDepth);
//		m_ntexDepth = 0;
//	}
//
//	if (m_pDissimilarityBuffer != NULL)
//	{
//		delete [] m_pDissimilarityBuffer;
//		m_pDissimilarityBuffer = NULL;
//	}
//}

//void CGlobalImageMatch::initShaders()
//{
//	GLint success;
//
//	m_nvsProject = loadVertexShader("ProjectRefMatch_3x4.vs");
//	if (m_nvsProject == 0)
//	{
//		printf("Error loading ProjectRefMatch_3x4.vs\n");
//		return;
//	}
//	glCompileShader(m_nvsProject);
//	glGetShaderiv(m_nvsProject, GL_COMPILE_STATUS, &success);
//	if (!success)
//	{
//		shaderError("Error compiling ProjectRefMatch_3x4.vs", m_nvsProject);
//		return;
//	}
//
//#if 0
//	m_nfsDisparity = loadFragmentShader("TemplateDisparityMap.fs");
//	if (m_nfsDisparity == 0)
//	{
//		printf("Error loading TemplateDisparityMap.fs\n");
//		return;
//	}
//	glCompileShader(m_nfsDisparity);
//	glGetShaderiv(m_nfsDisparity, GL_COMPILE_STATUS, &success);
//	if (!success)
//	{
//		shaderError("Error compiling TemplateDisparityMap.fs", m_nfsDisparity);
//		return;
//	}
//#endif
//
//	m_nfsPredict = loadFragmentShader("PredictImage_3x4.fs");
//	if (m_nfsPredict == 0)
//	{
//		printf("Error loading PredictImage_3x4.fs\n");
//		return;
//	}
//	glCompileShader(m_nfsPredict);
//	glGetShaderiv(m_nfsPredict, GL_COMPILE_STATUS, &success);
//	if (!success)
//	{
//		shaderError("Error compiling PredictImage_3x4.fs", m_nfsPredict);
//		return;
//	}
//
//#if 0
//	m_nprogDisparity = glCreateProgram();
//	glAttachShader(m_nprogDisparity, m_nvsProject);
//	glAttachShader(m_nprogDisparity, m_nfsDisparity);
//	glLinkProgram(m_nprogDisparity);
//	glGetProgramiv(m_nprogDisparity, GL_LINK_STATUS, &success);
//	if (!success)
//	{
//		programError("Error linking m_nprogDisparity", m_nprogDisparity);
//		return;
//	}
//	glValidateProgram(m_nprogDisparity);
//	glGetProgramiv(m_nprogDisparity, GL_VALIDATE_STATUS, &success);
//	if (!success)
//	{
//		programError("Error validating m_nprogDisparity", m_nprogDisparity);
//		return;
//	}
//#else
//	m_nprogDepth = glCreateProgram();
//	glAttachShader(m_nprogDepth, m_nvsProject);
//	glLinkProgram(m_nprogDepth);
//	glGetProgramiv(m_nprogDepth, GL_LINK_STATUS, &success);
//	if (!success)
//	{
//		programError("Error linking m_nprogDepth", m_nprogDepth);
//		return;
//	}
//	glValidateProgram(m_nprogDepth);
//	glGetProgramiv(m_nprogDepth, GL_VALIDATE_STATUS, &success);
//	if (!success)
//	{
//		programError("Error validating m_nprogDepth", m_nprogDepth);
//		return;
//	}
//#endif
//
//	m_nprogPredict = glCreateProgram();
//	glAttachShader(m_nprogPredict, m_nvsProject);
//	glAttachShader(m_nprogPredict, m_nfsPredict);
//	glLinkProgram(m_nprogPredict);
//	glGetProgramiv(m_nprogPredict, GL_LINK_STATUS, &success);
//	if (!success)
//	{
//		programError("Error linking m_nprogPredict", m_nprogPredict);
//		return;
//	}
//	glValidateProgram(m_nprogPredict);
//	glGetProgramiv(m_nprogPredict, GL_VALIDATE_STATUS, &success);
//	if (!success)
//	{
//		programError("Error validating m_nprogPredict", m_nprogPredict);
//		return;
//	}
//
//	m_nfsDissim = loadFragmentShader("Dissimilarity_AD_1.fs");
////	m_nfsDissim = loadFragmentShader("Dissimilarity_SAD_1.fs");
//	if (m_nfsDissim == 0)
//	{
//		printf("Error loading Dissimilarity_AD_1.fs\n");
////		printf("Error loading Dissimilarity_SAD_1.fs\n");
//		return;
//	}
//	glCompileShader(m_nfsDissim);
//	glGetShaderiv(m_nfsDissim, GL_COMPILE_STATUS, &success);
//	if (!success)
//	{
//		shaderError("Error compiling Dissimilarity_AD_1.fs", m_nfsDissim);
////		shaderError("Error compiling Dissimilarity_SAD_1.fs", m_nfsDissim);
//		return;
//	}
//
//	m_nprogDissim = glCreateProgram();
//	glAttachShader(m_nprogDissim, m_nfsDissim);
//	glLinkProgram(m_nprogDissim);
//	glGetProgramiv(m_nprogDissim, GL_LINK_STATUS, &success);
//	if (!success)
//	{
//		programError("Error linking m_nprogDissim", m_nprogDissim);
//		return;
//	}
//	glValidateProgram(m_nprogDissim);
//	glGetProgramiv(m_nprogDissim, GL_VALIDATE_STATUS, &success);
//	if (!success)
//	{
//		programError("Error validating m_nprogDissim", m_nprogDissim);
//		return;
//	}
//}
//
//void CGlobalImageMatch::clearShaders()
//{
//#if 0
//	if (m_nprogDisparity != 0)
//	{
//		glDeleteProgram(m_nprogDisparity);
//		m_nprogDisparity = 0;
//	}
//#else
//	if (m_nprogDepth != 0)
//	{
//		glDeleteProgram(m_nprogDepth);
//		m_nprogDepth = 0;
//	}
//#endif
//	if (m_nprogPredict != 0)
//	{
//		glDeleteProgram(m_nprogPredict);
//		m_nprogPredict = 0;
//	}
//	if (m_nprogDissim != 0)
//	{
//		glDeleteProgram(m_nprogDissim);
//		m_nprogDissim = 0;
//	}
//	if (m_nfsDissim != 0)
//	{
//		glDeleteShader(m_nfsDissim);
//		m_nfsDissim = 0;
//	}
//	if (m_nfsPredict != 0)
//	{
//		glDeleteShader(m_nfsPredict);
//		m_nfsPredict = 0;
//	}
//	if (m_nfsDisparity != 0)
//	{
//		glDeleteShader(m_nfsDisparity);
//		m_nfsDisparity = 0;
//	}
//	if (m_nvsProject != 0)
//	{
//		glDeleteShader(m_nvsProject);
//		m_nvsProject = 0;
//	}
//}

//void CGlobalImageMatch::setImages(int nImages, IplImage **ppImages, double **ppProjectionMatrices)
//{
//	clearImages();
//
//	m_nImages = nImages;
//	if (m_nImages > 0)
//	{
//		m_ppImages = new IplImage*[m_nImages];
//		_ASSERT(m_ppImages != NULL);
//		m_ppProjectMat = new double*[m_nImages];
//		_ASSERT(m_ppProjectMat != NULL);
//
//		for (int i = 0; i < m_nImages; i++)
//		{
//			m_ppImages[i] = ppImages[i];
//			m_ppProjectMat[i] = ppProjectionMatrices[i];
//		}
//	}
//}

//void CGlobalImageMatch::clearImages()
//{
//	if (m_ppImages != NULL)
//	{
//		delete [] m_ppImages;
//		m_ppImages = NULL;
//	}
//	if (m_ppProjectMat != NULL)
//	{
//		delete [] m_ppProjectMat;
//		m_ppProjectMat = NULL;
//	}
//	m_nImages = 0;
//}

void CGlobalImageMatch::setSurface(int nVertices, int nIndices, GLenum eType, float *pVertices, GLuint *pIndices, bool bUseBufferObjects)
{
	int i;

//	_ASSERT(pVertices != NULL);
//	_ASSERT(pIndices != NULL);

	m_nVertices = nVertices;
	m_nIndices = nIndices;
	m_ePrimType = eType;
	m_pVertices = pVertices;
	m_pIndices = pIndices;
	m_bUseBufferObjects = bUseBufferObjects;

	if (m_bUseBufferObjects)
	{
		glEnableClientState(GL_VERTEX_ARRAY);

		glGenBuffers(1, &m_nbufVertex);
		glBindBuffer(GL_ARRAY_BUFFER, m_nbufVertex);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_nVertices * sizeof(float), m_pVertices, GL_DYNAMIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glGenBuffers(1, &m_nbufIndex);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_nbufIndex);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_nIndices * sizeof(GLuint), NULL, GL_STATIC_DRAW);

		nIndices = 0;
		for (i = 0; i < m_nIndices; i++)
		{
			if (m_pIndices[i] != gim_indexMaskValue)
			{
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, nIndices * sizeof(GLuint), sizeof(GLuint), m_pIndices + i);
				nIndices++;
			}
		}

		m_nIndices = nIndices;

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		m_bDeleteBufferObjects = true;
	}

}

void CGlobalImageMatch::setSurface(int nVertices, int nIndices, GLenum eType, float* pVertices, float* pNormals, GLuint* pIndices, bool bUseBufferObjects)
{
	int i;

	_ASSERT(pVertices != NULL);
	_ASSERT(pIndices != NULL);
	_ASSERT(pNormals != NULL);

	m_nVertices = nVertices;
	m_nIndices = nIndices;
	m_ePrimType = eType;
	m_pVertices = pVertices;
	m_pNormals = pNormals;
	m_pIndices = pIndices;
	m_bUseBufferObjects = bUseBufferObjects;

	if (m_bUseBufferObjects)
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glGenBuffers(1, &m_nbufVertex);
		glBindBuffer(GL_ARRAY_BUFFER, m_nbufVertex);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_nVertices * sizeof(float), m_pVertices, GL_DYNAMIC_DRAW);//GL_STREAM_COPY);//

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glGenBuffers(1, &m_nbufNormal);
		glBindBuffer(GL_ARRAY_BUFFER, m_nbufNormal);
		glBufferData(GL_ARRAY_BUFFER, 3 * m_nVertices * sizeof(float), m_pNormals, GL_DYNAMIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, 0);

		glGenBuffers(1, &m_nbufIndex);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_nbufIndex);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_nIndices * sizeof(GLuint), NULL, GL_STATIC_DRAW);

		nIndices = 0;
		for (i = 0; i < m_nIndices; i++)
		{
			if (m_pIndices[i] != gim_indexMaskValue)
			{
				glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, nIndices * sizeof(GLuint), sizeof(GLuint), m_pIndices + i);
				nIndices++;
			}
		}

		m_nIndices = nIndices;

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		m_bDeleteBufferObjects = true;
	}

}

//void CGlobalImageMatch::setSurface(int nVertices, int nIndices, GLenum eType, GLuint nbufVertices, GLuint nbufIndices, GLuint nbufNormals)
//{
//	m_nbufVertex = nbufVertices;
//	m_nbufIndex = nbufIndices;
//	m_nbufNormal = nbufNormals;
//	m_bUseBufferObjects = true;
//	glEnableClientState(GL_VERTEX_ARRAY);
//	glEnableClientState(GL_NORMAL_ARRAY);
//	m_bDeleteBufferObjects = false;
//}

void CGlobalImageMatch::clearSurface()
{
	m_pVertices = NULL;
	m_pNormals = NULL;
	m_pIndices = NULL;
	m_nVertices = 0;
	m_nIndices = 0;
	m_bUseBufferObjects = false;
	if (m_bDeleteBufferObjects)
	{
		if (m_nbufVertex != 0)
		{
			glDeleteBuffers(1, &m_nbufVertex);
			m_nbufVertex = 0;
		}
		if (m_nbufNormal != 0)
		{
			glDeleteBuffers(1, &m_nbufNormal);
			m_nbufNormal = 0;
		}
		if (m_nbufIndex != 0)
		{
			glDeleteBuffers(1, &m_nbufIndex);
			m_nbufIndex = 0;
		}
		m_bDeleteBufferObjects = false;
	}
}

//void CGlobalImageMatch::renderSurface()
//{
//	int i, j;
//
//	if (m_bUseBufferObjects)
//	{
//		glBindBuffer(GL_ARRAY_BUFFER, m_nbufVertex);
////		checkGLError("glBindBuffer(GL_ARRAY_BUFFER, g_nbufDispMeshGrid)");
//		glVertexPointer(3, GL_FLOAT, 0, (GLvoid*)0);
////		checkGLError("glVertexPointer(2, GL_FLOAT, 0, (GLvoid*)0)");
//
//		if (m_nbufNormal != 0)
//		{
//			glBindBuffer(GL_ARRAY_BUFFER, m_nbufNormal);
//			glNormalPointer(GL_FLOAT, 0, (GLvoid*)0);
//		}
//
//		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_nbufIndex);
////		checkGLError("glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_nbufDispMeshElem)");
////		glDrawElements(GL_QUADS, m_nIndices / 4, GL_UNSIGNED_INT, (GLvoid*)0);
//		glDrawElements(GL_QUADS, m_nIndices, GL_UNSIGNED_INT, (GLvoid*)0);
////		checkGLError("glDrawElements(GL_QUADS, (g_nDispMeshHeight - 1) * (g_nDispMeshWidth - 1), GL_UNSIGNED_INT, (GLvoid*)0)");
//
//		glBindBuffer(GL_ARRAY_BUFFER, 0);
//		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
//	}
//	else
//	{
//		_ASSERT(m_pIndices != NULL);
//		_ASSERT(m_pVertices != NULL);
//
//		if (m_pNormals != NULL)
//		{
//			glBegin(m_ePrimType);
//
//			for (i = 0; i < m_nIndices; i++)
//			{
//				j = m_pIndices[i];// * 3;
//				if (j != gim_indexMaskValue)
//				{
//					j *= 3;
//					glNormal3fv(m_pNormals + j);//glNormal3f(m_pNormals[j], m_pNormals[j + 1], m_pNormals[j + 2]);
//					glVertex3fv(m_pVertices + j);//glVertex3f(m_pVertices[j], m_pVertices[j + 1], m_pVertices[j + 2]);
//				}
//			}
//
//			glEnd();
//		}
//		else
//		{
//			glBegin(m_ePrimType);
//
//			for (i = 0; i < m_nIndices; i++)
//			{
//				j = m_pIndices[i];// * 3;
//				if (j != gim_indexMaskValue)
//				{
//					j *= 3;
//					glVertex3f(m_pVertices[j], m_pVertices[j + 1], m_pVertices[j + 2]);
//				}
//			}
//
//			glEnd();
//		}
//	}
//}

//void CGlobalImageMatch::setUniformParameters()
//{
//	GLint ntemp;
//	GLuint nprog;
//	float vecTemp[4], matTemp[16];
//
//	glGetIntegerv(GL_CURRENT_PROGRAM, &ntemp);
//	nprog = (GLuint)ntemp;
//
//	int iMatRef = glGetUniformLocation(nprog, "matRef");
//	int iGridRef = glGetUniformLocation(nprog, "gridRef");
//	int iMatMatch = glGetUniformLocation(nprog, "matMatch");
//	int iGridMatch = glGetUniformLocation(nprog, "gridMatch");
//	int iBaseline = glGetUniformLocation(nprog, "baseline");
//	int iFocal = glGetUniformLocation(nprog, "focal");
//
//	vecTemp[2] = (float)m_ppImages[m_iReference]->width;
//	vecTemp[3] = (float)m_ppImages[m_iReference]->height;
//	vecTemp[0] = 1.0 / vecTemp[2];
//	vecTemp[1] = 1.0 / vecTemp[3];
//	glUniform4fv(iGridRef, 1, vecTemp);
//
//#if 1
//	matTemp[0] = (float)m_ppProjectMat[m_iReference][0];
//	matTemp[1] = (float)m_ppProjectMat[m_iReference][4];
//	matTemp[2] = (float)m_ppProjectMat[m_iReference][8];
//	matTemp[3] = (float)m_ppProjectMat[m_iReference][1];
//	matTemp[4] = (float)m_ppProjectMat[m_iReference][5];
//	matTemp[5] = (float)m_ppProjectMat[m_iReference][9];
//	matTemp[6] = (float)m_ppProjectMat[m_iReference][2];
//	matTemp[7] = (float)m_ppProjectMat[m_iReference][6];
//	matTemp[8] = (float)m_ppProjectMat[m_iReference][10];
//	matTemp[9] = (float)m_ppProjectMat[m_iReference][3];
//	matTemp[10] = (float)m_ppProjectMat[m_iReference][7];
//	matTemp[11] = (float)m_ppProjectMat[m_iReference][11];
//	glUniformMatrix4x3fv(iMatRef, 1, GL_FALSE, matTemp);
////	glUniformMatrix4x3fv(iMatRef, 1, GL_TRUE, matTemp);
//#elif 0
//	matTemp[0] = (float)m_ppProjectMat[m_iReference][0];
//	matTemp[1] = (float)m_ppProjectMat[m_iReference][1];
//	matTemp[2] = (float)m_ppProjectMat[m_iReference][2];
//	matTemp[3] = (float)m_ppProjectMat[m_iReference][3];
//	matTemp[4] = (float)m_ppProjectMat[m_iReference][4];
//	matTemp[5] = (float)m_ppProjectMat[m_iReference][5];
//	matTemp[6] = (float)m_ppProjectMat[m_iReference][6];
//	matTemp[7] = (float)m_ppProjectMat[m_iReference][7];
//	matTemp[8] = (float)m_ppProjectMat[m_iReference][8];
//	matTemp[9] = (float)m_ppProjectMat[m_iReference][9];
//	matTemp[10] = (float)m_ppProjectMat[m_iReference][10];
//	matTemp[11] = (float)m_ppProjectMat[m_iReference][11];
//	glUniformMatrix4x3fv(iMatRef, 1, GL_TRUE, matTemp);
//#else
//	matTemp[0] = (float)m_ppProjectMat[m_iReference][0];
//	matTemp[1] = (float)m_ppProjectMat[m_iReference][4];
//	matTemp[2] = (float)m_ppProjectMat[m_iReference][8];
//	matTemp[3] = 0.f;
//	matTemp[4] = (float)m_ppProjectMat[m_iReference][1];
//	matTemp[5] = (float)m_ppProjectMat[m_iReference][5];
//	matTemp[6] = (float)m_ppProjectMat[m_iReference][9];
//	matTemp[7] = 0.f;
//	matTemp[8] = (float)m_ppProjectMat[m_iReference][2];
//	matTemp[9] = (float)m_ppProjectMat[m_iReference][6];
//	matTemp[10] = (float)m_ppProjectMat[m_iReference][10];
//	matTemp[11] = 0.f;
//	matTemp[12] = (float)m_ppProjectMat[m_iReference][3];
//	matTemp[13] = (float)m_ppProjectMat[m_iReference][7];
//	matTemp[14] = (float)m_ppProjectMat[m_iReference][11];
//	matTemp[15] = 1.f;
//	glUniformMatrix4fv(iMatRef, 1, GL_FALSE, matTemp);
////	glUniformMatrix4x3fv(iMatRef, 1, GL_FALSE, matTemp);
//#endif
//
//	vecTemp[2] = (float)m_ppImages[m_iMatching]->width;
//	vecTemp[3] = (float)m_ppImages[m_iMatching]->height;
//	vecTemp[0] = 1.0 / vecTemp[2];
//	vecTemp[1] = 1.0 / vecTemp[3];
//	glUniform4fv(iGridMatch, 1, vecTemp);
//
//	matTemp[0] = (float)m_ppProjectMat[m_iMatching][0];
//	matTemp[1] = (float)m_ppProjectMat[m_iMatching][4];
//	matTemp[2] = (float)m_ppProjectMat[m_iMatching][8];
//	matTemp[3] = (float)m_ppProjectMat[m_iMatching][1];
//	matTemp[4] = (float)m_ppProjectMat[m_iMatching][5];
//	matTemp[5] = (float)m_ppProjectMat[m_iMatching][9];
//	matTemp[6] = (float)m_ppProjectMat[m_iMatching][2];
//	matTemp[7] = (float)m_ppProjectMat[m_iMatching][6];
//	matTemp[8] = (float)m_ppProjectMat[m_iMatching][10];
//	matTemp[9] = (float)m_ppProjectMat[m_iMatching][3];
//	matTemp[10] = (float)m_ppProjectMat[m_iMatching][7];
//	matTemp[11] = (float)m_ppProjectMat[m_iMatching][11];
//	glUniformMatrix4x3fv(iMatMatch, 1, false, matTemp);
//
//	glUniform1f(iFocal, m_disparityFocal);
//	glUniform1f(iBaseline, m_disparityBaseline);
//}

//void CGlobalImageMatch::computeDisparity()
//void CGlobalImageMatch::computeDepth()
//{
//	//set program to use
////	glUseProgram(m_nprogDisparity);
//	glUseProgram(m_nprogDepth);
//
//	setUniformParameters();
//
//	//attach disparity texture to framebuffer
////	fbAttachTexture(m_frameBuffer, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_pntexDisparity[m_iReference]);
////	fbAttachTexture(m_frameBuffer, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, m_ntexDepth);
//	fbAttachTexture(m_frameBuffer, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, m_pntexDepth[m_iReference]);
////	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
//	glDrawBuffer(GL_NONE);
//	glViewport(0, 0, m_ppImages[m_iReference]->width, m_ppImages[m_iReference]->height);
////	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	glClear(GL_DEPTH_BUFFER_BIT);
//
//	//process geometry
//	renderSurface();
//
//	//clean-up
//	fbDetachAll(m_frameBuffer);
//}

//void CGlobalImageMatch::predictImage()
//{
//	GLint iloc;
//
//	//set program to use
//	glUseProgram(m_nprogPredict);
//
//	setUniformParameters();
//
//	//attach predicted image texture to framebuffer
//	fbAttachTexture(m_frameBuffer, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_ntexPredicted);
//	fbAttachTexture(m_frameBuffer, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, m_ntexDepth);
//	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
//	glViewport(0, 0, m_ppImages[0]->width, m_ppImages[0]->height);
//	glClearColor(0.0, 0.0, 0.0, 0.0);
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//
//	//bind matching image and matching depth/disparity
//	glActiveTexture(GL_TEXTURE1);
////	glBindTexture(GL_TEXTURE_2D, m_pntexDisparity[m_iMatching]);
//	glBindTexture(GL_TEXTURE_2D, m_pntexDepth[m_iMatching]);
//	glActiveTexture(GL_TEXTURE0);
//	glBindTexture(GL_TEXTURE_2D, m_pntexImages[m_iMatching]);
//
//	iloc = glGetUniformLocation(m_nprogPredict, "texMatch");
//	glUniform1i(iloc, 0);
////	iloc = glGetUniformLocation(m_nprogPredict, "texDisparity");
//	iloc = glGetUniformLocation(m_nprogPredict, "texDepth");
//	glUniform1i(iloc, 1);
//
//	//process geometry
//	renderSurface();
//
//	//clean-up
//	fbDetachAll(m_frameBuffer);
//	glActiveTexture(GL_TEXTURE1);
//	glBindTexture(GL_TEXTURE_2D, 0);
//	glActiveTexture(GL_TEXTURE0);
//	glBindTexture(GL_TEXTURE_2D, 0);
//}

//void CGlobalImageMatch::computeDissimilarity()
//{
//	GLint iloc;
//
//	//set program to use
//	glUseProgram(m_nprogDissim);
//
//	setUniformParameters();
//
//	//attach dissimilarity texture to framebuffer
//	fbAttachTexture(m_frameBuffer, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_ntexDissimilarity);
//	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
//	glViewport(0, 0, m_ppImages[0]->width, m_ppImages[0]->height);
//	glClear(GL_COLOR_BUFFER_BIT);
//
//	//bind reference image and predicted image
//	glActiveTexture(GL_TEXTURE1);
//	glBindTexture(GL_TEXTURE_2D, m_ntexPredicted);
//	glActiveTexture(GL_TEXTURE0);
//	glBindTexture(GL_TEXTURE_2D, m_pntexImages[m_iReference]);
//
//	iloc = glGetUniformLocation(m_nprogDissim, "texReference");
//	glUniform1i(iloc, 0);
//	iloc = glGetUniformLocation(m_nprogDissim, "texPredicted");
//	glUniform1i(iloc, 1);
//
//	//process pixels/fragments via screen-aligned geometry
//	renderScreenAligned();
//
//	//clean-up
//	fbDetachAll(m_frameBuffer);
//	glActiveTexture(GL_TEXTURE1);
//	glBindTexture(GL_TEXTURE_2D, 0);
//	glActiveTexture(GL_TEXTURE0);
//	glBindTexture(GL_TEXTURE_2D, 0);
//}

//double CGlobalImageMatch::readBackSumDissimilarity()
//{
//	int i, nElements;
//	double sum = 0.0;
//
//	glBindTexture(GL_TEXTURE_2D, m_ntexDissimilarity);
//	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, m_pDissimilarityBuffer);
//	glBindTexture(GL_TEXTURE_2D, 0);
//
//	nElements = m_ppImages[0]->height * m_ppImages[0]->width * 4;
//	for (i = 0; i < nElements; i+=4)
//	{
//		sum += (double)m_pDissimilarityBuffer[i];
//		sum += (double)m_pDissimilarityBuffer[i + 1];
//		sum += (double)m_pDissimilarityBuffer[i + 2];
//		if (m_pDissimilarityBuffer[i + 3] > 0.5)
//			m_nPixelCount++;
//	}
//
//	return sum;
//}

//double CGlobalImageMatch::computeGlobalMatchingCost()
//{
//	int i, j;
//
//	bool bReadbackSum = false;
//
//	double matchingCost = 0.0;
//
//	if (m_bUseBufferObjects && m_pVertices != NULL)
//	{
//		glBindBuffer(GL_ARRAY_BUFFER, m_nbufVertex);
//		glBufferSubData(GL_ARRAY_BUFFER, 0, m_nVertices * 3 * sizeof(float), m_pVertices);
//	}
//
//	fbBind(m_frameBuffer);
//
//	for (i = 0; i < m_nImages; i++)
//	{
//		setReferenceImage(i);
////		computeDisparity();
//		computeDepth();
//	}
//
//	m_nPixelCount = 0;
//
//	for (i = 0; i < m_nImages; i++)
//	{
//		setReferenceImage(i);
//
//		for (j = 0; j < i; j++)
//		{
//			setMatchingImage(j);
//			predictImage();
//			if (bReadbackSum)
//				matchingCost += readBackSumDissimilarity();
//			computeDissimilarity();
//			bReadbackSum = true;
//
//#if 0
//			//DEBUG
//			fbUnbind();
////			glBindTexture(GL_TEXTURE_2D, m_ntexPredicted);
//			glBindTexture(GL_TEXTURE_2D, m_ntexDissimilarity);
////			glBindTexture(GL_TEXTURE_2D, m_pntexDisparity[j]);
//			glUseProgram(0);
//			renderScreenAligned();
//			glutSwapBuffers();
//			fbBind();
////			Sleep(2000);
//			//END DEBUG
//#endif
//
////			matchingCost += readBackSumDissimilarity();
//		}
//
//		for (j = i + 1; j < m_nImages; j++)
//		{
//			setMatchingImage(j);
//			predictImage();
//			if (bReadbackSum)
//				matchingCost += readBackSumDissimilarity();
//			computeDissimilarity();
//			bReadbackSum = true;
//
//#if 0
//			//DEBUG
//			fbUnbind();
////			glBindTexture(GL_TEXTURE_2D, m_ntexPredicted);
//			glBindTexture(GL_TEXTURE_2D, m_ntexDissimilarity);
////			glBindTexture(GL_TEXTURE_2D, m_pntexDisparity[j]);
//			glUseProgram(0);
//			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//			renderScreenAligned();
//			glutSwapBuffers();
//			fbBind();
////			Sleep(2000);
//			//END DEBUG
//#endif
//
////			matchingCost += readBackSumDissimilarity();
//		}
//	}
//
//	matchingCost += readBackSumDissimilarity();
//	matchingCost /= (double)m_nPixelCount;
//
//	fbUnbind();
//
//	return matchingCost;
//}




