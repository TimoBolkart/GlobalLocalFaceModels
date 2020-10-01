#include "StereoFace.h"
#include "BSplineGridWavelet.h"
#include "BSGWFactory.h"
#include "NearestNeighborAssistant.h"
#include "WaveletShapeSampler.h"
#include "NearestNeighborEnergy.h"
#include "FaceData.h"
#include "GlobalImageMatch.h"

#include <fstream>
#include <vector>

//Model specific constants used for the projection plane
#define SF_EMPERICAL_VIEW_X                -0.6704
#define SF_EMPERICAL_VIEW_Y                -0.7420
#define SF_EMPERICAL_VIEW_Z                0.0

//constants for indices of some landmarks on mean face of neutral expression
//for high-res scaled align version of parameterized BUFE3D database
//used to align template with 2D grid for wavelet decomposition
#define SF_BUFE3D_MEAN_NE_NOSETIP			4303
#define SF_BUFE3D_MEAN_NE_NOSEBRIDGE		3692
#define SF_BUFE3D_MEAN_NE_LEFTEYEOC			4736
#define SF_BUFE3D_MEAN_NE_RIGHTEYEOC		2358

//Width and height of the initial base mesh used for subdivision
#define SF_MODEL_BASE_MESH_WIDTH			5
#define SF_MODEL_BASE_MESH_HEIGHT		7
//Number of subdivision levels
#define SF_MODEL_LEVELS						6

//Defines the numbe
#define NUMBER_SAMPLES						50
#define SAMPLE_REGION						1.0

//Define the number of fitting refinement steps (needs to be > 0 and <= SF_MODEL_LEVELS)
#define SF_REFINE_GEOMETRY_LEVEL			SF_MODEL_LEVELS

namespace clapack
{
	extern "C"
	{
		#include "blaswrap.h"
		#include "f2c.h"
		int dgemm_(char *transa, char *transb, integer *m, integer *
			n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
			doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
			integer *ldc);
		int dgels_(char *trans, integer *m, integer *n, integer *
			nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
			doublereal *work, integer *lwork, integer *info);
		int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
			doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
			ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
			integer *info);
	}
}

//Help functions:

//Define plane for 2D sterographic projection and compute the projection.
GLuint initTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf> *& wavelet, GLfloat *& pGeomMapMask, double *& pTemplate2DPositions,
							   double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop,
							   GLuint *& pGeomMapIndices, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double *& matTemplateStereoMap, 
							   IplImage *& pTemplate3DMap, abutil::C3Vectorf *& pGeomMapNormals);

IplImage* generateTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf> * wavelet, GLfloat * pGeomMapMask, double * pFaceVertices, GLuint & nprogTemplateGeometryGen, 
									  GLuint & nfb, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double * matTemplateStereoMap,
									  double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop, double * pTemplate2DPositions, 
									  IplImage * pTemplate3DMap, GLuint * pGeomMapIndices);

void initGeometryMapIndices(int nWidth, int nHeight, GLuint * pIndices);
void maskGeometryMapIndices(int nWidth, int nHeight, GLfloat * pMask, GLuint * pIndices);
void initDisplay();

//Output geometry file in OFF file format.
void saveGeometryMapAsOFF(char* szFilename, IplImage* pTemplate3DMap, CBSplineGridWavelet<C3Vectorf>* pWavelet, GLuint* pGeomMapIndices);

void initGeometryMapColors(CBSplineGridWavelet<C3Vectorf>* pStdDevShape, abutil::C3Vectorf* pStdDevRecon, abutil::C3Vectorf *& pGeomMapColors);
void initStereoPoints(std::vector<SStereoPatch> patches, float *& pStereoPoints, int & nStereoPoints, CNearestNeighborAssistant *& NNA);
void updateGeometryMap(bool bUpdateNormals, IplImage *& pTemplate3DMap, CBSplineGridWavelet<abutil::C3Vectorf> *& wavelet, GLuint & nbufGeomMap, CGlobalImageMatch & imageMatcher,
					   abutil::C3Vectorf *& pGeomMapNormals, GLuint & nbufGeomMapNormals, GLuint & ntexTemplateGeometryMap);
void initWaveletShapeSampler(CWaveletShapeSampler & WSS, CBSplineGridWavelet<abutil::C3Vectorf> *& pMeanShape, CBSplineGridWavelet<abutil::C3Vectorf> *& wavelet, IplImage *& pTemplate3DMap, 
							 GLuint & nbufGeomMap, CGlobalImageMatch & imageMatcher, abutil::C3Vectorf *& pGeomMapNormals, GLuint & nbufGeomMapNormals, GLuint & ntexTemplateGeometryMap,
							 CWaveletShapePrior *& pStatPrior, float *& pmatReductTranspose, CNearestNeighborEnergy *& pStereoDataNN, float *& pStereoPoints, int& nStereoPoints,
							 double *& matAlignModelToData, double *& matAlignDataToModel, CBSplineGridWavelet<C3Vectorf> *& pStdDevShape, GLfloat *& pGeomMapMask, CNearestNeighborAssistant *& NNA);

//Read landmarks
bool readLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded);

//Compute transformation to align observed points with the model
void computeLSAlignmentModelToData(int nPoints, double* pModelPoints, double* pObservedPoints, double* pmatTransform, double* pmatTransformInv);
bool computeLSAlignmentModelToData(const std::vector<double>& modelLandmarks, const std::vector<bool>& modelLandmarksLoaded, const std::vector<double>& dataLandmarks, const std::vector<bool>& dataLandmarksLoaded
											, double* pmatTransform, double* pmatTransformInv);
void computeLSAlignmentModelToDataFromNN(double *& pmatTransform, double *& pmatTransformInv, IplImage *& pTemplate3DMap, CNearestNeighborAssistant * NNA, GLfloat*& pGeomMapMask,
										 CBSplineGridWavelet<abutil::C3Vectorf> * wavelet, float *& pStereoPoints);

///////////////////////////////////////////////////////////////////////////////
//main - program starts here
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
	int j, k, windowId, nStereoPoints, nStep;
	double templateMapLeft, templateMapRight, templateMapBottom, templateMapTop, xVal, yVal, zVal;
	double * pTemplate2DPositions, * matTemplateStereoMap, * matAlignModelToData, * matAlignDataToModel;
	float * pmatReductTranspose, * pStereoPoints, * pVertData;
	char * outputFileName, * vrmlFile, * landmarksData, * landmarksModel;
	char tempHelpFile [1000];
	GLfloat * pGeomMapMask;
	GLuint * pGeomMapIndices;
	GLuint nfb, ntexTemplateGeometryMap, ntexTemplateDepthMap, nprogTemplateGeometryGen, nbufGeomMap, nbufGeomMapNormals;
	CBSplineGridWavelet<C3Vectorf> * wavelet, * pMeanShape, * pStdDevShape, * pEigenShape;
	abutil::C3Vectorf * pStdDevRecon, * pGeomMapColors, * pGeomMapNormals;
	IplImage * pTemplate3DMap, * bestTemplate3DMap;
	CNearestNeighborAssistant * NNA;
	CWaveletShapePrior * pStatPrior = NULL;
	CNearestNeighborEnergy * pStereoDataNN = NULL;

	//Create a Window (needed for parameterization):
	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize( 1024, 768 );
	windowId = glutCreateWindow( "StereoFace" );

	//Setup:
	initDisplay();
	//initialize the framebuffer
	fbInit(nfb);

	if(argc != 7)
	{
		printf("Wrong number of command line arguments.\n");
		return 1;
	}

	printf("Initializing face data...\n");
	int index(0);

	const char* cstrModelPath = argv[++index];

	const char* cstrFileName = argv[++index];
	if(!CFaceData::init(cstrFileName))
	{
		printf("Unable to load template mesh %s.\n", cstrFileName);
		return 1;
	}
	
	printf("Initializing template geometry resampling...\n");
	wavelet = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, SF_MODEL_BASE_MESH_WIDTH, SF_MODEL_BASE_MESH_HEIGHT, SF_MODEL_LEVELS);
	_ASSERT(wavelet != NULL);
	
	matTemplateStereoMap = new double[16];

	pTemplate3DMap = NULL;
	pGeomMapIndices = NULL;
	pTemplate2DPositions = NULL;
	pGeomMapMask = NULL;

	nprogTemplateGeometryGen = initTemplateGeometryMap(wavelet, pGeomMapMask, pTemplate2DPositions, 
		templateMapLeft, templateMapRight, templateMapBottom, templateMapTop, pGeomMapIndices, ntexTemplateGeometryMap, ntexTemplateDepthMap, 
		matTemplateStereoMap, pTemplate3DMap, pGeomMapNormals);

	if(nprogTemplateGeometryGen == -1)
	{
		printf("Problem in initTemplateGeometryMap\n");
		return 1;
	}

	CvSize mapSize;
	mapSize.width = wavelet->getFullResWidth();
	mapSize.height = wavelet->getFullResHeight();
	bestTemplate3DMap = cvCreateImage(mapSize, IPL_DEPTH_32F, 3);
	_ASSERT(bestTemplate3DMap != NULL);

	//Load model landmarks
	if(argc > ++index)
		landmarksModel = argv[index];
	else
	{
		printf("Missing command line argument on landmarks of model.\n");
		return 1;
	}

	std::vector<double> modelLandmarks;
	std::vector<bool> modelLandmarksLoaded;
	if(!readLandmarks(std::string(landmarksModel), modelLandmarks, modelLandmarksLoaded))
	{
		printf("Problems reading model landmarks %s.\n", landmarksModel);
		return 1;
	}

	//Read a 3D point cloud (possibly generated using stereo data)
	//Right now, this method supports only vrml meshes, but in the future, we may want to allow reading a point cloud here
	if(argc > ++index)
		vrmlFile = argv[index];
	else
	{
		printf("Missing command line argument on vrml file of input data.\n");
		return 1;
	}

	std::vector<SStereoPatch> patches;
	CFaceData::readOFFMeshAsPatch(vrmlFile, patches, 1.0);	

	//Read 3D landmarks (a) on the input point cloud
	if(argc > ++index)
		landmarksData = argv[index];
	else
	{
		printf("Missing command line argument on landmarks of input data.\n");
		return 1;
	}

	std::vector<double> dataLandmarks;
	std::vector<bool> dataLandmarksLoaded;
	if(!readLandmarks(std::string(landmarksData), dataLandmarks, dataLandmarksLoaded))
	{
		printf("Problems reading target landmarks %s.\n", landmarksData);
		return 1;
	}

 	matAlignModelToData = new double[16];
	matAlignDataToModel = new double[16];

	if(!computeLSAlignmentModelToData(modelLandmarks, modelLandmarksLoaded, dataLandmarks, dataLandmarksLoaded, matAlignModelToData, matAlignDataToModel))
	{
		printf("Problems computing rigid alignment");
		return 1;
	}

	//Read the output filename:
	if(argc > ++index)
		outputFileName = argv[index];
	else
	{
		printf("Missing command line argument on how many PCA spaces are used.\n");
		return 1;
	}

	//Fit the learned shape space to the data.
	//Initialize stereo point cloud and nearest neighbor assistant
	NNA = new CNearestNeighborAssistant();
	pStereoPoints = NULL;
	nStereoPoints = 0;
	initStereoPoints(patches, pStereoPoints, nStereoPoints, NNA);

	printf("Generating template geometry map...\n");
	generateTemplateGeometryMap(wavelet, pGeomMapMask, CFaceData::m_pAvgVertices, nprogTemplateGeometryGen, nfb, ntexTemplateGeometryMap, ntexTemplateDepthMap, matTemplateStereoMap,
									  templateMapLeft, templateMapRight, templateMapBottom, templateMapTop, pTemplate2DPositions, pTemplate3DMap, pGeomMapIndices);
	
	CGlobalImageMatch imageMatcher;
	//imageMatcher.setFrameBuffer(nfb);
	uchar* pData;
	cvGetRawData(pTemplate3DMap, &pData);
	imageMatcher.setSurface(pTemplate3DMap->height * pTemplate3DMap->width, (pTemplate3DMap->height - 1) * (pTemplate3DMap->width - 1) * 4, GL_QUADS, 
										(float*)pData, (float*)pGeomMapNormals, pGeomMapIndices, true);

	updateGeometryMap(true, pTemplate3DMap, wavelet, nbufGeomMap, imageMatcher, pGeomMapNormals,  nbufGeomMapNormals, ntexTemplateGeometryMap);

	pGeomMapColors = NULL;

	//Do a finer alignment (one step of ICP): 
	computeLSAlignmentModelToDataFromNN(matAlignModelToData, matAlignDataToModel, pTemplate3DMap, NNA, pGeomMapMask, wavelet, pStereoPoints);

	//Initialize for nearest neighbor computation:
	NNAReal* pModel = new NNAReal[wavelet->getNumCoefficients()*3];
	_ASSERT(pModel != NULL);
	NNAReal* pDistances = new NNAReal[wavelet->getNumCoefficients()];
	_ASSERT(pDistances != NULL);
	NNAIndex* pIndices = new NNAIndex[wavelet->getNumCoefficients()];
	_ASSERT(pIndices != NULL);

	//Start a new wavelet shape sampler for each learned shape space:
	CWaveletShapeSampler WSS;

	//Read the learned shape space:
	sprintf(tempHelpFile, "%s\\mean_shape_%d.bsgw", cstrModelPath, 0);
	if (CBSGWFactory::loadBSGW(tempHelpFile, (void**)&pMeanShape) != eBSGWFloat3)
	{
		printf("wrong data type in %s\n", tempHelpFile);
		return 1;
	}
	if (pMeanShape == NULL)
	{
		printf("unable to load %s\n", tempHelpFile);
		return 1;
	}

	sprintf(tempHelpFile, "%s\\std_shape_%d.bsgw", cstrModelPath, 0);
	if (CBSGWFactory::loadBSGW(tempHelpFile, (void**)&pStdDevShape) != eBSGWFloat3)
	{
		printf("wrong data type in %s\n", tempHelpFile);
		return 1;
	}
	if (pStdDevShape == NULL)
	{
		printf("unable to load %s\n", tempHelpFile);
		return 1;
	}
	else
	{
		pStdDevRecon = new abutil::C3Vectorf[pStdDevShape->getNumCoefficients()];
		_ASSERT(pStdDevRecon != NULL);
		initGeometryMapColors(pStdDevShape, pStdDevRecon, pGeomMapColors);
	}

	pEigenShape = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, pMeanShape->getBaseResWidth(), pMeanShape->getBaseResHeight(), pMeanShape->getNumLevels());
	if (pEigenShape == NULL)
	{
		printf("unable to allocate pEigenShape\n");
	}

	pmatReductTranspose = new float[pMeanShape->getNumCoefficients() * 9];
	_ASSERT(pmatReductTranspose != NULL);

	sprintf(tempHelpFile, "%s\\reduct_tran_shape_%d.dat", cstrModelPath, 0);
	FILE* pfReduct = fopen(tempHelpFile, "rb");
	if(pfReduct == NULL)
	{
		printf("unable to load %s\n", tempHelpFile);
		return 1;
	}
	else
	{
		fread(pmatReductTranspose, sizeof(float), pMeanShape->getNumCoefficients() * 9, pfReduct);
		fclose(pfReduct);
	}

	pStatPrior = NULL;
	pStereoDataNN = NULL;

	for (j = 0; j < SF_REFINE_GEOMETRY_LEVEL; j++)
	{
		initWaveletShapeSampler(WSS, pMeanShape, wavelet, pTemplate3DMap, nbufGeomMap, imageMatcher, pGeomMapNormals, nbufGeomMapNormals, ntexTemplateGeometryMap,
							pStatPrior, pmatReductTranspose, pStereoDataNN, pStereoPoints, nStereoPoints, matAlignModelToData, matAlignDataToModel,
							pStdDevShape, pGeomMapMask, NNA);

		printf("optimizing model level %d of shape space %d ...\n", WSS.getLevel(), 0);
		    
		WSS.optimizeActiveWaveletModel();

		WSS.copyIntermediateWaveletModel(wavelet);
		updateGeometryMap(true, pTemplate3DMap, wavelet, nbufGeomMap, imageMatcher, pGeomMapNormals,  nbufGeomMapNormals, ntexTemplateGeometryMap);

		//Transfer the result into the coordinate system of the data:
		cvGetRawData(pTemplate3DMap, &pData, &nStep);
		pVertData = (float*)pData;
		for(k = 0; k < wavelet->getNumCoefficients(); k++)
		{
			xVal = matAlignModelToData[0] * pVertData[k*3] + matAlignModelToData[1] * pVertData[k*3+1] + matAlignModelToData[2] * pVertData[k*3+2] + matAlignModelToData[3];
			yVal = matAlignModelToData[4] * pVertData[k*3] + matAlignModelToData[5] * pVertData[k*3+1] + matAlignModelToData[6] * pVertData[k*3+2] + matAlignModelToData[7];
			zVal = matAlignModelToData[8] * pVertData[k*3] + matAlignModelToData[9] * pVertData[k*3+1] + matAlignModelToData[10] * pVertData[k*3+2] + matAlignModelToData[11];
			pVertData[k*3] = xVal;
			pVertData[k*3+1] = yVal;
			pVertData[k*3+2] = zVal;
		}

		sprintf(tempHelpFile, "%s_level_%02i_%02i.off", outputFileName, 0, WSS.getLevel() - 1);
		saveGeometryMapAsOFF(tempHelpFile, pTemplate3DMap, wavelet, pGeomMapIndices);

		WSS.reportProfiling(stdout);
	}

	//Compute the cost of the solution
	cvGetRawData(pTemplate3DMap, &pData);
	C3Vectorf * pverts = (C3Vectorf*)pData;

	for (j = 0; j < wavelet->getNumCoefficients(); j++)
	{
		pModel[3*j]		= pverts[j].x;
		pModel[3*j + 1]	= pverts[j].y;
		pModel[3*j + 2]	= pverts[j].z;
		pDistances[j] = 0.0;
	}

	NNA->setQueryPoints(wavelet->getNumCoefficients(), pModel, pDistances, pIndices);
	NNA->compute();

	//Delete stuff
	delete pStatPrior;
	delete pStereoDataNN;
	delete pMeanShape;
	delete pStdDevShape;
	delete pEigenShape;
	delete [] pmatReductTranspose;

	sprintf(tempHelpFile, "%s.off", outputFileName);
	saveGeometryMapAsOFF(tempHelpFile, bestTemplate3DMap, wavelet, pGeomMapIndices);

	//Free space:
	glutDestroyWindow(windowId);
	delete [] pStdDevRecon;
	delete [] matTemplateStereoMap;
	delete [] pGeomMapColors;
	delete NNA;
	delete [] matAlignModelToData;
	delete [] matAlignDataToModel;
	cvReleaseImage(&bestTemplate3DMap);
	delete [] pModel;
	delete [] pDistances;
	delete [] pIndices;
	delete wavelet;

	//Free stuff allocated in initTemplateGeometryMap
	cvReleaseImage(&pTemplate3DMap);
	delete [] pGeomMapMask;
	delete [] pTemplate2DPositions;
	delete [] pGeomMapIndices;

	//Free stuff allocated in initStereoPoints
	delete [] pStereoPoints;

	return 0;
}

GLuint initTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf> *& wavelet, GLfloat *& pGeomMapMask, double *& pTemplate2DPositions,
							   double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop,
							   GLuint *& pGeomMapIndices, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double *& matTemplateStereoMap, 
							   IplImage *& pTemplate3DMap, abutil::C3Vectorf *& pGeomMapNormals)
{
	GLint success;
	CvSize mapSize;
	double empNorm[3];

	//init image
	mapSize.width = wavelet->getFullResWidth();
	mapSize.height = wavelet->getFullResHeight();

	pTemplate3DMap = cvCreateImage(mapSize, IPL_DEPTH_32F, 3);
	_ASSERT(pTemplate3DMap != NULL);

	pGeomMapNormals = new C3Vectorf[mapSize.height * mapSize.width];
	_ASSERT(pGeomMapNormals != NULL);

	pGeomMapMask = new GLfloat[mapSize.height * mapSize.width]; 
	_ASSERT(pGeomMapMask != NULL);

	//init 2D mapping
	pTemplate2DPositions = new double[CFaceData::m_nVertices * 2];
	_ASSERT(pTemplate2DPositions != NULL);

	//select emperically derived plane for stereographic projection
	empNorm[0] = SF_EMPERICAL_VIEW_X;
	empNorm[1] = SF_EMPERICAL_VIEW_Y;
	empNorm[2] = SF_EMPERICAL_VIEW_Z;
	normalize(empNorm);

	//use landmark indices to fix horizontal and vertical grid directions
	int landmarkIndices[4] = { SF_BUFE3D_MEAN_NE_NOSETIP, SF_BUFE3D_MEAN_NE_NOSEBRIDGE, SF_BUFE3D_MEAN_NE_LEFTEYEOC, SF_BUFE3D_MEAN_NE_RIGHTEYEOC };
	CFaceData::empericalPlaneFromLandmarkIndices(empNorm, landmarkIndices, 4, matTemplateStereoMap);
	matTemplateStereoMap[12] = 0.0;
	matTemplateStereoMap[13] = 0.0;
	matTemplateStereoMap[14] = 0.0;
	matTemplateStereoMap[15] = 1.0;

	//compute 2D stereographic projections of vertices, and extents of projection
	CFaceData::stereographicProject(CFaceData::m_pAvgVertices, CFaceData::m_nVertices, matTemplateStereoMap, pTemplate2DPositions, &templateMapLeft, &templateMapBottom, &templateMapRight, &templateMapTop);
	
	//init indices for display
	int nQuadsX = pTemplate3DMap->width - 1;
	int nQuadsY = pTemplate3DMap->height - 1;
	int nQuads = nQuadsX * nQuadsY;
	int nIndices = nQuads * 4;
	pGeomMapIndices = new GLuint[nIndices];
	_ASSERT(pGeomMapIndices != NULL);

	initGeometryMapIndices(pTemplate3DMap->width, pTemplate3DMap->height, pGeomMapIndices);

	//init texture(s)
	ntexTemplateGeometryMap = 0;
	glGenTextures(1, &ntexTemplateGeometryMap);
	glBindTexture(GL_TEXTURE_2D, ntexTemplateGeometryMap);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, mapSize.width, mapSize.height, 0, GL_RGBA, GL_FLOAT, NULL);

	glBindTexture(GL_TEXTURE_2D, 0);

	ntexTemplateDepthMap = 0;
	glGenTextures(1, &ntexTemplateDepthMap);
	glBindTexture(GL_TEXTURE_2D, ntexTemplateDepthMap);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32_ARB, mapSize.width, mapSize.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	checkGLError("glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT32_ARB, mapSize.width, mapSize.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL)");

	glBindTexture(GL_TEXTURE_2D, 0);

	//init shaders
	GLuint nvsTemplateGeometryGen = loadVertexShader("TemplateGeometryProject.vs");
	glCompileShader(nvsTemplateGeometryGen);
	glGetShaderiv(nvsTemplateGeometryGen, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		shaderError("Error compiling TemplateGeometryProject.vs", nvsTemplateGeometryGen);
		return -1;
	}

	GLuint nfsTemplateGeometryGen = loadFragmentShader("TemplateGeometryMap.fs");
	glCompileShader(nfsTemplateGeometryGen);
	glGetShaderiv(nfsTemplateGeometryGen, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		shaderError("Error compiling TemplateGeometryMap.fs", nfsTemplateGeometryGen);
		return -1;
	}

	GLuint nprogTemplateGeometryGen = glCreateProgram();
	glAttachShader(nprogTemplateGeometryGen, nvsTemplateGeometryGen);
	glAttachShader(nprogTemplateGeometryGen, nfsTemplateGeometryGen);
	glLinkProgram(nprogTemplateGeometryGen);
	glGetProgramiv(nprogTemplateGeometryGen, GL_LINK_STATUS, &success);
	if (!success)
	{
		programError("Error linking nprogTemplateGeometryGen", nprogTemplateGeometryGen);
		return -1;
	}
	glValidateProgram(nprogTemplateGeometryGen);
	glGetProgramiv(nprogTemplateGeometryGen, GL_VALIDATE_STATUS, &success);
	if (!success)
	{
		programError("Error validating nprogTemplateGeometryGen", nprogTemplateGeometryGen);
		return -1;
	}

	return nprogTemplateGeometryGen;
}

void initGeometryMapIndices(int nWidth, int nHeight, GLuint* pIndices)
{
	int i, j, k;
	int nQuadsX = nWidth - 1;
	int nQuadsY = nHeight - 1;
	int nQuads = nQuadsX * nQuadsY;
	int nIndices = nQuads * 4;

	k = 0;
	for (i = 0; i < nQuadsY; i++)
	{
		for (j = 0; j < nQuadsX; j++)
		{
			pIndices[k]		= i * nWidth + j;
			pIndices[k + 1]	= (i + 1) * nWidth + j;
			pIndices[k + 2]	= (i + 1) * nWidth + j + 1;
			pIndices[k + 3]	= i * nWidth + j + 1;
			k += 4;
		}
	}
}

IplImage* generateTemplateGeometryMap(CBSplineGridWavelet<C3Vectorf>* wavelet, GLfloat * pGeomMapMask, double* pFaceVertices, GLuint & nprogTemplateGeometryGen, 
									  GLuint & nfb, GLuint & ntexTemplateGeometryMap, GLuint & ntexTemplateDepthMap, double * matTemplateStereoMap,
									  double & templateMapLeft, double & templateMapRight, double & templateMapBottom, double & templateMapTop, double * pTemplate2DPositions,
									  IplImage * pTemplate3DMap, GLuint * pGeomMapIndices)
{
	int i;

	double * matRigidFaceAlign = new double[16];

	for(i = 0; i < 16; i++) matRigidFaceAlign[i] = 0;
	matRigidFaceAlign[0] = matRigidFaceAlign[5] = matRigidFaceAlign[10] = matRigidFaceAlign[15] = 1.0;

	float matRigid[16];

	for(i = 0; i < 16; i++)
		matRigid[i] = (float)matRigidFaceAlign[i];

	fbBind(nfb);
	fbAttachTexture(nfb, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, ntexTemplateGeometryMap);
	fbAttachTexture(nfb, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, ntexTemplateDepthMap);
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	glViewport(0, 0, wavelet->getFullResWidth(), wavelet->getFullResHeight());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUseProgram(nprogTemplateGeometryGen);

	GLint iloc = glGetUniformLocation(nprogTemplateGeometryGen, "matRigidAlign");
	glUniformMatrix4fv(iloc, 1, true, matRigid);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadTransposeMatrixd(matTemplateStereoMap);
	glMultTransposeMatrixd(matRigidFaceAlign);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(templateMapLeft, templateMapRight, templateMapBottom, templateMapTop);

	//draw face with texture coordinates
	CFaceData::drawFaceWithTexCoords(pFaceVertices, pTemplate2DPositions);

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glFinish();
	glUseProgram(0);
	fbDetachAll(nfb);
	fbUnbind();

	uchar* pData;
	int nStep;
	cvGetRawData(pTemplate3DMap, &pData, &nStep);

	glBindTexture(GL_TEXTURE_2D, ntexTemplateGeometryMap);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGB, GL_FLOAT, pData);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_ALPHA, GL_FLOAT, pGeomMapMask);
	glBindTexture(GL_TEXTURE_2D, 0);

	//flip the image vertically (i.e. about x-axis)
//	cvFlip(pTemplate3DMap);

	float* pGeomMaskCopy = new float[pTemplate3DMap->height * pTemplate3DMap->width];
	_ASSERT(pGeomMaskCopy != NULL);
	memcpy(pGeomMaskCopy, pGeomMapMask, pTemplate3DMap->height * pTemplate3DMap->width * sizeof(float));

	const float mask_blend_eps = 0.2f;
	int x, y, blendedVerts, nNbrs;
	float w, weightTotal;
	abutil::C3Vectorf* pVerts = (abutil::C3Vectorf*)pData;
	abutil::C3Vectorf vHome;
	do
	{
		blendedVerts = 0;
		i = 0;
		for (y = 0; y < pTemplate3DMap->height; y++)
		{
			for (x = 0; x < pTemplate3DMap->width; x++)
			{
				if (pGeomMaskCopy[i] < 1.f - mask_blend_eps)
				{
					weightTotal = 0.f;
					vHome.set(0.f, 0.f, 0.f);
					nNbrs = 0;
					if (y > 0)
					{
						w = pGeomMaskCopy[i - pTemplate3DMap->width];
						vHome += pVerts[i - pTemplate3DMap->width] * w;
						weightTotal += w;
						nNbrs++;
					}
					if (x > 0)
					{
						w = pGeomMaskCopy[i - 1];
						vHome += pVerts[i - 1] * w;
						weightTotal += w;
						nNbrs++;
					}
					if (y < pTemplate3DMap->height - 1)
					{
						w = pGeomMaskCopy[i + pTemplate3DMap->width];
						vHome += pVerts[i + pTemplate3DMap->width] * w;
						weightTotal += w;
						nNbrs++;
					}
					if (x < pTemplate3DMap->width - 1)
					{
						w = pGeomMaskCopy[i + 1];
						vHome += pVerts[i + 1] * w;
						weightTotal += w;
						nNbrs++;
					}

					if (weightTotal > mask_blend_eps)
					{
						vHome /= weightTotal;
						weightTotal /= (float)nNbrs;

						pVerts[i] = vHome;
						pGeomMaskCopy[i] = 1.f;//weightTotal;
					}

					blendedVerts++;
				}

				i++;
			}
		}
	} while (blendedVerts > 0);

	//reset geometry map indices, then apply mask
	initGeometryMapIndices(pTemplate3DMap->width, pTemplate3DMap->height, pGeomMapIndices);
	maskGeometryMapIndices(pTemplate3DMap->width, pTemplate3DMap->height, pGeomMapMask, pGeomMapIndices);

	//set full resolution data of the wavelet
	wavelet->setFullResolutionData((C3Vectorf*)pData);

	delete [] matRigidFaceAlign;
	delete [] pGeomMaskCopy;

	return pTemplate3DMap;
}

void maskGeometryMapIndices(int nWidth, int nHeight, GLfloat* pMask, GLuint* pIndices)
{
	int i, j, k;
	int iv0, iv1, iv2, iv3;
	int nQuadsX = nWidth - 1;
	int nQuadsY = nHeight - 1;
	int nQuads = nQuadsX * nQuadsY;
	int nIndices = nQuads * 4;
	GLuint iMaskValue = CGlobalImageMatch::gim_indexMaskValue;

	k = 0;
	for (i = 0; i < nQuadsY; i++)
	{
		for (j = 0; j < nQuadsX; j++)
		{
			iv0 = pIndices[k];
			iv1 = pIndices[k + 1];
			iv2 = pIndices[k + 2];
			iv3 = pIndices[k + 3];
			if (pMask[iv0] < 0.5 || pMask[iv1] < 0.5 || pMask[iv2] < 0.5f || pMask[iv3] < 0.5f)
			{
				pIndices[k]		= iMaskValue;
				pIndices[k + 1]	= iMaskValue;
				pIndices[k + 2]	= iMaskValue;
				pIndices[k + 3] = iMaskValue;
			}
			k += 4;
		}
	}
}

void initDisplay()
{
	glewInit();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
	glDepthFunc(GL_LEQUAL);
	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0);

	checkGLError("initDisplay pre-shaders");

	if (!glewIsSupported( "GL_VERSION_2_0 GL_ARB_fragment_program GL_ARB_vertex_program GL_EXT_framebuffer_object" ))
	{
		printf("Minimal extensions not supported\n");
	}
}

void saveGeometryMapAsOFF(char* szFilename, IplImage* pTemplate3DMap, CBSplineGridWavelet<C3Vectorf>* pWavelet, GLuint* pGeomMapIndices)
{
	_ASSERT(pTemplate3DMap != NULL);
	_ASSERT(pWavelet != NULL);
	_ASSERT(pGeomMapIndices != NULL);
	_ASSERT(szFilename != NULL);

	float* pVertData;
	uchar* pData;
	int nStep, nWidth, nHeight, nVerts, nIndices;

	cvGetRawData(pTemplate3DMap, &pData, &nStep);
	pVertData = (float*)pData;

	nWidth = pWavelet->getFullResWidth();
	nHeight = pWavelet->getFullResHeight();
	nVerts = nHeight * nWidth;
	nIndices = (nHeight - 1) * (nWidth - 1) * 4;

	CFaceData::saveFaceOFFIndices(szFilename, nVerts, pVertData, nIndices, pGeomMapIndices, 4);
}

void initGeometryMapColors(CBSplineGridWavelet<C3Vectorf>* pStdDevShape, abutil::C3Vectorf* pStdDevRecon, abutil::C3Vectorf *& pGeomMapColors)
{
	_ASSERT(pStdDevShape != NULL);
	_ASSERT(pStdDevRecon != NULL);

	int i;
	int nPoints = pStdDevShape->getNumCoefficients();

	C3Vectorf std, color;
	float mag, magNorm, valMin = 1e30, valMax = -1e30, rcpRange;
	float logMax, logMin;

	if (pGeomMapColors == NULL)
	{
		pGeomMapColors = new C3Vectorf[nPoints];
		_ASSERT(pGeomMapColors != NULL);
	}

	pStdDevShape->reconstructCopy(pStdDevShape->getNumLevels() - 1, pStdDevRecon, pStdDevShape->getFullResWidth(), 1);

	for (i = 0; i < nPoints; i++)
	{
		std = pStdDevRecon[i];
		mag = std.length();
		if (mag < valMin)
			valMin = mag;
		if (mag > valMax)
			valMax = mag;
	}

	std::cout << "std dev range [" << valMin << ", " << valMax << "]\n";
	std::cout.flush();

	rcpRange = 1.f / (valMax - valMin);
	logMax = log10(valMax);
	logMin = log10(max(valMin, 1e-5));

	abutil::C3Vectorf red(1.f, 0.f, 0.f), green(0.f, 1.f, 0.f), blue(0.f, 0.f, 1.f);

	for (i = 0; i < nPoints; i++)
	{
		std = pStdDevRecon[i];
		mag = std.length();
		magNorm = rcpRange * (mag - valMin);

		color.set(0.f, 0.f, 0.f);
		if (magNorm < 0.5)
		{
			color = green * (2.f * magNorm) + blue * (1.f - 2.f * magNorm);
		}
		else
		{
			color = red * (2.f * (magNorm - 0.5f)) + green * (1.f - 2.f * (magNorm - 0.5f));
		}

		pGeomMapColors[i] = color;
	}
}

void initStereoPoints(std::vector<SStereoPatch> patches, float *& pStereoPoints, int & nStereoPoints, CNearestNeighborAssistant *& NNA)
{
	int i;

	if (pStereoPoints != NULL)
	{
		delete [] pStereoPoints;
		pStereoPoints = NULL;
		nStereoPoints = 0;
	}

	//this is where the check needs to go for what method was used to
	//get the initial stereo results, to decide where to initialize from
	//for now just assume Furukawa's method and initialize from patch set

	nStereoPoints = (int)patches.size();
	pStereoPoints = new float[nStereoPoints * 3];
	_ASSERT(pStereoPoints != NULL);

	for (i = 0; i < nStereoPoints; i++)
	{
		pStereoPoints[i * 3]		= patches[i].px;
		pStereoPoints[i * 3 + 1]	= patches[i].py;
		pStereoPoints[i * 3 + 2]	= patches[i].pz;
	}

	//now reset nearest neighbor assistant to use these points as reference points
	NNA->setReferencePoints(nStereoPoints, pStereoPoints);
}

void updateGeometryMap(bool bUpdateNormals, IplImage *& pTemplate3DMap, CBSplineGridWavelet<abutil::C3Vectorf> *& wavelet, GLuint & nbufGeomMap, CGlobalImageMatch & imageMatcher,
					   abutil::C3Vectorf *& pGeomMapNormals, GLuint & nbufGeomMapNormals, GLuint & ntexTemplateGeometryMap)
{
	int i, j, k, kup, kdown, kleft, kright, iprev, jprev, inext, jnext;
	uchar* pData;
	int nStep;
	C3Vectorf* pVerts;
	C3Vectorf v, vedge1, vedge2, n;

	cvGetRawData(pTemplate3DMap, &pData, &nStep);
	pVerts = (C3Vectorf*)pData;

	wavelet->reconstructCopy(wavelet->getNumLevels() - 1, pVerts, nStep / sizeof(C3Vectorf), 1);
	nbufGeomMap = imageMatcher.getVertexBuffer();
	glBindBuffer(GL_ARRAY_BUFFER, nbufGeomMap);
	glBufferSubData(GL_ARRAY_BUFFER, 0, wavelet->getNumCoefficients() * sizeof(C3Vectorf), pVerts);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	if (bUpdateNormals)
	{
		iprev = 0;
		inext = 1;
		for (i = 0; i < pTemplate3DMap->height; i++)
		{
			jprev = 0;
			jnext = 1;
			for (j = 0; j < pTemplate3DMap->width; j++)
			{
				k = i * pTemplate3DMap->width + j;
				kup = iprev * pTemplate3DMap->width + j;
				kleft = i * pTemplate3DMap->width + jprev;
				kdown = inext * pTemplate3DMap->width + j;
				kright = i * pTemplate3DMap->width + jnext;

				pGeomMapNormals[k].set(0.f, 0.f, 0.f);
				v = pVerts[k];
				if (kright != k && kup != k)
				{
					vedge1 = pVerts[kright] - v;
					vedge2 = pVerts[kup] - v;
					vedge1.crossFast(vedge2, n);
					n.normalize();
					pGeomMapNormals[k] += n;
				}
				if (kup != k && kleft != k)
				{
					vedge1 = pVerts[kup] - v;
					vedge2 = pVerts[kleft] - v;
					vedge1.crossFast(vedge2, n);
					n.normalize();
					pGeomMapNormals[k] += n;
				}
				if (kleft != k && kdown != k)
				{
					vedge1 = pVerts[kleft] - v;
					vedge2 = pVerts[kdown] - v;
					vedge1.crossFast(vedge2, n);
					n.normalize();
					pGeomMapNormals[k] += n;
				}
				if (kdown != k && kright != k)
				{
					vedge1 = pVerts[kdown] - v;
					vedge2 = pVerts[kright] - v;
					vedge1.crossFast(vedge2, n);
					n.normalize();
					pGeomMapNormals[k] += n;
				}
				pGeomMapNormals[k].normalize();

				jprev = j;
				if (++jnext >= pTemplate3DMap->width) jnext--;
			}

			iprev = i;
			if (++inext >= pTemplate3DMap->height) inext--;
		}

		nbufGeomMapNormals = imageMatcher.getNormalBuffer();
		glBindBuffer(GL_ARRAY_BUFFER, nbufGeomMapNormals);
		glBufferSubData(GL_ARRAY_BUFFER, 0, wavelet->getNumCoefficients() * sizeof(C3Vectorf), pGeomMapNormals);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	glBindTexture(GL_TEXTURE_2D, ntexTemplateGeometryMap);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, pTemplate3DMap->width, pTemplate3DMap->height, GL_RGB, GL_FLOAT, pData);
}

void initWaveletShapeSampler(CWaveletShapeSampler & WSS, CBSplineGridWavelet<abutil::C3Vectorf> *& pMeanShape, CBSplineGridWavelet<abutil::C3Vectorf> *& wavelet, IplImage *& pTemplate3DMap, 
							 GLuint & nbufGeomMap, CGlobalImageMatch & imageMatcher, abutil::C3Vectorf *& pGeomMapNormals, GLuint & nbufGeomMapNormals, GLuint & ntexTemplateGeometryMap,
							 CWaveletShapePrior *& pStatPrior, float *& pmatReductTranspose, CNearestNeighborEnergy *& pStereoDataNN, float *& pStereoPoints, int& nStereoPoints,
							 double *& matAlignModelToData, double *& matAlignDataToModel, CBSplineGridWavelet<C3Vectorf> *& pStdDevShape, GLfloat *& pGeomMapMask, CNearestNeighborAssistant *& NNA)
{
	//double matTemp[16];
	abutil::C4x4Matrixf matCopy, matAlign;
	
	if (pStatPrior == NULL)
	{
		wavelet->copy(pMeanShape);
		updateGeometryMap(true, pTemplate3DMap, wavelet, nbufGeomMap, imageMatcher, pGeomMapNormals,  nbufGeomMapNormals, ntexTemplateGeometryMap);

		pStatPrior = new CWaveletShapePrior();
		_ASSERT(pStatPrior != NULL);

		pStatPrior->m_pMean = pMeanShape;
		pStatPrior->m_pStdDev = pStdDevShape;
		pStatPrior->m_pmatReductTranspose = pmatReductTranspose;
		pStatPrior->m_weight = 1.f;
		pStatPrior->m_sampleRegion = SAMPLE_REGION;

		WSS.addPrior(pStatPrior);
	}

	WSS.setNumSamples(NUMBER_SAMPLES);
	WSS.setVertexMask(pGeomMapMask);

	if (pStereoDataNN == NULL)
	{
		pStereoDataNN = new CNearestNeighborEnergy();
		_ASSERT(pStereoDataNN != NULL);

		pStereoDataNN->m_pNNA = NNA;
		pStereoDataNN->m_pObservedPoints = pStereoPoints;
		pStereoDataNN->m_nObservedPoints = nStereoPoints;

		//NOTE: Distances here are computed in the coordinate system of the DATA, and NOT the model!
		//Here, consider points with distance up to a limit (in this case, the unit is in millimeters)
		pStereoDataNN->m_weight = 1.f;
		pStereoDataNN->m_truncThresh = 10.0*10.0; //the threshold is compared to the squared distance

		pStereoDataNN->m_matTransform[0]	= matAlignModelToData[0];
		pStereoDataNN->m_matTransform[1]	= matAlignModelToData[1];
		pStereoDataNN->m_matTransform[2]	= matAlignModelToData[2];
		pStereoDataNN->m_matTransform[3]	= matAlignModelToData[3];
		pStereoDataNN->m_matTransform[4]	= matAlignModelToData[4];
		pStereoDataNN->m_matTransform[5]	= matAlignModelToData[5];
		pStereoDataNN->m_matTransform[6]	= matAlignModelToData[6];
		pStereoDataNN->m_matTransform[7]	= matAlignModelToData[7];
		pStereoDataNN->m_matTransform[8]	= matAlignModelToData[8];
		pStereoDataNN->m_matTransform[9]	= matAlignModelToData[9];
		pStereoDataNN->m_matTransform[10]	= matAlignModelToData[10];
		pStereoDataNN->m_matTransform[11]	= matAlignModelToData[11];
		pStereoDataNN->m_matTransform[12]	= matAlignModelToData[12];
		pStereoDataNN->m_matTransform[13]	= matAlignModelToData[13];
		pStereoDataNN->m_matTransform[14]	= matAlignModelToData[14];
		pStereoDataNN->m_matTransform[15]	= matAlignModelToData[15];

		matAlign._11 = matAlignModelToData[0];	matAlign._12 = matAlignModelToData[1];	matAlign._13 = matAlignModelToData[2];	matAlign._14 = matAlignModelToData[3];
		matAlign._21 = matAlignModelToData[4];	matAlign._22 = matAlignModelToData[5];	matAlign._23 = matAlignModelToData[6];	matAlign._24 = matAlignModelToData[7];
		matAlign._31 = matAlignModelToData[8];	matAlign._32 = matAlignModelToData[9];	matAlign._33 = matAlignModelToData[10];	matAlign._34 = matAlignModelToData[11];
		matAlign._41 = matAlignModelToData[12];	matAlign._42 = matAlignModelToData[13];	matAlign._43 = matAlignModelToData[14];	matAlign._44 = matAlignModelToData[15];
		pStereoDataNN->m_matModelToData = matAlign;

		matAlign._11 = matAlignDataToModel[0];	matAlign._12 = matAlignDataToModel[1];	matAlign._13 = matAlignDataToModel[2];	matAlign._14 = matAlignDataToModel[3];
		matAlign._21 = matAlignDataToModel[4];	matAlign._22 = matAlignDataToModel[5];	matAlign._23 = matAlignDataToModel[6];	matAlign._24 = matAlignDataToModel[7];
		matAlign._31 = matAlignDataToModel[8];	matAlign._32 = matAlignDataToModel[9];	matAlign._33 = matAlignDataToModel[10];	matAlign._34 = matAlignDataToModel[11];
		matAlign._41 = matAlignDataToModel[12];	matAlign._42 = matAlignDataToModel[13];	matAlign._43 = matAlignDataToModel[14];	matAlign._44 = matAlignDataToModel[15];
		pStereoDataNN->m_matDataToModel = matAlign;

		WSS.addObservation(pStereoDataNN);
	}
}

bool readLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded)
{
	std::fstream fileStream;
	fileStream.open(sstrFileName.c_str(), std::ios::in);

	if(!fileStream.is_open())
	{
		std::cout << "Unable to open file " << sstrFileName << std::endl;
		return false;
	}

	std::vector<double> tmpLandmarks;

	while(true)
	{
		double value(0.0);
		fileStream >> value;
		tmpLandmarks.push_back(value);

		if(fileStream.eof() || fileStream.fail() || fileStream.bad())
		{
			break;
		}
	}

	fileStream.close();

	const size_t numLandmarks = (tmpLandmarks.size() - tmpLandmarks.size()%3)/3;

	loaded.clear();
	loaded.resize(numLandmarks);
	landmarks.clear();
	landmarks.resize(3*numLandmarks);

	for(int i = 0; i < numLandmarks; ++i)
	{
		loaded[i] = true;

		for(int j = 0; j < 3; ++j)
		{
			landmarks[3*i+j] = tmpLandmarks[3*i+j];
		}
	}

	return true;
}

void computeLSAlignmentModelToData(int nPoints, double* pModelPoints, double* pObservedPoints, double* pmatTransform, double* pmatTransformInv)
{
	const double SQRT3 = sqrt(3.0);

	int i, j, nCoords = nPoints * 3;
	long int dim, num, nrhs, nWork, info;
	char job, trans;
	double alpha, beta;

	job = 'A';
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;
	num = nPoints;
	dim = 4;
	nrhs = 3;
	nWork = 5 * nPoints;//this is enough for either dgels or dgesvd

	double* pNormModelPoints = new double[3 * nPoints];
	_ASSERT(pNormModelPoints != NULL);
	double* pNormObservedPoints = new double[3 * nPoints];
	_ASSERT(pNormObservedPoints != NULL);

	double* pA = new double[4 * nPoints];
	_ASSERT(pA != NULL);
	double* pB = new double[3 * nPoints];
	_ASSERT(pB != NULL);
	double* pWork = new double[nWork];
	_ASSERT(pWork != NULL);

	double R[9], S[3], U[9], VT[9];
	double matSimilModel[16], matSimilObservedInv[16], matSimilModelInv[16], matSimilObserved[16];
	double matTemp1[16], matTemp2[16], matTemp3[16], matTemp4[16];
	double x, y, z, cx, cy, cz, rmsdist, rmsscale;
	double rcpPoints = 1.0 / (double)nPoints;

	//normalize data points s.t. centroids are at origin and RMS distance (from origin) is sqrt(3)
	//save similarity transforms that do this
	//first model points
	cx = cy = cz = rmsdist = 0.0;

	//compute centroid
	for (i = 0; i < nCoords; i+=3)
	{
		cx += pModelPoints[i];
		cy += pModelPoints[i + 1];
		cz += pModelPoints[i + 2];
	}
	cx *= rcpPoints;
	cy *= rcpPoints;
	cz *= rcpPoints;

	//subtract centroid and compute RMS distance
	for (i = 0; i < nCoords; i+=3)
	{
		pNormModelPoints[i]		= x = pModelPoints[i] - cx;
		pNormModelPoints[i + 1]	= y = pModelPoints[i + 1] - cy;
		pNormModelPoints[i + 2]	= z = pModelPoints[i + 2] - cz;
		rmsdist += x*x + y*y + z*z;
	}
	rmsdist *= rcpPoints;
	rmsdist = sqrt(rmsdist);
	rmsscale = SQRT3 / rmsdist;

	//scale points
	for (i = 0; i < nCoords; i+=3)
	{
		pNormModelPoints[i]		*= rmsscale;
		pNormModelPoints[i + 1]	*= rmsscale;
		pNormModelPoints[i + 2]	*= rmsscale;
	}

	//"record" similarity transform
	matSimilModel[0]	= rmsscale;
	matSimilModel[1]	= 0.0;
	matSimilModel[2]	= 0.0;
	matSimilModel[3]	= 0.0;
	matSimilModel[4]	= 0.0;
	matSimilModel[5]	= rmsscale;
	matSimilModel[6]	= 0.0;
	matSimilModel[7]	= 0.0;
	matSimilModel[8]	= 0.0;
	matSimilModel[9]	= 0.0;
	matSimilModel[10]	= rmsscale;
	matSimilModel[11]	= 0.0;
	matSimilModel[12]	= -rmsscale * cx;
	matSimilModel[13]	= -rmsscale * cy;
	matSimilModel[14]	= -rmsscale * cz;
	matSimilModel[15]	= 1.0;

	matSimilModelInv[0]		= 1.0 / rmsscale;
	matSimilModelInv[1]		= 0.0;
	matSimilModelInv[2]		= 0.0;
	matSimilModelInv[3]		= 0.0;
	matSimilModelInv[4]		= 0.0;
	matSimilModelInv[5]		= 1.0 / rmsscale;
	matSimilModelInv[6]		= 0.0;
	matSimilModelInv[7]		= 0.0;
	matSimilModelInv[8]		= 0.0;
	matSimilModelInv[9]		= 0.0;
	matSimilModelInv[10]	= 1.0 / rmsscale;
	matSimilModelInv[11]	= 0.0;
	matSimilModelInv[12]	= cx;
	matSimilModelInv[13]	= cy;
	matSimilModelInv[14]	= cz;
	matSimilModelInv[15]	= 1.0;

	//now observed points
	cx = cy = cz = rmsdist = 0.0;

	//compute centroid
	for (i = 0; i < nCoords; i+=3)
	{
		cx += pObservedPoints[i];
		cy += pObservedPoints[i + 1];
		cz += pObservedPoints[i + 2];
	}
	cx *= rcpPoints;
	cy *= rcpPoints;
	cz *= rcpPoints;

	//subtract centroid and compute RMS distance
	for (i = 0; i < nCoords; i+=3)
	{
		pNormObservedPoints[i]		= x = pObservedPoints[i] - cx;
		pNormObservedPoints[i + 1]	= y = pObservedPoints[i + 1] - cy;
		pNormObservedPoints[i + 2]	= z = pObservedPoints[i + 2] - cz;
		rmsdist += x*x + y*y + z*z;
	}
	rmsdist *= rcpPoints;
	rmsdist = sqrt(rmsdist);
	rmsscale = SQRT3 / rmsdist;

	//scale points
	for (i = 0; i < nCoords; i+=3)
	{
		pNormObservedPoints[i]		*= rmsscale;
		pNormObservedPoints[i + 1]	*= rmsscale;
		pNormObservedPoints[i + 2]	*= rmsscale;
	}

	//"record" inverse similarity transform
	matSimilObservedInv[0]	= 1.0 / rmsscale;
	matSimilObservedInv[1]	= 0.0;
	matSimilObservedInv[2]	= 0.0;
	matSimilObservedInv[3]	= 0.0;
	matSimilObservedInv[4]	= 0.0;
	matSimilObservedInv[5]	= 1.0 / rmsscale;
	matSimilObservedInv[6]	= 0.0;
	matSimilObservedInv[7]	= 0.0;
	matSimilObservedInv[8]	= 0.0;
	matSimilObservedInv[9]	= 0.0;
	matSimilObservedInv[10]	= 1.0 / rmsscale;
	matSimilObservedInv[11]	= 0.0;
	matSimilObservedInv[12]	= cx;
	matSimilObservedInv[13]	= cy;
	matSimilObservedInv[14]	= cz;
	matSimilObservedInv[15]	= 1.0;

	matSimilObserved[0]		= rmsscale;
	matSimilObserved[1]		= 0.0;
	matSimilObserved[2]		= 0.0;
	matSimilObserved[3]		= 0.0;
	matSimilObserved[4]		= 0.0;
	matSimilObserved[5]		= rmsscale;
	matSimilObserved[6]		= 0.0;
	matSimilObserved[7]		= 0.0;
	matSimilObserved[8]		= 0.0;
	matSimilObserved[9]		= 0.0;
	matSimilObserved[10]	= rmsscale;
	matSimilObserved[11]	= 0.0;
	matSimilObserved[12]	= -rmsscale * cx;
	matSimilObserved[13]	= -rmsscale * cy;
	matSimilObserved[14]	= -rmsscale * cz;
	matSimilObserved[15]	= 1.0;

	//fill matrices for least-squares estimation
	for (i = 0; i < nPoints; i++)
	{
		j = i * 3;
		pA[i]				= pNormModelPoints[j];
		pA[nPoints + i]		= pNormModelPoints[j + 1];
		pA[2 * nPoints + i]	= pNormModelPoints[j + 2];
		pB[i]				= pNormObservedPoints[j];
		pB[nPoints + i]		= pNormObservedPoints[j + 1];
		pB[2 * nPoints + i]	= pNormObservedPoints[j + 2];
	}

	dim = 3;
	clapack::dgels_(&trans, &num, &dim, &nrhs, pA, &num, pB, &num, pWork, &nWork, &info);
	if (info != 0)
	{
		printf("computeLSAlignmentModelToData(...): clapack::dgels_(...) returned %i\n", info);
		goto computeLSAlignmentModelToData_EXIT;
	}

	//copy rotation separately for SVD
	R[0] = pB[0];			R[3] = pB[1];				R[6] = pB[2];
	R[1] = pB[nPoints];		R[4] = pB[nPoints + 1];		R[7] = pB[nPoints + 2];
	R[2] = pB[2 * nPoints];	R[5] = pB[2 * nPoints + 1];	R[8] = pB[2 * nPoints + 2];

	clapack::dgesvd_(&job, &job, &nrhs, &nrhs, R, &nrhs, S, U, &nrhs, VT, &nrhs, pWork, &nWork, &info);
	if (info != 0)
	{
		printf("computeLSAlignmentModelToData(...): clapack::dgesvd_(...) returned %i\n", info);
		goto computeLSAlignmentModelToData_EXIT;
	}

	//pure rotation R = U*VT
	clapack::dgemm_(&trans, &trans, &nrhs, &nrhs, &nrhs, &alpha, U, &nrhs, VT, &nrhs, &beta, R, &nrhs);

	matTemp1[0]		= R[0];
	matTemp1[1]		= R[1];
	matTemp1[2]		= R[2];
	matTemp1[3]		= 0.0;
	matTemp1[4]		= R[3];
	matTemp1[5]		= R[4];
	matTemp1[6]		= R[5];
	matTemp1[7]		= 0.0;
	matTemp1[8]		= R[6];
	matTemp1[9]		= R[7];
	matTemp1[10]	= R[8];
	matTemp1[11]	= 0.0;
	matTemp1[12]	= 0.0;
	matTemp1[13]	= 0.0;
	matTemp1[14]	= 0.0;
	matTemp1[15]	= 1.0;

	if (pmatTransformInv != NULL)
	{
		//transpose of rotation
		matTemp3[0]		= R[0];
		matTemp3[1]		= R[3];
		matTemp3[2]		= R[6];
		matTemp3[3]		= 0.0;
		matTemp3[4]		= R[1];
		matTemp3[5]		= R[4];
		matTemp3[6]		= R[7];
		matTemp3[7]		= 0.0;
		matTemp3[8]		= R[2];
		matTemp3[9]		= R[5];
		matTemp3[10]	= R[8];
		matTemp3[11]	= 0.0;
		matTemp3[12]	= 0.0;
		matTemp3[13]	= 0.0;
		matTemp3[14]	= 0.0;
		matTemp3[15]	= 1.0;
	}

	//right multiply by matSimilModel: matTemp2 = matTemp1 * matSimilModel
	clapack::integer p, r, q;
	p = 4;
	r = 4;
	q = 4;
	clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matTemp1, &q, matSimilModel, &q, &beta, matTemp2, &p);

	if (pmatTransformInv != NULL)
	{
		//right multiply transpose by matSimilObserved: matTemp4 = matTemp3 * matSimilObserved
		clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matTemp3, &q, matSimilObserved, &q, &beta, matTemp4, &p);
	}

	//left multiply by matSimilObservedInv: matTemp1 = matSimilObservedInv * matTemp2
	clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matSimilObservedInv, &p, matTemp2, &q, &beta, matTemp1, &p);

	if (pmatTransformInv != NULL)
	{
		//left multiply by matSimilModelInv: matTemp3 = matSimilModelInv * matTemp4
		clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, matSimilModelInv, &p, matTemp4, &q, &beta, matTemp3, &p);
	}

	//put into row-major order (not sure why I use row-major when OpenGL and CLAPACK use column-major, but...)
	pmatTransform[0]	= matTemp1[0];
	pmatTransform[1]	= matTemp1[4];
	pmatTransform[2]	= matTemp1[8];
	pmatTransform[3]	= matTemp1[12];
	pmatTransform[4]	= matTemp1[1];
	pmatTransform[5]	= matTemp1[5];
	pmatTransform[6]	= matTemp1[9];
	pmatTransform[7]	= matTemp1[13];
	pmatTransform[8]	= matTemp1[2];
	pmatTransform[9]	= matTemp1[6];
	pmatTransform[10]	= matTemp1[10];
	pmatTransform[11]	= matTemp1[14];
	pmatTransform[12]	= matTemp1[3];
	pmatTransform[13]	= matTemp1[7];
	pmatTransform[14]	= matTemp1[11];
	pmatTransform[15]	= matTemp1[15];
	if (pmatTransformInv != NULL)
	{
		pmatTransformInv[0]		= matTemp3[0];
		pmatTransformInv[1]		= matTemp3[4];
		pmatTransformInv[2]		= matTemp3[8];
		pmatTransformInv[3]		= matTemp3[12];
		pmatTransformInv[4]		= matTemp3[1];
		pmatTransformInv[5]		= matTemp3[5];
		pmatTransformInv[6]		= matTemp3[9];
		pmatTransformInv[7]		= matTemp3[13];
		pmatTransformInv[8]		= matTemp3[2];
		pmatTransformInv[9]		= matTemp3[6];
		pmatTransformInv[10]	= matTemp3[10];
		pmatTransformInv[11]	= matTemp3[14];
		pmatTransformInv[12]	= matTemp3[3];
		pmatTransformInv[13]	= matTemp3[7];
		pmatTransformInv[14]	= matTemp3[11];
		pmatTransformInv[15]	= matTemp3[15];
	}

	printf("%10lf %10lf %10lf %10lf\n", pmatTransform[0], pmatTransform[1], pmatTransform[2], pmatTransform[3]);
	printf("%10lf %10lf %10lf %10lf\n", pmatTransform[4], pmatTransform[5], pmatTransform[6], pmatTransform[7]);
	printf("%10lf %10lf %10lf %10lf\n", pmatTransform[8], pmatTransform[9], pmatTransform[10], pmatTransform[11]);
	printf("%10lf %10lf %10lf %10lf\n", pmatTransform[12], pmatTransform[13], pmatTransform[14], pmatTransform[15]);

	if (pmatTransformInv != NULL)
	{
		printf("%10lf %10lf %10lf %10lf\n", pmatTransformInv[0], pmatTransformInv[1], pmatTransformInv[2], pmatTransformInv[3]);
		printf("%10lf %10lf %10lf %10lf\n", pmatTransformInv[4], pmatTransformInv[5], pmatTransformInv[6], pmatTransformInv[7]);
		printf("%10lf %10lf %10lf %10lf\n", pmatTransformInv[8], pmatTransformInv[9], pmatTransformInv[10], pmatTransformInv[11]);
		printf("%10lf %10lf %10lf %10lf\n", pmatTransformInv[12], pmatTransformInv[13], pmatTransformInv[14], pmatTransformInv[15]);

		trans = 'T';
		clapack::dgemm_(&trans, &trans, &p, &r, &q, &alpha, pmatTransform, &p, pmatTransformInv, &q, &beta, matTemp3, &p);
		printf("%10lf %10lf %10lf %10lf\n", matTemp3[0], matTemp3[1], matTemp3[2], matTemp3[3]);
		printf("%10lf %10lf %10lf %10lf\n", matTemp3[4], matTemp3[5], matTemp3[6], matTemp3[7]);
		printf("%10lf %10lf %10lf %10lf\n", matTemp3[8], matTemp3[9], matTemp3[10], matTemp3[11]);
		printf("%10lf %10lf %10lf %10lf\n", matTemp3[12], matTemp3[13], matTemp3[14], matTemp3[15]);
	}

computeLSAlignmentModelToData_EXIT:
	delete [] pNormModelPoints;
	delete [] pNormObservedPoints;
	delete [] pA;
	delete [] pB;
	delete [] pWork;
}

bool computeLSAlignmentModelToData(const std::vector<double>& modelLandmarks, const std::vector<bool>& modelLandmarksLoaded, const std::vector<double>& dataLandmarks, const std::vector<bool>& dataLandmarksLoaded
											, double* pmatTransform, double* pmatTransformInv)
{
	int numLandmarks(0);

	const size_t numAlignmentLandmarks = std::min<size_t>(8, std::min<size_t>(modelLandmarksLoaded.size(), dataLandmarksLoaded.size()));
	for(size_t i = 0; i < numAlignmentLandmarks; ++i)
	{
		if(modelLandmarksLoaded[i] && dataLandmarksLoaded[i])
		{
			++numLandmarks;
		}
	}

	double* pModelPoints = new double[3*numLandmarks]; 
	double* pObservedPoints = new double[3*numLandmarks];

	int currIndex(0);
	for(size_t i = 0; i < numAlignmentLandmarks; ++i)
	{
		if(!modelLandmarksLoaded[i] || !dataLandmarksLoaded[i])
		{
			continue;
		}

		pModelPoints[3*currIndex] = modelLandmarks[3*i];
		pModelPoints[3*currIndex+1] = modelLandmarks[3*i+1];
		pModelPoints[3*currIndex+2] = modelLandmarks[3*i+2];

		pObservedPoints[3*currIndex] = dataLandmarks[3*i];
		pObservedPoints[3*currIndex+1] = dataLandmarks[3*i+1];
		pObservedPoints[3*currIndex+2] = dataLandmarks[3*i+2];

		++currIndex;
	}

	computeLSAlignmentModelToData(numLandmarks, pModelPoints, pObservedPoints, pmatTransform, pmatTransformInv);

	delete [] pModelPoints;
	delete [] pObservedPoints;
	return true;
}

void computeLSAlignmentModelToDataFromNN(double *& pmatTransform, double *& pmatTransformInv, IplImage *& pTemplate3DMap, CNearestNeighborAssistant * NNA, GLfloat*& pGeomMapMask,
										 CBSplineGridWavelet<abutil::C3Vectorf> * wavelet, float *& pStereoPoints)
{
	int i, j, k;
	double matInit[16];
	uchar* pData;
	C3Vectorf* pverts;
	
	memcpy(matInit, pmatTransform, 16 * sizeof(double));

	int nPoints = wavelet->getNumCoefficients();
	int nCoords = nPoints * 3;
	int nUsePoints;

	NNAReal* pModel = new NNAReal[nCoords];
	_ASSERT(pModel != NULL);
	NNAReal* pDistances = new NNAReal[nPoints];
	_ASSERT(pDistances != NULL);
	NNAIndex* pIndices = new NNAIndex[nPoints];
	_ASSERT(pIndices != NULL);

	double* pModel_d = new double[nCoords];
	_ASSERT(pModel_d != NULL);
	double* pObserved = new double[nCoords];
	_ASSERT(pObserved != NULL);

	cvGetRawData(pTemplate3DMap, &pData);
	pverts = (C3Vectorf*)pData;

	for (i = 0, j = 0; i < nPoints; i++, j+=3)
	{
		pModel[j]		= matInit[0] * pverts[i].x + matInit[1] * pverts[i].y + matInit[2] * pverts[i].z + matInit[3];
		pModel[j + 1]	= matInit[4] * pverts[i].x + matInit[5] * pverts[i].y + matInit[6] * pverts[i].z + matInit[7];
		pModel[j + 2]	= matInit[8] * pverts[i].x + matInit[9] * pverts[i].y + matInit[10] * pverts[i].z + matInit[11];
	}

	NNA->setQueryPoints(nPoints, pModel, pDistances, pIndices);
	NNA->compute();
	NNA->distanceStatsHistogram(1000, "nndist.txt");

	nUsePoints = 0;
	j = 0;
	for (i = 0; i < nPoints; i++)
	{
		//NOTE: Distances here are computed in the coordinate system of the DATA, and NOT the model!
		//Here, consider points with distance up to 1cm (in this case, the unit is in millimeters)
		if (pGeomMapMask[i] > 0.5 && pDistances[i] < 10.0) // 0.001)
		{
			pModel_d[j]		= pverts[i].x;
			pModel_d[j + 1]	= pverts[i].y;
			pModel_d[j + 2]	= pverts[i].z;

			k = pIndices[i];
			pObserved[j]		= pStereoPoints[k * 3];
			pObserved[j + 1]	= pStereoPoints[k * 3 + 1];
			pObserved[j + 2]	= pStereoPoints[k * 3 + 2];

			nUsePoints++;
			j += 3;
		}
	}

	printf("Aligning model using %i points nearest neighbors...\n", nUsePoints);
	computeLSAlignmentModelToData(nUsePoints, pModel_d, pObserved, pmatTransform, pmatTransformInv);

	delete [] pModel;
	delete [] pDistances;
	delete [] pIndices;
	delete [] pModel_d;
	delete [] pObserved;
}
