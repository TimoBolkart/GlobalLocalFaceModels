///////////////////////////////////////////////////////////////////////////////
//
//	BSGWFactory.cpp
//
//	Source file for the CBSGWFactory class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#include "StereoFace.h"
#include "BSGWFactory.h"

void* CBSGWFactory::createBSGW(EBSGWType eType, int nBaseWidth, int nBaseHeight, int nLevels)
{
	switch (eType)
	{
	case eBSGWFloat:
		return (void*)(new CBSplineGridWavelet<float>(nBaseWidth, nBaseHeight, nLevels, eType));
	case eBSGWDouble:
		return (void*)(new CBSplineGridWavelet<double>(nBaseWidth, nBaseHeight, nLevels, eType));
	case eBSGWFloat3:
		return (void*)(new CBSplineGridWavelet<abutil::C3Vectorf>(nBaseWidth, nBaseHeight, nLevels, eType));
	}
	
	return NULL;
}

EBSGWType CBSGWFactory::loadBSGW(char *szFilename, void **pbsgw)
{
	EBSGWType eType = eBSGWFloat;
	int nbw, nbh, nl;

	CBSplineGridWavelet<float>* pbsgwf;
	CBSplineGridWavelet<double>* pbsgwd;
	CBSplineGridWavelet<abutil::C3Vectorf>* pbsgwf3;

	FILE* pf = fopen(szFilename, "rb");
	if (pf == NULL)
	{
		*pbsgw = NULL;
		return eType;
	}

	fread(&nbw, sizeof(int), 1, pf);
	fread(&nbh, sizeof(int), 1, pf);
	fread(&nl, sizeof(int), 1, pf);
	fread(&eType, sizeof(EBSGWType), 1, pf);

	switch (eType)
	{
	case eBSGWFloat:
		pbsgwf = new CBSplineGridWavelet<float>(nbw, nbh, nl, eType);
		if (pbsgwf != NULL)
			fread(pbsgwf->getCoefficientStore(), sizeof(float), pbsgwf->getFullResHeight() * pbsgwf->getFullResWidth(), pf);
		*pbsgw = (void*)pbsgwf;
		break;
	case eBSGWDouble:
		pbsgwd = new CBSplineGridWavelet<double>(nbw, nbh, nl, eType);
		if (pbsgwd != NULL)
			fread(pbsgwd->getCoefficientStore(), sizeof(double), pbsgwd->getFullResHeight() * pbsgwd->getFullResWidth(), pf);
		*pbsgw = (void*)pbsgwd;
		break;
	case eBSGWFloat3:
		pbsgwf3 = new CBSplineGridWavelet<abutil::C3Vectorf>(nbw, nbh, nl, eType);
		if (pbsgwf3 != NULL)
			fread(pbsgwf3->getCoefficientStore(), sizeof(abutil::C3Vectorf), pbsgwf3->getFullResHeight() * pbsgwf3->getFullResWidth(), pf);
		*pbsgw = (void*)pbsgwf3;
		break;
	}

	fclose(pf);

	return eType;
}

//void* CBSGWFactory::createBSGW_cuda(EBSGWType eType, int nBaseWidth, int nBaseHeight, int nLevels)
//{
//#if STEREOFACE_USECUDA
//	switch (eType)
//	{
//	case eBSGWFloat3:
//		return (void*)(new CCUDABSplineGridWavelet3f(nBaseWidth, nBaseHeight, nLevels));
//	default:
//		printf("Type not supported\n");
//	}
//#endif //STEREOFACE_USECUDA
//
//	return NULL;
//}
//
//EBSGWType CBSGWFactory::loadBSGW_cuda(char* szFilename, void** pbsgw)
//{
//	EBSGWType eType = eBSGWFloat3;
//
//#if STEREOFACE_USECUDA
//#endif //STEREOFACE_USECUDA
//
//	return eType;
//}



