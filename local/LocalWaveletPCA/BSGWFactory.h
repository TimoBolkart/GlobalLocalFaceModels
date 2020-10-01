///////////////////////////////////////////////////////////////////////////////
//
//	BSGWFactory.h
//
//	Header file for the CBSGWFactory class
//	responsible for creates CBSplineGridWavelet objects
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __BSGWFACTORY_H__
#define __BSGWFACTORY_H__


#include "BSplineGridWavelet.h"

class CBSGWFactory
{
public:
	static void* createBSGW(EBSGWType eType, int nBaseWidth, int nBaseHeight, int nLevels);
	static EBSGWType loadBSGW(char* szFilename, void** pbsgw);
};


#endif //__BSGWFACTORY_H__


