///////////////////////////////////////////////////////////////////////////////
//
//	WaveletShapeSampler.h
//
//	Header file for the CWaveletShapeSampler class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __WAVELETSHAPESAMPLER_H__
#define __WAVELETSHAPESAMPLER_H__


#include "WaveletShapeFitter.h"


#define WSS_PROFILE							1

enum EWaveletShapeSamplingType
{
	eWSSUniform = 0,
	eWSSRandom,
	eWSSNumSamplingTypes
};


//main class in sampling-based optimization framework
class CWaveletShapeSampler : public CWaveletShapeFitter
{
protected:

	CBSplineGridWavelet<C3Vectorf>*			m_pActiveWaveletModel;

	int										m_nSamples;
	EWaveletShapeSamplingType				m_eSamplingType;

	//profiling objects
	CProfile								m_profOptimize;
	CProfile								m_profCoefficient;
	CProfile								m_profPriors;
	CProfile								m_profWaveletRecon;

	void initActiveWaveletModel();
	void clear();

	//void reconstructIntermediateModel()
	//{
	//	int j;
	//	for (j = 0; j < m_pActiveWaveletModel->getNumCoefficients(); j++)
	//	{
	//		m_priors[0]->inverseTransformCoefficient(j, m_pActiveWaveletModel->getCoefficient(j), m_pIntermediateModel->getCoefficientStore()[j]);
	//	}
	//}
	//void reconstructIntermediateModel(int iCoeff)
	//{
	//	m_priors[0]->inverseTransformCoefficient(iCoeff, m_pActiveWaveletModel->getCoefficient(iCoeff), m_pIntermediateModel->getCoefficientStore()[iCoeff]);
	//}


public:

	CWaveletShapeSampler(): CWaveletShapeFitter()
	{
		m_pActiveWaveletModel = NULL;
		m_pIntermediateModel = NULL;
		m_pReconVerts = NULL;
		m_pVertexMask = NULL;
		m_nSamples = 0;
		m_eSamplingType = eWSSUniform;

		m_profOptimize.setName("CWaveletShapeSampler::optimizeActiveWaveletModel()");
		m_profCoefficient.setName("CWaveletShapeSampler:: optimize individual coefficient");
		m_profPriors.setName("CWaveletShapeSampler:: evaluate priors");
		m_profWaveletRecon.setName("CWaveletShapeSampler::m_pIntermediateModel->reconstructCopy(...)");
	}

	~CWaveletShapeSampler()
	{
		clear();
	}

	void setNumSamples(int nSamples)		{ m_nSamples = nSamples; }

	void setVertexMask(float* pMask)		{ m_pVertexMask = pMask; }

	int getLevel()							{ return m_level; }

	void addObservation(CWaveletShapeObservation* pObs)
	{
		_ASSERT(pObs != NULL);
		pObs->m_nVertices = m_pActiveWaveletModel->getNumCoefficients();
		pObs->m_nWidth = m_pActiveWaveletModel->getFullResWidth();
		pObs->m_nHeight = m_pActiveWaveletModel->getFullResHeight();
		pObs->m_pGeometry = m_pReconVerts;
		pObs->m_pGeometryMask = m_pVertexMask;
		pObs->init();
		m_observations.push_back(pObs);
	}

	void reportProfiling(FILE* pfOutput)
	{
		int i;

		m_profOptimize.computeStatistics();
		m_profCoefficient.computeStatistics();
		m_profPriors.computeStatistics();
		m_profWaveletRecon.computeStatistics();
		m_profOptimize.reportStatistics(pfOutput);
		m_profCoefficient.reportStatistics(pfOutput);
		m_profPriors.reportStatistics(pfOutput);
		m_profWaveletRecon.reportStatistics(pfOutput);

		for (i = 0; i < (int)m_observations.size(); i++)
			m_observations[i]->reportProfiling(pfOutput);
		for (i = 0; i < (int)m_refinements.size(); i++)
			m_refinements[i]->reportProfiling(pfOutput);
	}

	void optimizeActiveWaveletModel();
};


#endif //__WAVELETSHAPESAMPLER_H__


