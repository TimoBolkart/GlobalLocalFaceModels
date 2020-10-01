///////////////////////////////////////////////////////////////////////////////
//
//	NearestNeighborEnergy.h
//
//	Header file for the CNearestNeighborEnergy class
//	subclass of CWaveletShapeObesrvation
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __NEARESTNEIGHBORENERGY_H__
#define __NEARESTNEIGHBORENERGY_H__


#include "WaveletShapeSampler.h"
#include "NearestNeighborAssistant.h"


//subclass for nearest neighbor distance-based energy
class CNearestNeighborEnergy : public CWaveletShapeObservation
{
public:

	CNearestNeighborAssistant*				m_pNNA;
	NNAReal*								m_pObservedPoints;
	int										m_nObservedPoints;

	C3Vectorf*								m_pTransformedGeometry;

	NNAReal*								m_pDistances;
	NNAIndex*								m_pIndices;

	float									m_truncThresh;

	float									m_matTransform[16];
	abutil::C4x4Matrixf						m_matModelToData;
	abutil::C4x4Matrixf						m_matDataToModel;

	CProfile								m_profUpdate;
	CProfile								m_profCompute;

	CNearestNeighborEnergy(): CWaveletShapeObservation()
	{
		m_pNNA = NULL;
		m_pObservedPoints = NULL;
		m_nObservedPoints = 0;
		m_pTransformedGeometry = NULL;
		m_pDistances = NULL;
		m_pIndices = NULL;
		m_truncThresh = 10.f;
		m_matTransform[0]	= 1.f;
		m_matTransform[1]	= 0.f;
		m_matTransform[2]	= 0.f;
		m_matTransform[3]	= 0.f;
		m_matTransform[4]	= 0.f;
		m_matTransform[5]	= 1.f;
		m_matTransform[6]	= 0.f;
		m_matTransform[7]	= 0.f;
		m_matTransform[8]	= 0.f;
		m_matTransform[9]	= 0.f;
		m_matTransform[10]	= 1.f;
		m_matTransform[11]	= 0.f;
		m_matTransform[12]	= 0.f;
		m_matTransform[13]	= 0.f;
		m_matTransform[14]	= 0.f;
		m_matTransform[15]	= 1.f;
		m_matModelToData.setIdentity();
		m_matDataToModel.setIdentity();

		m_profUpdate.setName("CNearestNeighborEnergy::updateObserveration()");
		m_profCompute.setName("CNearestNeighborEnergy::computeDataCost()");
	}

	~CNearestNeighborEnergy()
	{
		if (m_pTransformedGeometry != NULL)
		{
			delete [] m_pTransformedGeometry;
			m_pTransformedGeometry = NULL;
		}
		if (m_pDistances != NULL)
		{
			delete [] m_pDistances;
			m_pDistances = NULL;
		}
		if (m_pIndices != NULL)
		{
			delete [] m_pIndices;
			m_pIndices = NULL;
		}
	}

	void init()
	{
		_ASSERT(m_pGeometry != NULL);

		int i;
		int nVerts = m_nVertices;

		m_pTransformedGeometry = new C3Vectorf[nVerts];
		_ASSERT(m_pTransformedGeometry != NULL);

		m_pDistances = new NNAReal[nVerts];
		_ASSERT(m_pDistances != NULL);

		m_pIndices = new NNAIndex[nVerts];
		_ASSERT(m_pIndices != NULL);

		for (i = 0; i < nVerts; i++)
		{
			m_pTransformedGeometry[i].x = m_matTransform[0] * m_pGeometry[i].x + m_matTransform[1] * m_pGeometry[i].y + m_matTransform[2] * m_pGeometry[i].z + m_matTransform[3];
			m_pTransformedGeometry[i].y = m_matTransform[4] * m_pGeometry[i].x + m_matTransform[5] * m_pGeometry[i].y + m_matTransform[6] * m_pGeometry[i].z + m_matTransform[7];
			m_pTransformedGeometry[i].z = m_matTransform[8] * m_pGeometry[i].x + m_matTransform[9] * m_pGeometry[i].y + m_matTransform[10] * m_pGeometry[i].z + m_matTransform[11];
		}

		m_pNNA->setQueryPoints(nVerts, (NNAReal*)m_pTransformedGeometry, m_pDistances, m_pIndices);
		m_pNNA->compute();
	}

	void updateObserveration()
	{
		_ASSERT(m_pNNA != NULL);
		_ASSERT(m_pObservedPoints != NULL);

#if WSS_PROFILE
		m_profUpdate.start();
#endif

		int i;
		int nVerts = m_nVertices;
		int nWidth = m_nWidth;
		int nHeight = m_nHeight;

		for (i = 0; i < nVerts; i++)
		{
			m_pTransformedGeometry[i].x = m_matTransform[0] * m_pGeometry[i].x + m_matTransform[1] * m_pGeometry[i].y + m_matTransform[2] * m_pGeometry[i].z + m_matTransform[3];
			m_pTransformedGeometry[i].y = m_matTransform[4] * m_pGeometry[i].x + m_matTransform[5] * m_pGeometry[i].y + m_matTransform[6] * m_pGeometry[i].z + m_matTransform[7];
			m_pTransformedGeometry[i].z = m_matTransform[8] * m_pGeometry[i].x + m_matTransform[9] * m_pGeometry[i].y + m_matTransform[10] * m_pGeometry[i].z + m_matTransform[11];
		}

		m_pNNA->setQueryPoints(nVerts, (NNAReal*)m_pTransformedGeometry, m_pDistances, m_pIndices);
		m_pNNA->compute();

#if WSS_PROFILE
		m_profUpdate.stop();
#endif 
	}

	float computeDataCost()
	{
		_ASSERT(m_pGeometry != NULL);
		_ASSERT(m_pObservedPoints != NULL);
		_ASSERT(m_pDistances != NULL);
		_ASSERT(m_pIndices != NULL);

#if WSS_PROFILE
		m_profCompute.start();
#endif

		int i, j, nPoints = m_nVertices;
		float dataCost = 0.f;

		C3Vectorf modelPoint, obsPoint, displacement;
		int counter(0);

		if (m_pGeometryMask != NULL)
		{
			for (i = 0; i < nPoints; i++)
			{
				if (m_pGeometryMask[i] > 0.5)
				{
					++counter;
					j = m_pIndices[i];
#if 0
					modelPoint = m_pGeometry[i];
					m_pTransformedGeometry[i].x = m_matTransform[0] * modelPoint.x + m_matTransform[1] * modelPoint.y + m_matTransform[2] * modelPoint.z + m_matTransform[3];
					m_pTransformedGeometry[i].y = m_matTransform[4] * modelPoint.x + m_matTransform[5] * modelPoint.y + m_matTransform[6] * modelPoint.z + m_matTransform[7];
					m_pTransformedGeometry[i].z = m_matTransform[8] * modelPoint.x + m_matTransform[9] * modelPoint.y + m_matTransform[10] * modelPoint.z + m_matTransform[11];
					modelPoint = m_pTransformedGeometry[i];

					obsPoint = *((C3Vectorf*)(m_pObservedPoints + j * 3));
					displacement = modelPoint - obsPoint;
					dataCost += min(m_truncThresh, displacement.lengthSquared());
#else
					obsPoint = *((C3Vectorf*)(m_pObservedPoints + j * 3));
					obsPoint = m_matDataToModel * obsPoint;
					modelPoint = m_pGeometry[i];
					displacement = modelPoint - obsPoint;
					dataCost += min(m_truncThresh, displacement.length());
#endif
				}
			}
		}
		else
		{
			for (i = 0; i < nPoints; i++)
			{
				j = m_pIndices[i];
				modelPoint = m_pGeometry[i];
				m_pTransformedGeometry[i].x = m_matTransform[0] * modelPoint.x + m_matTransform[1] * modelPoint.y + m_matTransform[2] * modelPoint.z + m_matTransform[3];
				m_pTransformedGeometry[i].y = m_matTransform[4] * modelPoint.x + m_matTransform[5] * modelPoint.y + m_matTransform[6] * modelPoint.z + m_matTransform[7];
				m_pTransformedGeometry[i].z = m_matTransform[8] * modelPoint.x + m_matTransform[9] * modelPoint.y + m_matTransform[10] * modelPoint.z + m_matTransform[11];
				modelPoint = m_pTransformedGeometry[i];

				obsPoint = *((C3Vectorf*)(m_pObservedPoints + j * 3));
				displacement = modelPoint - obsPoint;
				dataCost += min(m_truncThresh, displacement.lengthSquared());
			}
		}

#if WSS_PROFILE
		m_profCompute.stop();
#endif

		return dataCost;
	}

	void reportProfiling(FILE* pfOutput)
	{
		m_profUpdate.computeStatistics();
		m_profCompute.computeStatistics();
		m_profUpdate.reportStatistics(pfOutput);
		m_profCompute.reportStatistics(pfOutput);
	}
};


#endif //__NEARESTNEIGHBORENERGY_H__


