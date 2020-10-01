///////////////////////////////////////////////////////////////////////////////
//
//	NearestNeighborAssistant.cpp
//
//	Source file for the CNearestNeighborAssistant class
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#include "NearestNeighborAssistant.h"

///////////////////////////////////////////////////////////////////////////////
//external functions 
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//member functions
///////////////////////////////////////////////////////////////////////////////

void CNearestNeighborAssistant::clearInternalMemory()
{
	int i;

	if (m_pkdtree != NULL)
	{
		delete m_pkdtree;
		m_pkdtree = NULL;
	}
	if (m_pRefPoints != NULL)
	{
		for (i = 0; i < m_nRefPoints; i++)
		{
			if (m_pRefPoints[i] != NULL)
			{
				delete [] m_pRefPoints[i];
			}
		}
		delete [] m_pRefPoints;
		m_pRefPoints = NULL;
	}
	if (m_pQueryPoints != NULL)
	{
		delete [] m_pQueryPoints;
		m_pQueryPoints = NULL;
	}
}

void CNearestNeighborAssistant::setReferencePoints(int nPoints, NNAReal *pPoints)
{
	_ASSERT(pPoints != NULL);

	int i, j, iread, iwrite;

	//setting new reference points has the effect of essentially resetting the class
	//setQueryPoints must be called again
	clearInternalMemory();

	m_nRefPoints = nPoints;

	m_pRefPoints = new ANNpoint[m_nRefPoints];
	_ASSERT(m_pRefPoints != NULL);

	for (i = 0; i < m_nRefPoints; i++)
	{
		m_pRefPoints[i] = new ANNcoord[m_dim];
		_ASSERT(m_pRefPoints[i] != NULL);

		for (j = 0; j < m_dim; j++)
		{
			m_pRefPoints[i][j] = pPoints[i * m_dim + j];
		}
	}

	m_pkdtree = new ANNkd_tree(m_pRefPoints, m_nRefPoints, m_dim);
	_ASSERT(m_pkdtree != NULL);
}

void CNearestNeighborAssistant::setQueryPoints(int nPoints, NNAReal *pPoints, NNAReal* pDist, NNAIndex* pIndices, int nNNPerQuery)
{
	_ASSERT(pDist != NULL);
	_ASSERT(pIndices != NULL);
	_ASSERT(pPoints != NULL);

	int i, j, iread, iwrite;

	m_pDistanceOutput = pDist;
	m_pIndexOutput = pIndices;

	m_nNNPerQuery = nNNPerQuery;

	if (m_pQueryPoints != NULL && nPoints != m_nQueryPoints)
	{
		delete [] m_pQueryPoints;
		m_pQueryPoints = NULL;
	}

	if (m_pQueryPoints == NULL)
	{
		m_nQueryPoints = nPoints;
		m_pQueryPoints = new ANNcoord[m_dim * m_nQueryPoints];
		_ASSERT(m_pQueryPoints != NULL);
	}

	for (i = 0; i < m_nQueryPoints; i++)
	{
		for (j = 0; j < m_dim; j++)
		{
			m_pQueryPoints[i * m_dim + j] = pPoints[i * m_dim + j];
		}
	}
}

void CNearestNeighborAssistant::compute()
{
	_ASSERT(m_pRefPoints != NULL);
	_ASSERT(m_pQueryPoints != NULL);
	_ASSERT(m_pDistanceOutput != NULL);
	_ASSERT(m_pIndexOutput != NULL);

	ANNdist* pDistances = new ANNdist[m_nQueryPoints * m_nNNPerQuery];
	_ASSERT(pDistances != NULL);

	int i, k;
	for (i = 0; i < m_nQueryPoints; i++)
	{
		m_pkdtree->annkSearch(m_pQueryPoints + i * m_dim, m_nNNPerQuery, m_pIndexOutput + i * m_nNNPerQuery, pDistances + i * m_nNNPerQuery);

		for (k = 0; k < m_nNNPerQuery; k++)
			m_pDistanceOutput[i * m_nNNPerQuery + k] = sqrt(pDistances[i * m_nNNPerQuery + k]);
	}

	delete [] pDistances;
}

void CNearestNeighborAssistant::distanceStatsHistogram(NNAReal *pmean, NNAReal *pstd, NNAReal *pmin, NNAReal *pmax, int nBins, NNAReal* pBinDist, int *pBinCount)
{
	_ASSERT(pmean != NULL);
	_ASSERT(pstd != NULL);
	_ASSERT(pmin != NULL);
	_ASSERT(pmax != NULL);
	_ASSERT(pBinDist != NULL);
	_ASSERT(pBinCount != NULL);

	int i, iBin;
	int nPoints;
	NNAReal distMean, distStd, distMin, distMax;
	NNAReal diff, var, binSize, rcpBinSize;

	nPoints = m_nQueryPoints;

	distMean = m_pDistanceOutput[0];
	distMin = m_pDistanceOutput[0];
	distMax = m_pDistanceOutput[0];
	for (i = 1; i < m_nQueryPoints; i++)
	{
		if (_isnan((double)m_pDistanceOutput[i]))
		{
			nPoints--;
			continue;
		}

		distMean += m_pDistanceOutput[i];
		if (m_pDistanceOutput[i] < distMin)
			distMin = m_pDistanceOutput[i];
		if (m_pDistanceOutput[i] > distMax)
			distMax = m_pDistanceOutput[i];
	}

	distMean /= (NNAReal)nPoints;

	var = 0.0;
	for (i = 0; i < m_nQueryPoints; i++)
	{
		if (_isnan((double)m_pDistanceOutput[i]))
			continue;
		diff = m_pDistanceOutput[i] - distMean;
		var += diff*diff;
	}

	var /= (NNAReal)(nPoints - 1);
	distStd = sqrt(var);

	*pmean = distMean;
	*pstd = distStd;
	*pmin = distMin;
	*pmax = distMax;

	binSize = distMax * 1.001 / (float)nBins;
	rcpBinSize = 1.0 / binSize;
	
	for (i = 0; i < nBins; i++)
	{
		pBinDist[i] = binSize * (float)(i + 1);
		pBinCount[i] = 0;
	}
	for (i = 0; i < m_nQueryPoints; i++)
	{
		if (_isnan((double)m_pDistanceOutput[i]))
			continue;
		iBin = (int)(m_pDistanceOutput[i] * rcpBinSize);
		pBinCount[iBin]++;
	}
}

void CNearestNeighborAssistant::distanceStatsHistogram(int nBins, char *szFilename)
{
	int i;
	NNAReal distMean, distStd, distMin, distMax;
	NNAReal* pBinDist;
	int* pHisto;

	pBinDist = new NNAReal[nBins];
	_ASSERT(pBinDist != NULL);

	pHisto = new int[nBins];
	_ASSERT(pHisto != NULL);

	distanceStatsHistogram(&distMean, &distStd, &distMin, &distMax, nBins, pBinDist, pHisto);

	FILE* pf = fopen(szFilename, "w");
	if (pf == NULL)
		goto distanceStatsHistorgam_EXIT;

	fprintf(pf, "distances statistics:\n");
	fprintf(pf, "\tmean = %.10lf\n", distMean);
	fprintf(pf, "\tstd  = %.10lf\n", distStd);
	fprintf(pf, "\tmin  = %.10lf\n", distMin);
	fprintf(pf, "\tmax  = %.10lf\n", distMax);
	fprintf(pf, "\n");

	fprintf(pf, "distance histogram:\n");
	for (i = 0; i < nBins; i++)
	{
		fprintf(pf, "%.10lf\t%i\n", pBinDist[i], pHisto[i]);
	}
	fprintf(pf, "\n");

distanceStatsHistorgam_EXIT:
	if (pf != NULL)
		fclose(pf);
	delete [] pHisto;
	delete [] pBinDist;
}

