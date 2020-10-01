///////////////////////////////////////////////////////////////////////////////
//
//	WaveletShapeSampler.cpp
//
//	Source file for the CWaveletShapeSampler class and helper classes
//
//	Alan Brunton
//
///////////////////////////////////////////////////////////////////////////////


#include "WaveletShapeSampler.h"


///////////////////////////////////////////////////////////////////////////////
//member functions
///////////////////////////////////////////////////////////////////////////////

void CWaveletShapeSampler::initActiveWaveletModel()
{
	_ASSERT(m_pActiveWaveletModel == NULL);
	_ASSERT(m_pIntermediateModel == NULL);
	_ASSERT(m_pReconVerts == NULL);

	m_pActiveWaveletModel = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, m_priors[0]->m_pMean->getBaseResWidth(), m_priors[0]->m_pMean->getBaseResHeight(), m_priors[0]->m_pMean->getNumLevels());
	_ASSERT(m_pActiveWaveletModel != NULL);
	m_pIntermediateModel = (CBSplineGridWavelet<C3Vectorf>*)CBSGWFactory::createBSGW(eBSGWFloat3, m_priors[0]->m_pMean->getBaseResWidth(), m_priors[0]->m_pMean->getBaseResHeight(), m_priors[0]->m_pMean->getNumLevels());
	_ASSERT(m_pIntermediateModel != NULL);
	m_pReconVerts = new C3Vectorf[m_pActiveWaveletModel->getNumCoefficients()];
	_ASSERT(m_pReconVerts != NULL);

	m_pIntermediateModel->copy(m_priors[0]->m_pMean);
	reconstructGeometry();
}

void CWaveletShapeSampler::clear()
{
	delete m_pActiveWaveletModel;
	m_pActiveWaveletModel = NULL;
	delete m_pIntermediateModel;
	m_pIntermediateModel = NULL;
	delete m_pReconVerts;
	m_pReconVerts = NULL;
}

//Setting the prior cost to zero avoids the problem of the model "pulling" towards the mean

void CWaveletShapeSampler::optimizeActiveWaveletModel()
{
	_ASSERT(m_priors.size() > 0);
	_ASSERT(m_pActiveWaveletModel != NULL);

#if WSS_PROFILE
	m_profOptimize.start();
#endif

	int i, j, k;
	int iPrior, iObs;
	int nCoeff = m_pActiveWaveletModel->getLevelHeight(m_level) * m_pActiveWaveletModel->getLevelWidth(m_level);

	float dataCost, priorCost, priorCostMin, priorCostTotal, totalCost, totalCostMin, sampleValMin;
	float primarySigma;
	C3Vectorf vb, vc, vstd;

	float* pSamples = new float[m_nSamples];
	_ASSERT(pSamples != NULL);
	float* pTranSamples = new float[m_nSamples];
	_ASSERT(pTranSamples != NULL);
	float* pPriorCost = new float[m_nSamples];
	_ASSERT(pPriorCost != NULL);
	float* pSecondaryPrior = new float[m_nSamples];
	_ASSERT(pSecondaryPrior != NULL);

	priorCostTotal = 0.f;
	totalCostMin = 0.f;

	printf("updating observations...\n");
//	m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
	reconstructGeometry();
	for (iObs = 0; iObs < m_observations.size(); iObs++)
		m_observations[iObs]->updateObserveration();

	for (i = m_iStartCoefficient; i < nCoeff; i++)
	{
		printf("coefficient %i: total prior energy = %f; total energy = %f\n", i, priorCostTotal, totalCostMin);

#if WSS_PROFILE
		m_profCoefficient.start();
#endif //WSS_PROFILE

		j = m_pActiveWaveletModel->getSerialIndex(i);

		//get current value for this coefficient
		vb = m_pActiveWaveletModel->getCoefficient(j);
		
		//get primary standard deviation for this coefficient
		vstd = m_priors[0]->m_pStdDev->getCoefficient(j);
		printf("%.10f %.10f %.10f\n", vstd.x, vstd.y, vstd.z);

		//first local principal component
		primarySigma = vstd.x;
		if (primarySigma > (float)1e-12)
		{
			//compute initial cost before modifying model parameter
			pSamples[0] = sampleValMin = vb.x;
			m_priors[0]->evaluateSamples(j, 0, 1, pSamples, pPriorCost);
			for (iPrior = 1; iPrior < m_priors.size(); iPrior++)
			{
				pTranSamples[0] = pSamples[0];
				m_priors[iPrior]->forwardTransformSamples(j, 0, 1, pTranSamples);
				m_priors[iPrior]->evaluateSamples(j, 0, 1, pTranSamples, pSecondaryPrior, pPriorCost);
			}
			priorCost = 0; //pPriorCost[0];
			m_priors[0]->inverseTransformCoefficient(j, vb, vc);
			m_pIntermediateModel->setCoefficient(j, vc);
//			m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
			reconstructGeometry();
			dataCost = 0.f;
			for (iObs = 0; iObs < m_observations.size(); iObs++)
				dataCost += m_observations[iObs]->computeDataCost() * m_observations[iObs]->m_weight;
			totalCostMin = priorCost + dataCost;
			priorCostMin = priorCost;

#if WSS_PROFILE
			m_profPriors.start();
#endif //WSS_PROFILE

			//sample from primary prior
			switch (m_eSamplingType)
			{
			case eWSSUniform:
				m_priors[0]->drawNSamplesUniform(j, 0, m_nSamples, pSamples, pPriorCost);
				break;
			case eWSSRandom:
				m_priors[0]->drawNSamplesGaussian(j, 0, m_nSamples, pSamples, pPriorCost);
				break;
			}

			//compute priors
			for (iPrior = 1; iPrior < m_priors.size(); iPrior++)
			{
				memcpy(pTranSamples, pSamples, m_nSamples * sizeof(float));
				m_priors[iPrior]->forwardTransformSamples(j, 0, m_nSamples, pTranSamples);
				m_priors[iPrior]->evaluateSamples(j, 0, m_nSamples, pTranSamples, pSecondaryPrior, pPriorCost);
			}

#if WSS_PROFILE
			m_profPriors.stop();
#endif //WSS_PROFILE

			//loop over samples
			for (k = 0; k < m_nSamples; k++)
			{
				vb.x = pSamples[k];
				m_priors[0]->inverseTransformCoefficient(j, vb, vc);
				m_pIntermediateModel->setCoefficient(j, vc);

#if WSS_PROFILE
				m_profWaveletRecon.start();
#endif WSS_PROFILE
				//reconstruct surface
//				m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
				reconstructGeometry();
#if WSS_PROFILE
				m_profWaveletRecon.stop();
#endif WSS_PROFILE

				//compute observations
				dataCost = 0.f;
				for (iObs = 0; iObs < m_observations.size(); iObs++)
					dataCost += m_observations[iObs]->computeDataCost() * m_observations[iObs]->m_weight;

				priorCost = 0; //pPriorCost[k];
				totalCost = priorCost + dataCost;
				if (totalCost < totalCostMin)
				{
					sampleValMin = pSamples[k];
					totalCostMin = totalCost;
					priorCostMin = priorCost;
				}
			}

			vb.x = sampleValMin;
			priorCostTotal += priorCostMin;
		}

		//second local principal component
		primarySigma = vstd.y;
		if (primarySigma > (float)1e-12)
		{
			//compute initial cost before modifying model parameter
			pSamples[0] = sampleValMin = vb.y;
			m_priors[0]->evaluateSamples(j, 1, 1, pSamples, pPriorCost);
			for (iPrior = 1; iPrior < m_priors.size(); iPrior++)
			{
				pTranSamples[0] = pSamples[0];
				m_priors[iPrior]->forwardTransformSamples(j, 1, 1, pTranSamples);
				m_priors[iPrior]->evaluateSamples(j, 1, 1, pTranSamples, pSecondaryPrior, pPriorCost);
			}
			priorCost = 0; //pPriorCost[0];
			m_priors[0]->inverseTransformCoefficient(j, vb, vc);
			m_pIntermediateModel->setCoefficient(j, vc);
//			m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
			reconstructGeometry();
			dataCost = 0.f;
			for (iObs = 0; iObs < m_observations.size(); iObs++)
				dataCost += m_observations[iObs]->computeDataCost() * m_observations[iObs]->m_weight;
			totalCostMin = priorCost + dataCost;
			priorCostMin = priorCost;

#if WSS_PROFILE
			m_profPriors.start();
#endif //WSS_PROFILE

			//sample from primary prior
			switch (m_eSamplingType)
			{
			case eWSSUniform:
				m_priors[0]->drawNSamplesUniform(j, 1, m_nSamples, pSamples, pPriorCost);
				break;
			case eWSSRandom:
				m_priors[0]->drawNSamplesGaussian(j, 1, m_nSamples, pSamples, pPriorCost);
				break;
			}

			//compute priors
			for (iPrior = 1; iPrior < m_priors.size(); iPrior++)
			{
				memcpy(pTranSamples, pSamples, m_nSamples * sizeof(float));
				m_priors[iPrior]->forwardTransformSamples(j, 1, m_nSamples, pTranSamples);
				m_priors[iPrior]->evaluateSamples(j, 1, m_nSamples, pTranSamples, pSecondaryPrior, pPriorCost);
			}

#if WSS_PROFILE
			m_profPriors.stop();
#endif

			//loop over samples
			for (k = 0; k < m_nSamples; k++)
			{
				vb.y = pSamples[k];
				m_priors[0]->inverseTransformCoefficient(j, vb, vc);
				m_pIntermediateModel->setCoefficient(j, vc);

#if WSS_PROFILE
				m_profWaveletRecon.start();
#endif WSS_PROFILE
				//reconstruct surface
//				m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
				reconstructGeometry();
#if WSS_PROFILE
				m_profWaveletRecon.stop();
#endif WSS_PROFILE

				//compute observations
				dataCost = 0.f;
				for (iObs = 0; iObs < m_observations.size(); iObs++)
					dataCost += m_observations[iObs]->computeDataCost() * m_observations[iObs]->m_weight;

				priorCost = 0; //pPriorCost[k];
				totalCost = priorCost + dataCost;
				if (totalCost < totalCostMin)
				{
					sampleValMin = pSamples[k];
					totalCostMin = totalCost;
					priorCostMin = priorCost;
				}
			}

			vb.y = sampleValMin;
			priorCostTotal += priorCostMin;
		}

		//third local principal component
		primarySigma = vstd.z;
		if (primarySigma > (float)1e-12)
		{
			//compute initial cost before modifying model parameter
			pSamples[0] = sampleValMin = vb.z;
			m_priors[0]->evaluateSamples(j, 2, 1, pSamples, pPriorCost);
			for (iPrior = 1; iPrior < m_priors.size(); iPrior++)
			{
				pTranSamples[0] = pSamples[0];
				m_priors[iPrior]->forwardTransformSamples(j, 2, 1, pTranSamples);
				m_priors[iPrior]->evaluateSamples(j, 2, 1, pTranSamples, pSecondaryPrior, pPriorCost);
			}
			priorCost = 0; //pPriorCost[0];
			m_priors[0]->inverseTransformCoefficient(j, vb, vc);
			m_pIntermediateModel->setCoefficient(j, vc);
//			m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
			reconstructGeometry();
			dataCost = 0.f;
			for (iObs = 0; iObs < m_observations.size(); iObs++)
				dataCost += m_observations[iObs]->computeDataCost() * m_observations[iObs]->m_weight;
			totalCostMin = priorCost + dataCost;
			priorCostMin = priorCost;

#if WSS_PROFILE
			m_profPriors.start();
#endif

			//sample from primary prior
			switch (m_eSamplingType)
			{
			case eWSSUniform:
				m_priors[0]->drawNSamplesUniform(j, 2, m_nSamples, pSamples, pPriorCost);
				break;
			case eWSSRandom:
				m_priors[0]->drawNSamplesGaussian(j, 2, m_nSamples, pSamples, pPriorCost);
				break;
			}

			//compute priors
			for (iPrior = 1; iPrior < m_priors.size(); iPrior++)
			{
				memcpy(pTranSamples, pSamples, m_nSamples * sizeof(float));
				m_priors[iPrior]->forwardTransformSamples(j, 2, m_nSamples, pTranSamples);
				m_priors[iPrior]->evaluateSamples(j, 2, m_nSamples, pTranSamples, pSecondaryPrior, pPriorCost);
			}

#if WSS_PROFILE
			m_profPriors.stop();
#endif

			//loop over samples
			for (k = 0; k < m_nSamples; k++)
			{
				vb.z = pSamples[k];
				m_priors[0]->inverseTransformCoefficient(j, vb, vc);
				m_pIntermediateModel->setCoefficient(j, vc);

#if WSS_PROFILE
				m_profWaveletRecon.start();
#endif WSS_PROFILE
				//reconstruct surface
//				m_pIntermediateModel->reconstructCopy(m_pIntermediateModel->getNumLevels() - 1, m_pReconVerts, m_pIntermediateModel->getFullResWidth(), 1);
				reconstructGeometry();
#if WSS_PROFILE
				m_profWaveletRecon.stop();
#endif WSS_PROFILE

				//compute observations
				dataCost = 0.f;
				for (iObs = 0; iObs < m_observations.size(); iObs++)
					dataCost += m_observations[iObs]->computeDataCost() * m_observations[iObs]->m_weight;

				priorCost = 0; //pPriorCost[k];
				totalCost = priorCost + dataCost;
				if (totalCost < totalCostMin)
				{
					sampleValMin = pSamples[k];
					totalCostMin = totalCost;
					priorCostMin = priorCost;
				}
			}

			vb.z = sampleValMin;
			priorCostTotal += priorCostMin;
		}

		//save coefficients
		m_pActiveWaveletModel->setCoefficient(j, vb);
		m_priors[0]->inverseTransformCoefficient(j, vb, vc);
		m_pIntermediateModel->setCoefficient(j, vc);

#if WSS_PROFILE
		m_profCoefficient.stop();
#endif //WSS_PROFILE
	}
	
	m_pIntermediateModel->reconstruct(m_level);

	m_iStartCoefficient = nCoeff;
	m_level++;

	delete [] pSamples;
	delete [] pPriorCost;
	delete [] pTranSamples;
	delete [] pSecondaryPrior;

#if WSS_PROFILE
	m_profOptimize.stop();
#endif //WSS_PROFILE
}