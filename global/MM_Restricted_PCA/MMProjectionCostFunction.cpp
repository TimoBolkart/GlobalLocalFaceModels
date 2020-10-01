#include "MMProjectionCostFunction.h"
#include "MultilinearModel.h"
#include "Definitions.h"
#include "KDTree3.h"

#include <iostream>

MMProjectionCostFunction::MMProjectionCostFunction(const Tensor* pMM, const std::vector<double>& dataMean, const std::vector<double>& targetPoints, const double maxPointDistance)
: vnl_cost_function(pMM->getModeDimension(2)+pMM->getModeDimension(3))
, m_pMultilinearModel(pMM)
, m_dataMean(dataMean)
, m_targetData(targetPoints)
, m_maxPointDistance(maxPointDistance)
{
	m_pKDTargetKDTree = new KDTree3(targetPoints);
}

MMProjectionCostFunction::~MMProjectionCostFunction()
{
	delete m_pKDTargetKDTree;
}

void MMProjectionCostFunction::compute(const vnl_vector<double>& x, double *f, vnl_vector<double>* g)
{
	const size_t d1(m_pMultilinearModel->getModeDimension(1));
	const size_t m2(m_pMultilinearModel->getModeDimension(2));
	const size_t m3(m_pMultilinearModel->getModeDimension(3));

	const double sqrMaxFrameDist = pow(m_maxPointDistance, 2);

	std::vector<double> w2;
	w2.resize(m2);

	for(size_t i = 0; i < m2; ++i)
	{
		w2[i] = x[i];
	}

	Tensor M3;
	m_pMultilinearModel->modeMultiply(w2, "T", m2, 1, 2, M3);

	std::vector<double> w3;
	w3.resize(m3);

	for(size_t i = 0; i < m3; ++i)
	{
		w3[i] = x[m2+i];
	}

	Tensor M2;
	m_pMultilinearModel->modeMultiply(w3, "T", m3, 1, 3, M2);

	Tensor currData;
	M3.modeMultiply(w3, "T", m3, 1, 3, currData);

	double E_data(0.0);
	std::vector<double> gradE_data;
	initVector(gradE_data, m2+m3, 0.0);

	std::vector<double> sourcePoints;
	sourcePoints.resize(d1);

	for(size_t i = 0; i < d1; ++i)
	{
		sourcePoints[i] = currData.getElement(i, 0, 0, 0)+m_dataMean[i];
	}

	std::vector<double> nnTargetPoints;
	std::vector<double> nnSqrPointDists;
	getNearestNeighbors(sourcePoints, nnTargetPoints, nnSqrPointDists);

	const size_t numPoints = d1/3;
	if((sourcePoints.size()!=nnTargetPoints.size())
		|| (nnTargetPoints.size() != d1)
		|| (nnSqrPointDists.size() != numPoints))
	{
		std::cout << "compute() - point set dimensions not correct" << std::endl;
		return;
	}

	std::vector<double> diffVec;
	diffVec.resize(d1);

	for(size_t i = 0; i < numPoints; ++i)
	{
		const double currSqrPointDist = nnSqrPointDists[i];
		E_data += std::min<double>(currSqrPointDist, sqrMaxFrameDist);

		const double diffx = sourcePoints[3*i]-nnTargetPoints[3*i];
		const double diffy = sourcePoints[3*i+1]-nnTargetPoints[3*i+1];
		const double diffz = sourcePoints[3*i+2]-nnTargetPoints[3*i+2];

		const double tmpPointScaling = currSqrPointDist < sqrMaxFrameDist ? 1.0 : m_maxPointDistance/sqrt(currSqrPointDist);
		diffVec[3*i] = tmpPointScaling*diffx;
		diffVec[3*i+1] = tmpPointScaling*diffy;
		diffVec[3*i+2] = tmpPointScaling*diffz;
	}

	// calc mode 2 gradient
	for(size_t i = 0; i < m2; ++i)
	{
		double out = 0.0;
		for(size_t j = 0; j < d1; ++j)
		{
			out += M2.getElement(j, i, 0, 0)*diffVec[j]; 
		}

		gradE_data[i] += 2.0*out;
	}

	// calc mode 3 gradient
	for(size_t i = 0; i < m3; ++i)
	{
		double out = 0.0;
		for(size_t j = 0; j < d1; ++j)
		{
			out += M3.getElement(j, 0, i, 0)*diffVec[j]; 
		}

		gradE_data[m2+i] += 2.0*out;
	}


	for(size_t i = 0; i < m2+m3; ++i)
	{
		(*g)[i] = gradE_data[i];
	}

	*f = E_data;
}
	
void MMProjectionCostFunction::getNearestNeighbors(const std::vector<double>& sourcePoints, std::vector<double>& nnTargetPoints, std::vector<double>& sqrPointDists)
{
	const size_t numSourcePoints = sourcePoints.size()/3;

	nnTargetPoints.clear();
	nnTargetPoints.resize(3*numSourcePoints);

	sqrPointDists.clear();
	sqrPointDists.resize(numSourcePoints);

	for(size_t i = 0; i < numSourcePoints; ++i)
	{
		std::vector<double> querySourcePoint;
		querySourcePoint.push_back(sourcePoints[3*i]);
		querySourcePoint.push_back(sourcePoints[3*i+1]);
		querySourcePoint.push_back(sourcePoints[3*i+2]);

		int nnPointIndex(0);
		double nnSqrPointDist(0);
		m_pKDTargetKDTree->getNearestPoint(querySourcePoint, nnPointIndex, nnSqrPointDist);

		nnTargetPoints[3*i] = m_targetData[3*nnPointIndex];
		nnTargetPoints[3*i+1] =m_targetData[3*nnPointIndex+1];
		nnTargetPoints[3*i+2] =m_targetData[3*nnPointIndex+2];

		sqrPointDists[i] = nnSqrPointDist;
	}
}