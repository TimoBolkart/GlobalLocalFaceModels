#include "MultilinearModelHandler.h"	
#include "FileLoader.h"
#include "FileWriter.h"
#include "MathHelper.h"
#include "MMProjectionCostFunction.h"
#include "KDTree3.h"

#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgsb.h>


MultilinearModelHandler::MultilinearModelHandler()
{

}

MultilinearModelHandler::~MultilinearModelHandler()
{

}

void MultilinearModelHandler::clear()
{
	m_multilinearModel.clear();
	m_modeDimensions.clear();
	m_truncModeDimensions.clear();
	m_mean.clear();
	m_meanWeights.clear();
	m_neutralMeanWeights.clear();
}

void MultilinearModelHandler::reconstructForVariations(const std::vector<double>& variations, std::vector<double>& outPoints)
{
	if(m_truncModeDimensions.size()!=3)
	{
		std::cout << "Wrong number of mode diemensions " << m_truncModeDimensions.size() << " != " << 3 << std::endl;
		return;
	}

	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];
	if(variations.size() != m2+m3)
	{
		std::cout << "Wrong input variation vector dimension " << variations.size() << " != " << m2+m3 << std::endl;
		return;
	}

	std::vector<double> weightVector;
	for(size_t i = 0; i < m2; ++i)
	{
		double tmpWeight = getWeightForVariation(variations[i], m_meanWeights[i], 2);
		weightVector.push_back(tmpWeight);
	}

	for(size_t i = 0; i < m3; ++i)
	{
		double tmpWeight = getWeightForVariation(variations[m2+i], m_meanWeights[m2+i], 3);
		weightVector.push_back(tmpWeight);
	}

	reconstructForWeights(weightVector, outPoints);
}

void MultilinearModelHandler::reconstructForWeights(const std::vector<double>& weightVector, std::vector<double>& outPoints)
{
	if(m_truncModeDimensions.size()!=3)
	{
		std::cout << "Wrong number of mode diemensions " << m_truncModeDimensions.size() << " != " << 3 << std::endl;
		return;
	}

	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];
	if(weightVector.size() != m2+m3)
	{
		std::cout << "Wrong input weight vector dimension " << weightVector.size() << " != " << m2+m3 << std::endl;
		return;
	}

	std::vector<double> w2;
	w2.reserve(m_truncModeDimensions[1]);

	for(size_t i = 0; i < m_truncModeDimensions[1]; ++i)
	{
		w2.push_back(weightVector[i]);
	}

	Tensor t2;
	m_multilinearModel.modeMultiply(w2, "T", m_truncModeDimensions[1], 1, 2, t2);

	std::vector<double> w3;
	w3.reserve(m_truncModeDimensions[2]);

	for(size_t i = 0; i < m_truncModeDimensions[2]; ++i)
	{
		const size_t index = m_truncModeDimensions[1]+i;
		w3.push_back(weightVector[index]);
	}

	Tensor result;
	t2.modeMultiply(w3, "T", m_truncModeDimensions[2], 1, 3, result);

	outPoints.clear();
	outPoints.reserve(m_truncModeDimensions[0]);

	for(size_t i = 0; i < m_truncModeDimensions[0]; ++i)
	{
		const double value = result.getElement(i, 0, 0, 0);
		outPoints.push_back(value);
	}

	for(size_t i = 0; i < m_mean.size(); ++i)
	{
		outPoints[i] += m_mean[i];
	}
}

void MultilinearModelHandler::getDataMeanVertices(std::vector<double>& dataMean)
{
	reconstructForWeights(m_meanWeights, dataMean);
}

void MultilinearModelHandler::getNeutralDataMeanVertices(std::vector<double>& neutralDataMean)
{
	reconstructForWeights(m_neutralMeanWeights, neutralDataMean);
}

void MultilinearModelHandler::projectToMM(const std::vector<double>& targetPoints, const double maxPointDistance, std::vector<double>& outWeights)
{
	projectToMultilinearModelNonLinear(targetPoints, maxPointDistance, PROJECTION_MODE2_PRIOR_BOX_SIZE, PROJECTION_MODE3_PRIOR_BOX_SIZE, m_meanWeights, outWeights);
}

bool MultilinearModelHandler::importMultilinearModel(const std::string& sstrFileName)
{
	clear();

	std::vector<double> multModel;

	FileLoader loader;
	if(!loader.loadRestrictedMultilinearModel(sstrFileName, m_modeDimensions, m_truncModeDimensions, multModel, m_meanWeights, m_neutralMeanWeights, m_mean))
	{
		return false;
	}

	const size_t d1 = m_modeDimensions[0];
	const size_t d2 = m_modeDimensions[1];
	const size_t d3 = m_modeDimensions[2];

	const size_t m1 = m_truncModeDimensions[0];
	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];

	m_multilinearModel.init(multModel, d1, m2, m3);

	std::cout << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << "Data dimension" << std::endl;
	std::cout << "Mode1: " << d1 << std::endl;
	std::cout << "Mode2: " << d2 << std::endl;
	std::cout << "Mode3: " << d3 << std::endl;
	std::cout << std::endl;
	std::cout << "Multilinear Model dimension (truncated)" << std::endl;
	std::cout << "Mode1: " << m1 << std::endl;
	std::cout << "Mode2: " << m2 << std::endl;
	std::cout << "Mode3: " << m3 << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << std::endl;

	return true;
}

void MultilinearModelHandler::projectToMultilinearModelNonLinear(const std::vector<double>& targetPoints, const double maxPointDistance, const double maxMode2Variation, const double maxMode3Variation
																					, const std::vector<double>& initWeights, std::vector<double>& outWeights)
{
	const size_t d1(m_modeDimensions[0]);
	const size_t m2(m_truncModeDimensions[1]);
	const size_t m3(m_truncModeDimensions[2]);

	const size_t numWeights = m2+m3;
	if(initWeights.size() != numWeights)
	{
		std::cout << "MultilinearModelHandler::projectToMultilinearModelNonLinear(...) - Initial weights of wrong size " <<  initWeights.size() << " != " << numWeights << std::endl;
		return;
	}

	vnl_vector<long> boundSelection(numWeights, 2);

	vnl_vector<double> lowerBounds(numWeights, 0.0);
	vnl_vector<double> upperBounds(numWeights, 0.0);

	for(size_t i = 0; i < m2; ++i)
	{
		const double lowerBound = getWeightForVariation(-maxMode2Variation, m_meanWeights[i], 2);
		const double upperBound = getWeightForVariation(maxMode2Variation, m_meanWeights[i], 2);

		lowerBounds[i] = lowerBound;
		upperBounds[i] = upperBound;
	}

	for(size_t i = 0; i < m3; ++i)
	{
		const size_t sourceIndex = m2+i;

		const double lowerBound = getWeightForVariation(-maxMode3Variation, m_meanWeights[sourceIndex], 3);
		const double upperBound = getWeightForVariation(maxMode3Variation, m_meanWeights[sourceIndex], 3);

		lowerBounds[m2+i] = lowerBound;
		upperBounds[m2+i] = upperBound;
	}

	vnl_vector<double> x(numWeights, 0.0);
	for(size_t i = 0; i < numWeights; ++i)
	{
		x[i] = initWeights[i];
	}

	MMProjectionCostFunction fkt(&m_multilinearModel, m_mean, targetPoints, maxPointDistance);

	vnl_lbfgsb minimizer(fkt);
	minimizer.set_cost_function_convergence_factor(10000000);
	minimizer.set_projected_gradient_tolerance(0.00001);
	minimizer.set_max_function_evals(1000);
	minimizer.set_bound_selection(boundSelection);
	minimizer.set_lower_bound(lowerBounds);
	minimizer.set_upper_bound(upperBounds);
#ifdef TRACE_NONLINEAR_PROJECTION_METHOD
	minimizer.set_trace(true);
#endif
	minimizer.minimize(x);

	outWeights.clear();
	outWeights.resize(numWeights);

	for(size_t i = 0; i < numWeights; ++i)
	{
		outWeights[i] = x[i];
	}
}