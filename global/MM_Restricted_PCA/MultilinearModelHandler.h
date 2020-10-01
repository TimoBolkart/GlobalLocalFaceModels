#ifndef MULTILINEARMODELHANLDER
#define MULTILINEARMODELHANLDER

#include "Definitions.h"
#include "DataContainer.h"
#include "MultilinearModel.h"

#include <vector>
#include <map>

class KDTree3;

class MultilinearModelHandler
{
public:
	//! Handler for all multilinear model related stuff computation stuff. 
	MultilinearModelHandler();

	~MultilinearModelHandler();

	//! Clear all data structures.
	void clear();

	//! Compute the vertices for variations. The variation parameters are centered at 0 and and scaled. 
	//! E.g. To get the variations where some principal components are r*stdev far from the mean, call this function with a vector, where corresponding entries are r.
	//! \param variations		(m2+m3)-dimensional input vector where the first m2 entries represent the variations of the second mode, and the following m3 entries the variation of third mode
	//! \param outPoints			vertices of the reconstruction
	void reconstructForVariations(const std::vector<double>& variations, std::vector<double>& outPoints);

	//! Compute the vertices for mode coefficients.
	//! \param weightVector		(m2+m3)-dimensional input vector where the first m2 entries represent the coefficients of the second mode, and the following m3 entries the cofficients of third mode
	//! \param outPoints			vertices of the reconstruction
	void reconstructForWeights(const std::vector<double>& weightVector, std::vector<double>& outPoints);

	//! Get the vertices of the learned data mean.
	//! \param dataMean			learned data mean
	void getDataMeanVertices(std::vector<double>& dataMean);

	//! Get the vertices of the learned neutral data mean.
	//! \param neutralDataMean	learned neutral data mean
	void getNeutralDataMeanVertices(std::vector<double>& neutralDataMean);

	//! Project data to the multilinear model space.
	//! \param targetPoints				target vertices
	//! \param maxPointDistance		maximum distance of corresponding points
	//! \param outWeight					coefficients fo the projection				
	void projectToMM(const std::vector<double>& targetPoints, const double maxPointDistance, std::vector<double>& outWeights);

	//! Import a multilinear model.
	//! \param sstrFileName			full file name of the model
	//! \return true if import was successful
	bool importMultilinearModel(const std::string& sstrFileName);

private:

	//! Convert from centered and scaled variation parameter to model coefficient.
	//! \param variation			parameter in the centered and scaled variation space
	//! \param meanWeight		mean weight of the model of corresponding coefficient
	//! \param mode				mode of the model
	//! \return	the converted model coefficient
	double getWeightForVariation(const double variation, const double meanWeight, const size_t mode)
	{
		if((mode<1) || (m_modeDimensions.size() < mode-1))
		{
			return 0.0;
		}

		const double tmp = sqrt(static_cast<double>(m_modeDimensions[mode-1]));
		return variation/tmp+meanWeight;
	}

	//! Convert from the model coefficient to the centered and scaled variation parameter.
	//! \param weight				model coefficient
	//! \param meanWeight		mean weight of the model of corresponding coefficient
	//! \param mode				mode of the model
	//! \return	the converted parameter in the centered and scaled variation space
	double getVariationForWeight(const double weight, const double meanWeight, const size_t mode)
	{
		if((mode<1) || (m_modeDimensions.size() < mode-1))
		{
			return 0.0;
		}

		const double tmp = sqrt(static_cast<double>(m_modeDimensions[mode-1]));
		return (weight-meanWeight)*tmp;
	}

	//! Helper function to project data to the multilinear model space.
	//! \param targetPoints					target vertices
	//! \param maxPointDistance			maximum distance of corresponding points
	//! \param maxMode2Variation			size of the prior-box of the second mode (maximum distance from the model mode mean of the second mode)
	//! \param maxMode3Variaton			size of the prior-box of the third mode (maximum distance from the model mode mean of the third mode)
	//! \param initWeight					(m2+m3)-dimensional vector with initial coefficients in multilinear model space
	//! \param outWeight						coefficients fo the projection
	void projectToMultilinearModelNonLinear(const std::vector<double>& targetPoints, const double maxPointDistance, const double maxMode2Variation, const double maxMode3Variation, const std::vector<double>& initWeights, std::vector<double>& outWeights);

	//! Multilinear model tensor of dimension d1 x m2 x m3
	Tensor m_multilinearModel;
	//! Dimension of the original data d1 d2 d3
	std::vector<size_t> m_modeDimensions;
	//! Dimensions of the truncated model d1 m2 m3
	std::vector<size_t> m_truncModeDimensions;
	//! Learned data mean
	std::vector<double> m_mean;
	//! (m2+m3)-dimensional vector with mean coefficients of mode 2 and mode 3
	std::vector<double> m_meanWeights;
	//! (m2+m3)-dimensional vector with mean coefficients of mode 2 and netural mean coefficients of mode 3
	std::vector<double> m_neutralMeanWeights;
};
#endif