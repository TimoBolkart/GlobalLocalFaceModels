#ifndef MMPROJECTIONCOSTFUNCTION_H
#define MMPROJECTIONCOSTFUNCTION_H

#include <vnl/vnl_vector.h>
#include <vnl/vnl_cost_function.h>

#include <vector>

class Tensor;
class KDTree3;

class MMProjectionCostFunction : public vnl_cost_function
{
public:
	//! Specific projection cost function.
	//! \param pMM								pointer to the multilinear model used for error minimization
	//! \param dataMean						learned data mean
	//! \param targetPoints					
	//! \param maxPointDistance			
	MMProjectionCostFunction(const Tensor* pMM, const std::vector<double>& dataMean, const std::vector<double>& targetPoints, const double maxPointDistance);

	~MMProjectionCostFunction();

	//! Derived function to compute function value and gradient.
	//! \param x		initial values for current iteration
	//! \param f		function value evaluated on x
	//! \param g		gradient on x
	virtual void compute(const vnl_vector<double>& x, double *f, vnl_vector<double>* g);

private:
	//! Computes nearest neighbors for a set of vertices.
	//! \param sourcePoints					set of vertices the nearest neighbors are computed for
	//! \param nnTargetPoints				set fo nearest neighbors
	//! \param sqrPointDists				set of squared Euclidean distances of the nearest neighbors
	void getNearestNeighbors(const std::vector<double>& sourcePoints, std::vector<double>& nnTargetPoints, std::vector<double>& sqrPointDists);

	template<typename T>
	void initVector(std::vector<T>& vec, const size_t dim, const T& value)
	{
		vec.clear();
		vec.resize(dim);

		for(size_t i = 0; i < dim; ++i)
		{
			vec[i] = value;
		}
	}

	const Tensor* m_pMultilinearModel;
	std::vector<double> m_dataMean;
	std::vector<double> m_targetData;
	KDTree3* m_pKDTargetKDTree;
	double m_maxPointDistance;
};

#endif