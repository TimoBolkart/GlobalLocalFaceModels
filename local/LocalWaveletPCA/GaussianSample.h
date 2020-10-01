///////////////////////////////////////////////////////////////////////////////
//
//	GaussianSample.h
//
//	Header file for the CGaussianSample class
//
//	Alan Brunton 2008
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __GAUSSIANSAMPLE_H__
#define __GAUSSIANSAMPLE_H__


#include <math.h>

#include "RandomSample.h"

//+
class CGaussianSample
{
private:

	///////////////////////////////////////////////////////////////////////////
	//constructors/destructor
	///////////////////////////////////////////////////////////////////////////

	double									m_mean;
	double									m_var;
	double									m_std;

	CRandomSample							m_uniform;


public:

	///////////////////////////////////////////////////////////////////////////
	//constructors/destructor
	///////////////////////////////////////////////////////////////////////////
	//+
	CGaussianSample(double mean = 0.0, double var = 1.0): m_mean(mean), m_var(var), m_uniform(0)  
	{
		m_std = ::sqrt(m_var);
	}


	///////////////////////////////////////////////////////////////////////////
	//actions/operations
	///////////////////////////////////////////////////////////////////////////

	double getMean()						{ return m_mean; }
	double getVariance()					{ return m_var; }
	double getStandardDeviation()			{ return m_std; }
	//+
	void setMean(double mean)				{ m_mean = mean; }
	void setVariance(double var)
	{
		m_var = var;
		m_std = ::sqrt(m_var);
	}
	//+
	void setStandardDeviation(double std)
	{
		m_std = std;
		m_var = m_std*m_std;
	}


	///////////////////////////////////////////////////////////////////////////
	//actions/operations
	///////////////////////////////////////////////////////////////////////////

	void init(long seed)					{ m_uniform.reset(seed); }

	double next()
	{
		int i;
		double result = 0;

		//By central limit theorem use 48 = 12 * 4 trials of uniform random variable:
		for (i = 0; i < 48; i++)
			result += m_uniform.next();

		//the variable result is now of a Gaussian distribution (24, 4).
		result *= 0.5;

		//the variable result is now of a Gaussian distribution (12, 1).
		result -= 12.0;

		//scale the result by the standard deviation
		result *= m_std;

		//shift by the mean
		result += m_mean;

		return result;
	}
};


#endif //__GAUSSIANSAMPLE_H__


