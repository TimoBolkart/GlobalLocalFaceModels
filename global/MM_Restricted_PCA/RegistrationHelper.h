#ifndef REGISTRATIONHELPER_H
#define REGISTRATIONHELPER_H

#include "MultilinearModel.h"

class RegistrationHelper
{
public:
	enum Resolution
	{
		HIGH_RESOLUTION,
		LOW_RESOLUTION
	};

	RegistrationHelper();

	~RegistrationHelper();

	void init(const std::vector<bool>& lowResValids);



	bool projectRegisteredToModel(const std::vector<double>& targetPoints, std::vector<double>& coefficients);

	bool projectRegisteredToModel(const std::vector<double>& targetPoints, const std::vector<bool>& targetValidValues, const std::vector<double>& initialCoefficients, std::vector<double> coefficients);



	void projectRegisteredSequenceToModel(const std::vector<double>& targetSequencePoints, );




	void projectToModel();

	void projectSequenceToModel();

private:

	Tensor m_lowResTensor;
};

#endif