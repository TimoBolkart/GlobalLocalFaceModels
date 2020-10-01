///////////////////////////////////////////////////////////////////////////////
//
//	TemplateRegistration.h
//
//	Header file for the CTemplateRegistration class
//	abstract base class
//
//	Alan Brunton, NRC-IIT VIT January 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __TEMPLATEREGISTRATION_H__
#define __TEMPLATEREGISTRATION_H__


#include "StereoFace.h"


class CTemplateRegistration
{
protected:

	double*									m_pTemplateVertices;
	int										m_nTemplateVertices;

	int										m_nLandmarks;
	int*									m_piTemplateLandmarks;
	double*									m_pTemplateLandmarks;
	double*									m_pTargetLandmarks;


public:

	CTemplateRegistration()
	{
		m_pTemplateVertices = NULL;
		m_nTemplateVertices = 0;
		m_nLandmarks = 0;
		m_piTemplateLandmarks = NULL;
		m_pTemplateLandmarks = NULL;
		m_pTargetLandmarks = NULL;
	}

	virtual void init() = 0;

	virtual void term()
	{
		if (m_piTemplateLandmarks != NULL)
		{
			delete [] m_piTemplateLandmarks;
			m_piTemplateLandmarks = NULL;
		}
		if (m_pTemplateLandmarks != NULL)
		{
			delete [] m_pTemplateLandmarks;
			m_pTemplateLandmarks = NULL;
		}
		if (m_pTargetLandmarks != NULL)
		{
			delete [] m_pTargetLandmarks;
			m_pTargetLandmarks = NULL;
		}
	}

	void setTemplateVertices(double* pVertices, int nVertices)
	{
		m_pTemplateVertices = pVertices;
		m_nTemplateVertices = nVertices;
	}

	void setLandmarks(int nLandmarks, double* pTemplateLandmarks, double* pTargetLandmarks)
	{
		if (m_pTemplateLandmarks != NULL)
			delete [] m_pTemplateLandmarks;
		if (m_pTargetLandmarks != NULL)
			delete [] m_pTargetLandmarks;

		m_nLandmarks = nLandmarks;
		int nAlloc = m_nLandmarks * 3;
		m_pTemplateLandmarks = new double[nAlloc];
		_ASSERT(m_pTemplateLandmarks != NULL);
		memcpy(m_pTemplateLandmarks, pTemplateLandmarks, nAlloc * sizeof(double));

		m_pTargetLandmarks = new double[nAlloc];
		_ASSERT(m_pTargetLandmarks != NULL);
		memcpy(m_pTargetLandmarks, pTargetLandmarks, nAlloc * sizeof(double));
	}

	void setLandmarks(int nLandmarks, int* piTemplateLandmarks, double* pTargetLandmarks)
	{
		int i, j, k;

		_ASSERT(m_pTemplateVertices != NULL);

		if (m_piTemplateLandmarks != NULL)
			delete [] m_piTemplateLandmarks;
		if (m_pTemplateLandmarks != NULL)
			delete [] m_pTemplateLandmarks;
		if (m_pTargetLandmarks != NULL)
			delete [] m_pTargetLandmarks;

		m_nLandmarks = nLandmarks;

		m_piTemplateLandmarks = new int[m_nLandmarks];
		_ASSERT(m_piTemplateLandmarks != NULL);
		memcpy(m_piTemplateLandmarks, piTemplateLandmarks, m_nLandmarks * sizeof(int));

		int nAlloc = m_nLandmarks * 3;
		m_pTemplateLandmarks = new double[nAlloc];
		_ASSERT(m_pTemplateLandmarks != NULL);

		j = 0;
		for (i = 0; i < m_nLandmarks; i++)
		{
			k = m_piTemplateLandmarks[i] * 3;
			m_pTemplateLandmarks[j]		= m_pTemplateVertices[k];
			m_pTemplateLandmarks[j + 1]	= m_pTemplateVertices[k + 1];
			m_pTemplateLandmarks[j + 2]	= m_pTemplateVertices[k + 2];
			j += 3;
		}

		m_pTargetLandmarks = new double[nAlloc];
		_ASSERT(m_pTargetLandmarks != NULL);
		memcpy(m_pTargetLandmarks, pTargetLandmarks, nAlloc * sizeof(double));
	}

	virtual void registerTemplate() = 0;

	virtual int getNumRegistrationParams() = 0;
	virtual void getRegistrationParams(double* pParams) = 0;

	virtual void getRegisteredVertices(double* pVertices) = 0;
};


#endif //__TEMPLATEREGISTRATION_H__


