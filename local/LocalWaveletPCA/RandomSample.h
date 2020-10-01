///////////////////////////////////////////////////////////////////////////////
//
//	RandomSample.h
//
//	Header file for the CRandomSample class
//	represents a randomly distributed sample from 0 to 1
//
//	Alan Brunton 2008
//
///////////////////////////////////////////////////////////////////////////////


#ifndef __AB_RANDOMSAMPLE_H__
#define __AB_RANDOMSAMPLE_H__


#define IA		16807
#define	IM		2147483647
#define AM		(1.0/IM)
#define IQ		127773
#define IR		2836
#define NTAB	32
#define NDIV	(1+(IM-1)/NTAB)
#define EPS		1.2e-7
#define RNMX	(1.0-EPS)


//Random number generator of Park and Miller with Bays-Durham shuffle
//Numerical Recipes in C (1992 Cambridge University Press) Chapter 7.1
class CRandomSample
{
private:

	long		idum;
	long		iy;
	long		iv[NTAB];

	void init()
	{
		int j;
		long k;

		if (idum < 0)
			idum = -idum;
		else if (idum == 0)
			idum = 1;
		for (j = NTAB + 7; j >= 0; j--)
		{
			k = idum / IQ;
			idum = IA * (idum - k * IQ) - IR * k;
			if (idum < 0)
				idum += IM;
			if (j < NTAB)
				iv[j] = idum;
		}
		iy = iv[0];
	}

public:

	CRandomSample(long seed)
	{
		reset(seed);
	}

	void reset(long seed)
	{
		idum = seed;
		iy = 0;
		init();
	}

	float next()
	{
		int j;
		long k;
		float temp;

		k = idum / IQ;
		idum = IA * (idum - k * IQ) - IR * k;
		if (idum < 0)
			idum += IM;
		j = iy / NDIV;
		iy = iv[j];
		iv[j] = idum;
		temp = AM * iy;
		if (temp > RNMX)
			return RNMX;
		return temp;
	}
};


#endif //__AB_RANDOMSAMPLE_H__


