#include <iostream>
#include <boost/timer.hpp>
#include <gsl/gsl_blas.h>
#include "MCFunctions_MeanReversion.h"

void quickTest() {
	unsigned long seed = 123;
	double kappa = 0.05, theta = 80, sigma = 0.15, deltaT = 1.0 / 252, startTime = 0, endTime = 300;
	double Y0 = 109;
	gsl_vector* Y = generateMCPath_MeanReversion(Y0, kappa, theta, sigma, startTime, endTime, deltaT, seed);
	
}


int main()
{
	unsigned long seed = 123;

	return 0;
}
