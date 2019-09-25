#include <iostream>
#include <boost/timer.hpp>
#include <gsl/gsl_blas.h>
#include "MCFunctions_MeanReversion.h"


void quickTest() {
	unsigned long seed = 123;
	double kappa = 0.05, theta = 80, sigma = 0.15, deltaT = 1.0 / 252, startTime = 0, endTime = 300;
	double Y0 = 109;
	gsl_vector* Y = generateMCPath_MeanReversion(Y0, kappa, theta, sigma, startTime, endTime, deltaT, seed);
	
	gsl_vector* initialGuess = gsl_vector_alloc(3);
	gsl_vector_set(initialGuess, 0, 0.1);
	gsl_vector_set(initialGuess, 1, 82);
	gsl_vector_set(initialGuess, 2, 0.1);
	Inputs_MR* inputs_ptr = new Inputs_MR;
	inputs_ptr->deltaT = deltaT;
	inputs_ptr->Y = Y;
	search_minimum(inputs_ptr, my_func_MeanReversion, initialGuess);
	delete inputs_ptr;
}


int main()
{
	quickTest();
	return 0;
}