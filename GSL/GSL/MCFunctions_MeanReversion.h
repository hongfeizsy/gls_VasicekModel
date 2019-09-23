#include <random>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include "MCFunctions.h"

double YtFunction_MeanReversion(double Yt, double kappa, double theta, double sigma, double deltaT, double BM) {
	return Yt * (1 - kappa * deltaT) + kappa * theta + sigma * sqrt(Yt * deltaT) * BM;
}


gsl_vector* generateMCPath_MeanReversion(double Y0, double kappa, double theta, double sigma, double startTime, double endTime, double deltaT, unsigned long seed) {
	int size = round((endTime - startTime) / deltaT);
	gsl_vector* BMs = generateRandomNumbers(seed, size);
	gsl_vector* res = gsl_vector_alloc(size + 1);
	gsl_vector_set(res, 0, Y0);
	for (int i = 1; i != res->size; i++) {
		gsl_vector_set(res, i, YtFunction_MeanReversion(gsl_vector_get(res, i - 1), kappa, theta, sigma, deltaT, gsl_vector_get(BMs, i - 1)));
	}

	return res;
}


