#include <random>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>


gsl_vector* generateRandomNumbers(unsigned long seed, size_t size) {
	gsl_vector* res = gsl_vector_alloc(size);
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i != size; i++) {
		gsl_vector_set(res, i, distribution(generator));
	}

	return res;
}


