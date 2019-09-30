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


double YtFunction_MeanReversion(double Yt, double kappa, double theta, double sigma, double deltaT, double BM) {
	return Yt * (1 - kappa * deltaT) + kappa * theta * deltaT + sigma * sqrt(deltaT) * BM;
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


struct Inputs_MR
{
	gsl_vector* Y; double deltaT;
};


template<class myStruct>
void search_minimum(myStruct params, double(*my_func)(const gsl_vector*, void*), const gsl_vector* initialGuess) {
	/*'myStruct' has to be defined in advance.*/
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc(initialGuess->size);
	gsl_vector_memcpy(x, initialGuess);

	/* Set initial step size to 1 */
	ss = gsl_vector_alloc(initialGuess->size);
	gsl_vector_set_all(ss, 1.0);

	/* Initialize method and iterate */
	minex_func.n = (initialGuess->size);
	minex_func.f = my_func;
	minex_func.params = params;

	s = gsl_multimin_fminimizer_alloc(T, initialGuess->size);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-2);

		if (status == GSL_SUCCESS) {
			printf("Converged to minimum at\n");
		}
		printf("%10.3e, %10.3e, %10.3e, f() = %7.3f, size = %.3f, iter = %d\n", gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), \
			gsl_vector_get(s->x, 2), s->fval, size, iter);
	} while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);
}


double logDensityFunc_MeanReversion(double ytp1, double yt, double deltaT, double kappa, double theta, double sigma) {
	return -(ytp1 - yt * (1 - kappa * deltaT) - kappa * theta * deltaT) * (ytp1 - yt * (1 - kappa * deltaT) - kappa * theta * deltaT) \
		/ (2 * sigma * sigma * deltaT) - 0.5 * log(2 * M_PI) - log(sigma * deltaT);
}


double generateLikelihoodFunction_MeanReversion(const gsl_vector* Y, const double& deltaT, const double& kappa, const double& theta, const double& sigma) {
	double res = 0;
	for (int i = 0; i != Y->size - 1; i++) {
		res += logDensityFunc_MeanReversion(gsl_vector_get(Y, i + 1), gsl_vector_get(Y, i), deltaT, kappa, theta, sigma);
	}

	return res;
}


double my_func_MeanReversion(const gsl_vector* variables, void* params) {
	double x, y, z, deltaT;
	Inputs_MR* inputs = (Inputs_MR*)params;

	x = gsl_vector_get(variables, 0);
	y = gsl_vector_get(variables, 1);
	z = gsl_vector_get(variables, 2);
	deltaT = inputs->deltaT;
	gsl_vector* Y = inputs->Y;

	return -generateLikelihoodFunction_MeanReversion(Y, deltaT, x, y, z);
}

