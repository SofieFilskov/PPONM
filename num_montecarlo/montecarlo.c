#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

void plainmc(int dim, double *a, double *b, double f(double* x), int N,
double* result, double* error){
  double V = 1;
  for (int i = 0; i < dim; i++) {
    V *= b[i]-a[i];
  }

  double sum = 0;
  double sum2 = 0;
  double fx;
  double x[dim];

  const gsl_rng_type* T;
  gsl_rng* r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (int i = 0; i < N; i++){
    for (int j = 0; j < dim; j++){
      double u = gsl_rng_uniform (r);
      x[j] = a[j] + u*(b[j]-a[j]);
    }
    fx = f(x);
    sum += fx;
    sum2 += fx*fx;
  }
  gsl_rng_free (r);

  double av_f = sum/N; /* <f> */
  double sigma = sum2/N - av_f*av_f; /*sigma^2 = <f^2> - <f>^2 */
  *result = av_f*V;
  *error = sqrt(sigma/N)*V;
}
