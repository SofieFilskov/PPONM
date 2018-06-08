#include<stdio.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <gsl/gsl_multimin.h>

int rosen_func(const gsl_vector * current_point, void *params){
  double x = gsl_vector_get(current_point,0);
  double y = gsl_vector_get(current_point,1);

  return exp(-x)+exp(-y)-1-x*y;
}

int main() {
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;

  const size_t dim = 2;

  gsl_multimin_function F; //Giving the function
  F.f = rosen_func;
  F.n = dim;
  F.params = NULL;

  T = gsl_multimin_fminimizer_nmsimplex2;
  s = gsl_multimin_fminimizer_alloc(T, dim);

  int iter = 0;
  int status;

  /* Starting point */
  gsl_vector* start = gsl_vector_alloc (dim);
  gsl_vector_set (start, 0, 2.0);
  gsl_vector_set (start, 1, 2.0);

  /* Set initial step sizes to 0.1 */
  gsl_vector* step = gsl_vector_alloc (dim);
  gsl_vector_set_all (step, 0.1);

  gsl_multimin_fminimizer_set (s, &F, start, step);

  do{
    iter++;

    status = gsl_multimin_fminimizer_iterate(s);
    if (status)
      break;

    status = gsl_multimin_test_size(s->size,1e-4);
  }
  while (status == GSL_CONTINUE && iter < 1000);

  printf("For first, # of iterations %i\n", iter);

  gsl_vector_free(start);
  gsl_vector_free(step);
  gsl_multimin_fminimizer_free(s);

  return 0;
}
