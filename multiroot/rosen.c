#include<stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <math.h>

int rosen_func(const gsl_vector * current_point, void *params, gsl_vector * f){
  // f(x,y) = (1-x)^2 + 100*(y-x^2)^2
  //gradient f = (-2*(1-x)-400*x*(y-x^2)) i + 200*(y-x^2) j

  double x = gsl_vector_get(current_point,0);
  double y = gsl_vector_get(current_point,1);

  double dfdx = -2*(1-x)-400*x*(y-x*x);
  double dfdy = 200*(y-x*x);

  gsl_vector_set(f, 0, dfdx);
  gsl_vector_set(f, 1, dfdy);

  return GSL_SUCCESS;
}

int
print_state (size_t iter, gsl_multiroot_fsolver * s){
  printf ("iterations = %3u nxy = % .3f % .3f"
          "F(xy) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));

  return 0;
}

int main(){
  gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid; //This is type
  gsl_multiroot_fsolver *s; //This is the workspace

  const size_t dim = 2; //Dimension

  gsl_multiroot_function F; //Giving the function
  F.f = &rosen_func;
  F.n = dim;
  F.params = NULL;

  int status;
  size_t iter = 0;

  double xy_init[2] = {3.0, -2.0};
  gsl_vector * xy = gsl_vector_alloc(dim);
  gsl_vector_set(xy, 0, xy_init[0]);
  gsl_vector_set(xy, 1, xy_init[1]);

  s = gsl_multiroot_fsolver_alloc(T, dim);
  gsl_multiroot_fsolver_set(s, &F, xy);

  printf("Before solver:\n");
  print_state(iter,s);
  printf("\nDuring solver:\n");

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    print_state(iter, s);
    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, 1e-7);
  }
  while((status == GSL_CONTINUE) & (iter < 10000));

  printf("\nAfter solver:\n");
  print_state(iter, s);
  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (xy);

  return 0;
}
