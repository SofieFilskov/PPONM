#include<stdio.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_odeiv2.h>
#include<math.h>
#include<gsl/gsl_errno.h>
#include<assert.h>

int hydro_func(double r, const double f[], double dfdr[], void* params){
  double eps = *(double* )params;
  dfdr[0] = f[1]; //f' = v
  dfdr[1] = -2*(eps*f[0] + f[0]/r); //v' = f'' = -2 eps f - 2/r f
  return GSL_SUCCESS;
}


double ode(double r, double eps){
  const double rmin = 1e-3;
  if (r<rmin)
    return r-r*r;

  //making an ODE system
  //double eps = 0.01;
  gsl_odeiv2_system sys;
  sys.function = hydro_func;
  sys.jacobian = NULL;
  sys.dimension = 2;
  sys.params =  &eps;

  //Make it ode45
  gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;

  //The initial step size. We do so the step is negative if x is negative
  double hstart = copysign(0.1, r);

  //Absolute and relative error
  double e_abs = 1e-7;
  double e_rel = 1e-7;

  //Makes the driver
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, \
    T, hstart, e_abs, e_rel);

  //Start conditions for applying the driver
  //double rmin = 1e-3;
  double r0 = rmin;
  double f[2] = {r0-r0*r0, 1-2*r0};
  gsl_odeiv2_driver_apply(driver, &r0, r, f);

  gsl_odeiv2_driver_free(driver);
  return f[0];
}

int M_func(const gsl_vector* current_point, void* params, gsl_vector* M){
  double rmax = *(double *) params;
  double eps = gsl_vector_get(current_point,0);
  assert(eps<0);
  gsl_vector_set(M, 0, ode(rmax, eps));
  return GSL_SUCCESS;
}


int main() {

  gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid; //This is type
  gsl_multiroot_fsolver *s; //This is the workspace

  int dim = 1;
  double rmax = 8.0;

  gsl_multiroot_function F;
	F.f=M_func;
	F.n=dim;
  F.params = (void*) &rmax;

  int status;
  int iter=0;

  double eps_init[1] = {-10};
  gsl_vector * v = gsl_vector_alloc(dim);
  gsl_vector_set(v, 0, eps_init[0]);

  s = gsl_multiroot_fsolver_alloc(T, dim);
  gsl_multiroot_fsolver_set(s, &F, v);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    if (status)
      break;

    status = gsl_multiroot_test_residual(s->f, 1e-3);
  }
  while((status == GSL_CONTINUE) & (iter < 10000));

  double eps=gsl_vector_get(s->x,0);
  printf("epsilon = %g\n", eps);

  for(double r=0; r<=rmax; r+=0.1)
    printf("%g\t%g\t%g\n",r,ode(r,eps),r*exp(-r));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (v);

  return 0;
}
