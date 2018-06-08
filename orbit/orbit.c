#include<stdio.h>
#include<gsl/gsl_odeiv2.h>
#include<math.h>
#include<gsl/gsl_errno.h>

int func_orbit(double phi, const double vec[], double dvecdphi[], void* params){
  double eps = *(double* )params;
  dvecdphi[0] = vec[1]; //u' = v
  dvecdphi[1] = 1-vec[0] + eps* pow(vec[0],2); //v' = u'' = 1 - u + eps*u^2
  return GSL_SUCCESS;
}

double ode(double phi, double eps, double start1, double start2){
  //making an ODE system
  //double eps = 0.01;
  gsl_odeiv2_system sys;
  sys.function = func_orbit;
  sys.jacobian = NULL;
  sys.dimension = 2;
  sys.params =  &eps;

  //Make it ode45
  gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;

  //The initial step size. We do so the step is negative if x is negative
  double hstart = copysign(0.1, phi);

  //Absolute and relative error
  double e_abs = 1e-6;
  double e_rel = 1e-6;

  //Makes the driver
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, \
    T, hstart, e_abs, e_rel);

  //Start conditions for applying the driver
  double phi0 = 0;
  double vec[2] = {start1, start2};
  gsl_odeiv2_driver_apply(driver, &phi0, phi, vec);

  gsl_odeiv2_driver_free(driver);
  return vec[0];
}

int main(){

  for (double phi = 0; phi < 3*3.1415; phi+=0.1){
    double eps = 0;
    double start1 = 1;
    double start2 = 0;
    printf("%g\t%g\n", phi, ode(phi, eps, start1, start2));
  }

  printf("\n\n");

  for (double phi = 0; phi < 3*3.1415; phi+=0.1){
    double eps = 0;
    double start1 = 1;
    double start2 = -0.5;
    printf("%g\t%g\n", phi, ode(phi, eps, start1, start2));
  }

  printf("\n\n");

  for (double phi = 0; phi < 30*3.1415; phi+=0.1){
    double eps = 0.01;
    double start1 = 1;
    double start2 = -0.5;
    printf("%g\t%g\n", phi, ode(phi, eps, start1, start2));
  }
  return 0;
}
