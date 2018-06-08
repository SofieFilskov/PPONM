#include<stdio.h>
#include<gsl/gsl_odeiv2.h>
#include<math.h>
#include<gsl/gsl_errno.h>

int func(double x, const double y[], double dydx[], void* params){
  dydx[0] = y[0] * (1-y[0]);
  return GSL_SUCCESS;
}

double ode(double x){
  //making an ODE system
  gsl_odeiv2_system sys;
  sys.function = func;
  sys.jacobian = NULL;
  sys.dimension = 1;
  sys.params =  NULL;

  //Make it ode45
  gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;

  //The initial step size. We do so the step is negative if x is negative
  double hstart = copysign(0.1, x);

  //Absolute and relative error
  double e_abs = 1e-6;
  double e_rel = 1e-6;

  //Makes the driver
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, \
    T, hstart, e_abs, e_rel);

  //Start conditions for applying the driver
  double x0 = 0;
  double y[1] = {0.5};
  gsl_odeiv2_driver_apply(driver, &x0, x, y);

  gsl_odeiv2_driver_free(driver);
  return y[0];
}


int main() {
  for (double x = 0; x < 3.1; x+=0.1)
    printf("%g\t%g\n", x, ode(x));

  printf("\n\n");

  for (double x = 0; x < 3.1; x+=0.1)
    printf("%g\t%g\n", x, 1/(1+exp(-x)));

  return 0;
}
