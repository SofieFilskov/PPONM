//My number is 42 (yay), so the number exercise is 7..
#include<stdio.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <math.h>

typedef struct {double x;} parameter;

double cube(double c, void *params){
  // f(c) = c^3-x
  parameter *par = (parameter*) params;
	double x = par->x;

  return (c*c*c+x);
}

int main(int argc, char const *argv[]){
  double number = atof(argv[1]);
  int status;
  int iter = 0, max_iter = 100;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;

  parameter par;
  par.x = -1*number;
  double c_low = 0.0;
  double c_high = number;

  gsl_function F; //Giving the function
  F.function = &cube;
  F.params = (void*) &par;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);

  gsl_root_fsolver_set(s, &F, c_low, c_high);
  double epsabs = 1e-5;
  double epsrel = 1e-5;

  double r = gsl_root_fsolver_root (s);
  double r_exp = pow(number,1.0/3.0);
  printf ("%5d %.7f %.7f\n",
          iter, r, r_exp);

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      c_low = gsl_root_fsolver_x_lower (s);
      c_high = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (c_low, c_high,
                                       epsabs, epsrel);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d %.7f %.7f\n",
              iter, r, r_exp);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);
  return 0;
}
