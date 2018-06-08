#include <stdio.h>
#include <gsl/gsl_integration.h>
//#include<gsl/gsl_errno.h>
#include <math.h>

double hamil (double x, void * params) {
  double alpha = *(double *) params;
  double f = (-alpha*alpha*x*x/2.0 + alpha/2.0 + x*x/2.0)*exp(-alpha*x*x);
  return f;
}

double norm (double x, void * params) {
  double alpha = *(double *) params;
  double f = exp(-alpha*x*x);
  return f;
}

int main() {
  int limit = 5000;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);

  double result_h, error_h, result_n, error_n;
  double err_abs = 1e-3;
  double err_rel = 1e-3;
  double alpha;

  gsl_function Hamil;
  Hamil.function = &hamil;
  Hamil.params = &alpha;
  gsl_function Norm;
  Norm.function = &norm;
  Norm.params = &alpha;

  for (double i = 0.4; i<1.40; i+=0.001){
    alpha = i;

    gsl_integration_qagi (&Hamil, err_abs, err_rel, limit, w, &result_h, &error_h);
    gsl_integration_qagi (&Norm, err_abs, err_rel, limit, w, &result_n, &error_n);

    double E = result_h/result_n;

    printf("%g\t%g\n", alpha, E);
  }

  gsl_integration_workspace_free (w);

  return 0;
}
