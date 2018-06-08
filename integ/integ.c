#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <math.h>

double func (double x, void * params) {
  double f = log(x)/sqrt(x);
  return f;
}

int main() {
  int limit = 100;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);

  double result;
  double error;
  double err_abs = 1e-6;
  double err_rel = 1e-6;

  gsl_function F;
  F.function = &func;
  F.params = NULL;

  gsl_integration_qags (&F, 0, 1, err_abs, err_rel, limit, w, &result, &error);

  printf("int_0^1 ln(x) / sqrt(x) = %g\n", result);
  gsl_integration_workspace_free (w);

  return 0;
}
