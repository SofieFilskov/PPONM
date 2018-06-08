#include<math.h>
#include<stdio.h>
#include <gsl/gsl_integration.h>

int num_int;
double pi = 3.14159265359;

double infty_integrator(double f(double), double a, double b, double acc, double eps,
  double *err);

double function_inf(double x){
  num_int ++;
  return exp(-(x*x));
}

double function_both(double f(double), double t){
  return (function_inf((1-t)/t)+function_inf(-(1-t)/t))/t/t;
}

double function_up(double f(double), double t){
  return function_inf(1+(1-t)/t)/t/t;
}

double function_low(double f(double), double t){
  return function_inf(1-(1-t)/t)/t/t;
}


int main(int argc, char const *argv[]) {
  double acc = 1e-4;
  double eps = 1e-4;
  double error;

  double a = -INFINITY;
  double b = INFINITY;
  num_int = 0;
  printf("exp(x^2) integrated from -infty to infty\n");
  double Q1 = infty_integrator(function_inf, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q1);
  printf("Exact value: %g\n", sqrt(pi));
  printf("Need to satisfy %g < %g\n", fabs(Q1-sqrt(pi)), acc+fabs(Q1)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  a = 1;
  b = INFINITY;
  num_int = 0;
  printf("exp(x^2) integrated from 1 to infty\n");
  double Q2 = infty_integrator(function_inf, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q2);
  printf("Exact value: %g\n", 0.139403);
  printf("Need to satisfy %g < %g\n", fabs(Q2-0.139403), acc+fabs(Q2)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  a = -INFINITY;
  b = 1;
  num_int = 0;
  printf("exp(x^2) integrated from -infty to 1\n");
  double Q3 = infty_integrator(function_inf, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q3);
  printf("Exact value: %g\n", 1.63305);
  printf("Need to satisfy %g < %g\n", fabs(Q3-1.63305), acc+fabs(Q3)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  /* Checking with GSL */
  int limit = 50000;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);

  gsl_function F;
  F.function = &function_both;
  F.params = NULL;
  double Q4;
  gsl_integration_qags (&F, 0, 1, acc, eps, limit, w, &Q4, &error);
  printf("With GSL integration from -infty to infty exp(-x^2) = %g\n", Q4);

  F.function = &function_up;
  F.params = NULL;
  double Q5;
  gsl_integration_qags (&F, 0, 1, acc, eps, limit, w, &Q5, &error);
  printf("With GSL integration from 1 to infty exp(-x^2) = %g\n", Q5);

  F.function = &function_low;
  F.params = NULL;
  double Q6;
  gsl_integration_qags (&F, 0, 1, acc, eps, limit, w, &Q6, &error);
  printf("With GSL integration from -infty to 1 exp(-x^2) = %g\n", Q6);

  gsl_integration_workspace_free (w);

  return 0;
}
