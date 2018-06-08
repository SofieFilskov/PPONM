#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void plainmc(int dim, double *a, double *b, double f(double* x), int N,
double* result, double* error);

double inv_sqrt(double *x){
  return 1/sqrt(x[0]);
}

double d2_gauss(double *x){
  return exp(-(x[0]*x[0])-(x[1]*x[1]));
}

double function(double *x){
  return 1.0/((1.0-cos(x[0])*cos(x[1])*cos(x[2]))*M_PI*M_PI*M_PI);
}

int main(int argc, char const *argv[]) {
  int dim = 1;
  double a[1] = {0.0};
  double b[1] = {1.0};
  double result, error;
  int N = 10000;

  plainmc(dim, a, b, inv_sqrt, N, &result, &error);

  printf("1/qrt(x) integrated from 0 to 1\n");
  printf("Calculated to %g with error %g\n", result, error);
  printf("Exact value: %g\n", 2.0);
  printf("Real error: %g\n", fabs(result-2.0));
  printf("\n\n");

  int dim2 = 2;
  double a2[2] = {-2.0, -2.0};
  double b2[2] = {2.0, 2.0};
  double result2, error2;
  int N2 = 10000;

  plainmc(dim2, a2, b2, d2_gauss, N2, &result2, &error2);

  printf("2D gauss integrated from 0 to 1\n");
  printf("Calculated to %g with error %g\n", result2, error2);
  printf("Exact value: %g\n", M_PI);
  printf("Real error: %g\n", fabs(result2-M_PI));
  printf("\n\n");

  int dim3 = 3;
  double a3[3] = {0, 0, 0};
  double b3[3] = {M_PI, M_PI, M_PI};
  double result3, error3;
  int N3 = 100000;

  plainmc(dim3, a3, b3, function, N3, &result3, &error3);

  printf("Long function integrated from 0 to pi\n");
  printf("Calculated to %g with error %g\n", result3, error3);
  printf("Exact value: %g\n", 1.3932039296856768591842462603255);
  printf("Real error: %g\n", fabs(result3-1.3932039296856768591842462603255));
  

  return 0;
}
