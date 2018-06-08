#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void plainmc(int dim, double *a, double *b, double f(double* x), int N,
double* result, double* error);

double d2_gauss(double *x){
  return exp(-(x[0]*x[0])-(x[1]*x[1]));
}

int main(int argc, char const *argv[]) {
  int dim = 2;
  double a[2] = {-4.0, -4.0};
  double b[2] = {4.0, 4.0};
  double result, error;

  for(int N = 10; N < 1e+8; N *= 10) {
    plainmc(dim, a, b, d2_gauss, N, &result, &error);
    printf("%i, %g, %g, %g\n", N, 1/sqrt(N), error, fabs(result-M_PI));
  }


  return 0;
}
