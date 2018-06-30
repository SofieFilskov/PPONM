#include<math.h>
#include<stdio.h>
#include <omp.h>

int num_int;
double pi = 3.14159265359;

double my_integrator(double f(double), double a, double b, double acc, double eps,
  double *err);

double squareroot(double x){
  num_int ++;
  return sqrt(x);
}

double inv_sqrt(double x){
  num_int ++;
  return 1/sqrt(x);
}

double ln_inv_sqrt(double x){
  num_int ++;
  return log(x)/sqrt(x);
}

double function(double x){
  num_int ++;
  return 4*sqrt(1-(1-x)*(1-x));
}

int main(int argc, char const *argv[]) {
  double a = 0;
  double b = 1;
  double acc = 1e-6;
  double eps = 1e-6;
  double error;
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      printf("first section started\n");
      num_int = 0;
      printf("Sqrt(x) integrated from 0 to 1\n");
      double Q1 = my_integrator(squareroot, a, b, acc, eps, &error);
      printf("Calculated: %g\n", Q1);
      printf("Exact value: %g\n", 2.0/3.0);
      printf("Need to satisfy %g < %g\n", fabs(Q1-2.0/3.0), acc+fabs(Q1)*eps);
      printf("\n\n");
      printf("first section finished\n");
    }
    #pragma omp section
    {
      printf("second section started\n");
      num_int = 0;
      printf("1/sqrt(x) integrated from 0 to 1\n");
      double Q2 = my_integrator(inv_sqrt, a, b, acc, eps, &error);
      printf("Calculated: %g\n", Q2);
      printf("Exact value: %g\n", 2.0);
      printf("Need to satisfy %g < %g\n", fabs(Q2-2.0), acc+fabs(Q2)*eps);
      printf("\n\n");
      printf("second section finished\n");
    }
    #pragma omp section
    {
      printf("third section started\n");
      num_int = 0;
      printf("ln(x)/sqrt(x) integrated from 0 to 1\n");
      double Q3 = my_integrator(ln_inv_sqrt, a, b, acc, eps, &error);
      printf("Calculated: %g\n", Q3);
      printf("Exact value: %g\n", -4.0);
      printf("Need to satisfy %g < %g\n", fabs(Q3+4.0), acc+fabs(Q3)*eps);
      printf("\n\n");
      printf("third section finished\n");
    }
  }

  num_int = 0;
  acc = 1e-10;
  eps = 1e-10;
  printf("4*sqrt(1-(1-x^2)) integrated from 0 to 1\n");
  double Q4 = my_integrator(function, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q4);
  printf("Exact value: %g\n", pi);
  printf("Need to satisfy %g < %g\n", fabs(Q4-pi), acc+fabs(Q4)*eps);
  printf("Integrand evaluated %i times\n", num_int);


  return 0;
}
