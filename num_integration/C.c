#include<math.h>
#include<stdio.h>

int num_int;
double pi = 3.14159265359;

double my_integrator(double f(double), double a, double b, double acc, double eps,
  double *err);
double clenshaw_curtis(double f(double), double a, double b, double acc, double eps,
    double *err);

double function(double x){
  num_int ++;
  return 1/(x*x+1);
}

double function2(double x){
  num_int ++;
  return x*sin(2*exp(2*sin(2*exp(2*x))));
}

int main(int argc, char const *argv[]) {
  double a = -1;
  double b = 1;
  double acc = 1e-6;
  double eps = 1e-6;
  double error;

  num_int = 0;
  printf("1/(x*x+1) integrated from -1 to 1\n");
  double Q1 = my_integrator(function, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q1);
  printf("Exact value: %g\n", pi/2.0);
  printf("Need to satisfy %g < %g\n", fabs(Q1-pi/2.0), acc+fabs(Q1)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  num_int = 0;
  printf("With CC 1/(x*x+1) integrated from -1 to 1\n");
  double Q2 = clenshaw_curtis(function, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q2);
  printf("Exact value: %g\n", pi/2.0);
  printf("Need to satisfy %g < %g\n", fabs(Q2-pi/2.0), acc+fabs(Q2)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  num_int = 0;
  printf("x*sin(2*exp(2*sin(2*exp(2*x)))) integrated from -1 to 1\n");
  double Q3 = my_integrator(function2, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q3);
  printf("Exact value: %g\n", 0.336733);
  printf("Need to satisfy %g < %g\n", fabs(Q3-0.336733), acc+fabs(Q3)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  num_int = 0;
  printf("With CC x*sin(2*exp(2*sin(2*exp(2*x)))) integrated from -1 to 1\n");
  double Q4 = clenshaw_curtis(function2, a, b, acc, eps, &error);
  printf("Calculated: %g\n", Q4);
  printf("Exact value: %g\n", 0.336733);
  printf("Need to satisfy %g < %g\n", fabs(Q4-0.336733), acc+fabs(Q4)*eps);
  printf("Integrand evaluated %i times\n", num_int);
  printf("\n\n");

  printf("No visible improvements with CC\n");

  return 0;
}
