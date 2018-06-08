#include<math.h>
#include<stdio.h>



double adapt24(double f(double), double a, double b, double acc, double eps,
double f2, double f3, double *interr, int num_rec) {
  if (num_rec > 1e+5){
    printf("Too many recursive calls, %i\n", num_rec);
    return 1;
  }

  double f1 = f(a+(b-a)/6.0);
  double f4 = f(a+5.0*(b-a)/6.0);

  double Q = (2.0*f1+f2+f3+2.0*f4)/6.0*(b-a);
  double q = (f1+f2+f3+f4)/4.0*(b-a);

  double tol = acc + eps * fabs(Q);
  double error = fabs(Q-q);
  if(error < tol){
    *interr = error;
    return Q;
  }
  else {
    double interr1, interr2;

    double Q1 = adapt24(f, a, (a+b)/2.0, acc/sqrt(2.0), eps, f1, f2, &interr1, num_rec+1.0);
    double Q2 = adapt24(f, (a+b)/2.0, b, acc/sqrt(2.0), eps, f3, f4, &interr2, num_rec+1.0);

    *interr = sqrt(interr1*interr1 + interr2*interr2);

    return Q1+Q2;
  }
}

double my_integrator(double f(double), double a, double b, double acc, double eps,
double *err) {
  int num_rec = 0;

  double f2 = f(a + 2.0*(b-a)/6.0);
  double f3 = f(a + 4.0*(b-a)/6.0);

  double Q = adapt24(f, a, b, acc, eps, f2, f3, err, num_rec);
  return Q;
}

double infty_integrator(double f(double), double a, double b, double acc, double eps,
double *err) {
  double Q;

  int a_inf = isinf(a); /* non-zero if argument is infinity */
  int b_inf = isinf(b);

  if(a_inf == 0 && b_inf == 0){
    Q = my_integrator(f, a, b, acc, eps, err);
  }
  else if(a_inf != 0 && b_inf != 0){
    double both_infty(double t) {
      return (f((1-t)/t)+f(-1*(1-t)/t))/t/t;
    }
    Q = my_integrator(both_infty, 0, 1, acc, eps, err);
  }
  else if(a_inf == 0 && b_inf != 0){
    double up_infty(double t) {
      return f(a+t/(1-t))/((1-t)*(1-t));
    }
    Q = my_integrator(up_infty, 0, 1, acc, eps, err);
  }
  else if(a_inf != 0 && b_inf == 0){
    double low_infty(double t) {
      return f(b-(1-t)/t)/(t*t);
    }
    Q = my_integrator(low_infty, 0, 1, acc, eps, err);
  }
  else {
    printf("Something is wrong with infinity...\n");
    Q = INFINITY;
  }
  return Q;
}

double clenshaw_curtis(double f(double), double a, double b, double acc, double eps,
double *err){
  double func(double t){
    return f(cos(t))*sin(t);
  }
  double pi = 3.14159265359;
  double Q = my_integrator(func, 0, pi, acc, eps, err);
}
