#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);

int qnewton(
            double F(gsl_vector* x, gsl_vector * df),
            gsl_vector * xstart,
            double epsilon);

double rosen(gsl_vector * x, gsl_vector * df) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);

  gsl_vector_set(df, 0, -2.0*(1.0-xx)-400.0*xx*(yy-xx*xx));
  gsl_vector_set(df, 1, 200.0*(yy-xx*xx));

  return (1.0-xx)*(1.0-xx)+100.0*(yy-xx*xx)*(yy-xx*xx);
}

double himmel(gsl_vector * x, gsl_vector * df) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);

  gsl_vector_set(df, 0, 4*(xx*xx+yy-11.0)*xx+2*(xx+yy*yy-7.0));
  gsl_vector_set(df, 1, 2*(xx*xx+yy-11.0) + 4*(xx+yy*yy-7.0)*yy);

  return (xx*xx+yy-11)*(xx*xx+yy-11)+(xx+yy*yy-7)*(xx+yy*yy-7);
}

double radio(gsl_vector * x, gsl_vector * df) {
  double A = gsl_vector_get(x, 0);
  double T = gsl_vector_get(x, 1);
  double B = gsl_vector_get(x, 2);

  double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
  double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
  double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
  int N = sizeof(t)/sizeof(t[0]);

  double sum_dA = 0;
  double sum_dT = 0;
  double sum_dB = 0;
  double sum_F = 0;
  double dFdA;
  double dFdT;
  double dFdB;
  double F;
  for (int i=0; i<N; i++) {
    dFdA = 2*((A*exp(-t[i]/T)+B)-y[i])*((A*exp(-t[i]/T)+B)-y[i])/(e[i]*e[i])\
          *exp(-t[i]/T);
    sum_dA = sum_dA + dFdA;
    dFdT = 2*((A*exp(-t[i]/T)+B)-y[i])*((A*exp(-t[i]/T)+B)-y[i])/(e[i]*e[i])\
          *A*exp(-t[i]/T)/(T*T);
    sum_dT = sum_dT + dFdT;
    dFdB = 2*((A*exp(-t[i]/T)+B)-y[i])*((A*exp(-t[i]/T)+B)-y[i])/(e[i]*e[i]);
    sum_dB = sum_dB + dFdB;

    F = ((A*exp(-t[i]/T)+B)-y[i])*((A*exp(-t[i]/T)+B)-y[i])/(e[i]*e[i]);
    sum_F = sum_F + F;
  }

  gsl_vector_set(df, 0, sum_dA);
  gsl_vector_set(df, 1, sum_dT);
  gsl_vector_set(df, 2, sum_dB);

  return sum_F;
}

int main(int argc, char const *argv[]) {
  printf("Rosenbrock function:\n");

  int m = 2; //because only 2 equations in system
  gsl_vector * x = gsl_vector_alloc(m);
  gsl_vector* df = gsl_vector_alloc(m);

  int steps;
  double f;

  gsl_vector_set(x,0,16);
  gsl_vector_set(x,1,12);
  f = rosen(x, df);
  printf("Starting x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);

  steps = qnewton(rosen, x, 1e-5);
  f = rosen(x, df);

  printf("Ending x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);
  printf("df(x) = [%g, %g]\n", gsl_vector_get(df,0), gsl_vector_get(df,1));
  printf("steps: %i\n", steps);

  printf("\nHimmelblau function:\n");

  gsl_vector_set(x,0,16);
  gsl_vector_set(x,1,12);
  f = himmel(x, df);
  printf("Starting x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);

  steps = qnewton(himmel, x, 1e-5);
  f = himmel(x, df);

  printf("Ending x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);
  printf("df(x) = [%g, %g]\n", gsl_vector_get(df,0), gsl_vector_get(df,1));
  printf("steps: %i\n", steps);

  gsl_vector_free(x);
  gsl_vector_free(df);

  /*    Radioactivity   */
  printf("\n\n");

  double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
  double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
  double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
  int N = sizeof(t)/sizeof(t[0]);

  for (int i=0; i<N; i++) {
    printf("%g\t%g\t%g\n", t[i], y[i], e[i]);
  }
  printf("\n\n");

  m = 3; //because only 2 equations in system
  gsl_vector * x2 = gsl_vector_alloc(m);
  gsl_vector* df2 = gsl_vector_alloc(m);

  gsl_vector_set(x2,0,3.5);
  gsl_vector_set(x2,1,2.8);
  gsl_vector_set(x2,2,1.5);
  f = radio(x2, df2);
  printf("Starting x = [%g, %g, %g]\n", gsl_vector_get(x2,0), gsl_vector_get(x2,1),\
          gsl_vector_get(x2,2));
  printf("f(x) = %g\n", f);

  steps = qnewton(radio, x2, 1e-5);
  f = radio(x2, df2);

  double A = gsl_vector_get(x2,0);
  double T = gsl_vector_get(x2,1);
  double B = gsl_vector_get(x2,2);

  printf("Ending x = [%g, %g, %g]\n", A, T, B);
  printf("f(x) = %g\n", f);
  printf("df(x) = [%g, %g, %g]\n", gsl_vector_get(df2,0), gsl_vector_get(df2,1), \
          gsl_vector_get(df2,2));
  printf("steps: %i\n", steps);

  printf("\n\n");

  double F;
  for (int i = 0; i < N; i++) {
    F = A*exp(-t[i]/T)+B;
    printf("%g\t%g\n", t[i], F);
  }


  gsl_vector_free(x2);
  gsl_vector_free(df2);


  return 0;
}
