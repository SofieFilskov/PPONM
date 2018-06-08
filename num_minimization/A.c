#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);
int newton(
            double F(gsl_vector* x, gsl_vector * df, gsl_matrix * H),
            gsl_vector * xstart,
            double epsilon);

double rosen(gsl_vector * x, gsl_vector * df, gsl_matrix* H) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);

  gsl_vector_set(df, 0, -2.0*(1.0-xx)-400.0*xx*(yy-xx*xx));
  gsl_vector_set(df, 1, 200.0*(yy-xx*xx));

  gsl_matrix_set(H, 0, 0, 2.0-400.0*(yy-3.0*xx*xx));
  gsl_matrix_set(H, 0, 1, -400.0*xx);
  gsl_matrix_set(H, 1, 0, -400.0*xx);
  gsl_matrix_set(H, 1, 1, 200.0);

  return (1.0-xx)*(1.0-xx)+100.0*(yy-xx*xx)*(yy-xx*xx);
}

double himmel(gsl_vector * x, gsl_vector * df, gsl_matrix* H) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);

  gsl_vector_set(df, 0, 4*(xx*xx+yy-11.0)*xx+2*(xx+yy*yy-7.0));
  gsl_vector_set(df, 1, 2*(xx*xx+yy-11.0) + 4*(xx+yy*yy-7.0)*yy);

  gsl_matrix_set(H, 0, 0, 4.0*(3.0*xx*xx+yy-11.0)+2.0);
  gsl_matrix_set(H, 0, 1, 4.0*xx+4.0*yy);
  gsl_matrix_set(H, 1, 0, 4.0*xx+4.0*yy);
  gsl_matrix_set(H, 1, 1, 2.0+4.0*(xx+3*yy*yy-7));

  return (xx*xx+yy-11)*(xx*xx+yy-11)+(xx+yy*yy-7)*(xx+yy*yy-7);
}

int main(int argc, char const *argv[]) {
  printf("Rosenbrock function:\n");

  int m = 2; //because only 2 equations in system
  gsl_vector * x = gsl_vector_alloc(m);
  gsl_vector* df = gsl_vector_alloc(m);
  gsl_matrix* H = gsl_matrix_alloc(m, m);

  int steps;
  double f;

  gsl_vector_set(x,0,16);
  gsl_vector_set(x,1,12);
  f = rosen(x, df, H);
  printf("Starting x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);

  steps = newton(rosen, x, 1e-5);
  f = rosen(x, df, H);

  printf("Ending x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);
  printf("df(x) = [%g, %g]\n", gsl_vector_get(df,0), gsl_vector_get(df,1));
  printf("steps: %i\n", steps);

  printf("\nHimmelblau function:\n");

  gsl_vector_set(x,0,16);
  gsl_vector_set(x,1,12);
  f = himmel(x, df, H);
  printf("Starting x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);

  steps = newton(himmel, x, 1e-5);
  f = himmel(x, df, H);

  printf("x = [%g, %g]\n", gsl_vector_get(x,0), gsl_vector_get(x,1));
  printf("f(x) = %g\n", f);
  printf("df(x) = [%g, %g]\n", gsl_vector_get(df,0), gsl_vector_get(df,1));
  printf("steps: %i\n", steps);

  gsl_vector_free(x);
  gsl_vector_free(df);
  gsl_matrix_free(H);


  return 0;
}
