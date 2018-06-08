#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

int func_calls = 0;

void qrdecomp(gsl_matrix * A, gsl_matrix * R);
void qrbacksub(gsl_matrix * Q, gsl_matrix * R, gsl_vector * b, gsl_vector * x);
void newton(
            void F(gsl_vector * x, gsl_vector * fx),
            gsl_vector * xstart,
            double dx,
            double epsilon);

void newton_jacobian(
            void FJ(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
            gsl_vector* xstart,
            double eps);

void F(gsl_vector * x, gsl_vector * fx) {
  double A = 10000;
  double f1 = A*gsl_vector_get(x,0)*gsl_vector_get(x,1)-1.0;
  double f2 = exp(-gsl_vector_get(x,0))+exp(-gsl_vector_get(x,1))-1.0-1.0/A;
  gsl_vector_set(fx, 0, f1);
  gsl_vector_set(fx, 1, f2);
  func_calls++;
}

void rosen(gsl_vector * x, gsl_vector * fx) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);
  gsl_vector_set(fx, 0, -2*(1-xx)-400*xx*(yy-xx*xx));
  gsl_vector_set(fx, 1, 200*(yy-xx*xx));
  func_calls++;
}

void himmel(gsl_vector * x, gsl_vector * fx) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);
  gsl_vector_set(fx, 0, 4*(xx*xx+yy-11.0)*xx+2*(xx+yy*yy-7.0));
  gsl_vector_set(fx, 1, 2*(xx*xx+yy-11.0) + 4*(xx+yy*yy-7.0)*yy);
  func_calls++;
}

void FJ(gsl_vector * x, gsl_vector * fx, gsl_matrix* J) {
  double A = 10000;
  double xx = gsl_vector_get(x,0);
  double yy = gsl_vector_get(x,1);
  double f1 = A*xx*yy-1.0;
  double f2 = exp(-xx)+exp(-yy)-1.0-1.0/A;
  gsl_vector_set(fx, 0, f1);
  gsl_vector_set(fx, 1, f2);

  gsl_matrix_set(J, 0, 0, A*yy);
  gsl_matrix_set(J, 0, 1, A*xx);
  gsl_matrix_set(J, 1, 0, -exp(-xx));
  gsl_matrix_set(J, 1, 1, -exp(-yy));
  func_calls++;
}

void rosen_jacobian(gsl_vector * x, gsl_vector * fx, gsl_matrix* J) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);
  gsl_vector_set(fx, 0, -2.0*(1.0-xx)-400.0*xx*(yy-xx*xx));
  gsl_vector_set(fx, 1, 200.0*(yy-xx*xx));

  gsl_matrix_set(J, 0, 0, 1200.0*xx*xx - 400.0 * yy + 2.0);
  gsl_matrix_set(J, 0, 1, -400.0*xx);
  gsl_matrix_set(J, 1, 0, -400.0*xx);
  gsl_matrix_set(J, 1, 1, 200.0);
  func_calls++;
}

void himmel_jacobian(gsl_vector * x, gsl_vector * fx, gsl_matrix* J) {
  double xx = gsl_vector_get(x, 0);
  double yy = gsl_vector_get(x, 1);
  gsl_vector_set(fx, 0, 4*(xx*xx+yy-11.0)*xx+2*(xx+yy*yy-7.0));
  gsl_vector_set(fx, 1, 2*(xx*xx+yy-11.0) + 4*(xx+yy*yy-7.0)*yy);

  gsl_matrix_set(J, 0, 0, 4*(3*xx*xx+yy-11)+2);
  gsl_matrix_set(J, 0, 1, 4*xx+4*yy);
  gsl_matrix_set(J, 1, 0, 4*xx+4*yy);
  gsl_matrix_set(J, 1, 1, 2+4*(xx+3*yy*yy-7));
  func_calls++;
}

int main(int argc, char const *argv[]) {
  int m = 2; //because only 2 equations in system
  gsl_vector * x = gsl_vector_calloc(m);
  gsl_vector * fx = gsl_vector_calloc(m);

  printf("The first function:\n");
  gsl_vector_set(x, 0, -2);
  gsl_vector_set(x, 1, 8);

  newton(F, x, 1e-6, 1e-3);
  F(x,fx);

  printf("x =\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x) = (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("# of functioncalls = %i \n",func_calls);
  func_calls=0;

  /*        Rosenbrock function       */
  printf("\nRosenbrock function:\n");
  gsl_vector_set(x,0,0.6);
  gsl_vector_set(x,1,1.4);

  newton(rosen, x, 1e-6,1e-3);
  rosen(x,fx);

  printf("x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("# of functioncalls = %i \n",func_calls);
  func_calls=0;

  /*        Himmelblau's function       */
  printf("\nHimmelblau's function:\n");
  gsl_vector_set(x,0,2.5);
  gsl_vector_set(x,1,2.5);

  newton(himmel, x, 1e-6,1e-3);
  himmel(x,fx);

  printf("x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("# of functioncalls = %i \n",func_calls);
  func_calls=0;

  /*        First function with jacobian       */
  printf("\nFirst function with jacobian:\n");
  gsl_matrix * J = gsl_matrix_calloc(m,m);
  gsl_vector_set(x, 0, -2);
  gsl_vector_set(x, 1, 8);

  newton_jacobian(FJ, x, 1e-3);
  FJ(x, fx, J);

  printf("x =\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x) = (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("# of functioncalls = %i \n",func_calls);
  func_calls=0;

  /*        Rosenbrock function with jacobian      */
  printf("\nRosenbrock function with jacobian:\n");
  gsl_vector_set(x,0,0.6);
  gsl_vector_set(x,1,1.4);

  newton_jacobian(rosen_jacobian, x, 1e-3);
  rosen_jacobian(x, fx, J);

  printf("x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("# of functioncalls = %i \n",func_calls);
  func_calls=0;

  /*        Himmelblau's function       */
  printf("\nHimmelblau's function with jacobian:\n");
  gsl_vector_set(x,0,2.5);
  gsl_vector_set(x,1,2.5);

  newton_jacobian(himmel_jacobian, x,1e-3);
  himmel_jacobian(x, fx, J);

  printf("x=\n");
  gsl_vector_fprintf(stdout,x,"%g");
  printf("f(x)=  (should be 0)\n");
  gsl_vector_fprintf(stdout,fx,"%g");
  printf("# of functioncalls = %i \n",func_calls);
  func_calls=0;

  gsl_vector_free(x);
  gsl_vector_free(fx);
  gsl_matrix_free(J);

  return 0;
}
