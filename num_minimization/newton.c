#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

void qr_gs_decomp(gsl_matrix * A, gsl_matrix * R);
void backsub(gsl_matrix * R, gsl_matrix * Q, gsl_vector * x, gsl_vector * b);

int newton(double F(gsl_vector* x, gsl_vector* df, gsl_matrix* H),
            gsl_vector* xstart,
            double eps) {
  int n = xstart->size;
  int steps = 0;
  int maxsteps = 1000;
  double alpha = 1e-4;
  double y, yplus;
  double lambda;
  double dotprod;

  gsl_vector* df = gsl_vector_alloc(n);
  gsl_matrix* H = gsl_matrix_alloc(n,n);
  gsl_matrix * R = gsl_matrix_alloc(n,n);
  gsl_vector * xnew = gsl_vector_alloc(n);
  gsl_vector* dx = gsl_vector_alloc(n);

  do {
    y = F(xstart, df, H);
    gsl_vector_scale(df,-1.0);
    qr_gs_decomp(H,R);
    backsub(R, H, dx, df);

    lambda = 1.0;
    do {
      gsl_vector_memcpy(xnew,xstart);
      gsl_vector_scale(dx, lambda); //s
      gsl_vector_add(xnew, dx); //x+s
      yplus = F(xnew, df, R); //f(x+s)

      if (yplus < y+alpha*gsl_blas_ddot(dx, df, &dotprod)) {
        gsl_vector_memcpy(xstart, xnew);
        break;
      }

      lambda/=2.0;
    }
    while(lambda > 1.0/248.0);
    steps++;
  }
  while(gsl_blas_dnrm2(df) > eps && steps < maxsteps);

  gsl_matrix_free(H);
  gsl_vector_free(df);
  gsl_vector_free(xnew);
  gsl_matrix_free(R);
  gsl_vector_free(dx);

  return steps;
}
