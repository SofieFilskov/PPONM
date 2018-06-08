#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

void qr_gs_decomp(gsl_matrix * A, gsl_matrix * R);
void backsub(gsl_matrix * R, gsl_matrix * Q, gsl_vector * x, gsl_vector * b);

void newton(void F(gsl_vector* x, gsl_vector* fx),
            gsl_vector* xstart,
            double dx,
            double eps) {
  int n = xstart->size;
  int steps = 0;
  gsl_matrix* J = gsl_matrix_alloc(n, n);
  gsl_vector* fx = gsl_vector_alloc(n);
  gsl_vector* xplus = gsl_vector_alloc(n);
  gsl_vector* fxplus = gsl_vector_alloc(n);
  gsl_matrix * R = gsl_matrix_alloc(n,n);
  gsl_vector* minusfx = gsl_vector_alloc(n);
  gsl_vector* deltax = gsl_vector_alloc(n);
  gsl_vector* lamdelx = gsl_vector_alloc(n);

  double lambda;
  double norm_fx;
  do {
    F(xstart, fx);
    for (int k= 0; k < n; k++) /*Jacobian*/ {
      gsl_vector_memcpy(xplus, xstart);
      gsl_vector_set(xplus, k, gsl_vector_get(xplus, k)+dx);
      F(xplus, fxplus);
      gsl_vector_sub(fxplus, fx);
      for( int i = 0; i < n; i++){
        gsl_matrix_set(J, i, k, gsl_vector_get(fxplus, i)/dx);
      }
    }
    gsl_vector_memcpy(minusfx, fx);
    gsl_vector_scale(minusfx,-1);
    qr_gs_decomp(J, R);
    backsub(R, J, deltax, minusfx);

    lambda = 1.0;
    norm_fx = gsl_blas_dnrm2(fx);
    gsl_vector_memcpy(lamdelx, deltax);
    do {
      lambda/=2;
      gsl_vector_scale(lamdelx, lambda);
      gsl_vector_add(xstart, lamdelx);
      F(xstart, fxplus);
    }
    while(gsl_blas_dnrm2(fxplus) > (1-lambda/2)*norm_fx && lambda > 1/128.0);
    steps++;
  }
  while(gsl_blas_dnrm2(fxplus) > eps && gsl_blas_dnrm2(deltax)>dx);
  printf("Newton done in %i steps\n", steps);

  gsl_matrix_free(J);
  gsl_vector_free(fx);
  gsl_vector_free(xplus);
  gsl_vector_free(fxplus);
  gsl_matrix_free(R);
  gsl_vector_free(minusfx);
  gsl_vector_free(deltax);
  gsl_vector_free(lamdelx);
}


void newton_jacobian(void FJ(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
            gsl_vector* xstart,
            double eps) {
  int n = xstart->size;
  int steps = 0;
  gsl_vector* fx = gsl_vector_alloc(n);
  gsl_vector* xplus = gsl_vector_alloc(n);
  gsl_vector* fxplus = gsl_vector_alloc(n);
  gsl_matrix * R = gsl_matrix_alloc(n,n);
  gsl_vector* minusfx = gsl_vector_alloc(n);
  gsl_vector* deltax = gsl_vector_alloc(n);
  gsl_vector* lamdelx = gsl_vector_alloc(n);
  gsl_matrix * J = gsl_matrix_alloc(n,n);

  double lambda;
  double norm_fx;

  do {
    FJ(xstart, fx, J);

    gsl_vector_memcpy(minusfx, fx);
    gsl_vector_scale(minusfx,-1);
    qr_gs_decomp(J, R);
    backsub(R, J, deltax, minusfx);

    lambda = 1.0;
    norm_fx = gsl_blas_dnrm2(fx);
    gsl_vector_memcpy(lamdelx, deltax);
    do {
      lambda/=2;
      gsl_vector_scale(lamdelx, lambda);
      gsl_vector_add(xstart, lamdelx);
      FJ(xstart, fxplus, J);
    }
    while(gsl_blas_dnrm2(fxplus) > (1-lambda/2)*norm_fx && lambda > 1/128.0);
    steps++;
  }
  while(gsl_blas_dnrm2(fxplus) > eps);
  printf("Newton with jacobian done in %i steps\n", steps);

  gsl_vector_free(fx);
  gsl_vector_free(xplus);
  gsl_vector_free(fxplus);
  gsl_matrix_free(R);
  gsl_vector_free(minusfx);
  gsl_vector_free(deltax);
  gsl_vector_free(lamdelx);
  gsl_matrix_free(J);
}
