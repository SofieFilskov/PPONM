#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#define RND (double)rand()/RAND_MAX

void print_matrix(gsl_matrix * A);
void qr_gs_decomp(gsl_matrix * A, gsl_matrix * R);
void backsub (gsl_matrix* R, gsl_matrix* Q, gsl_vector* x, gsl_vector* b);
void inverse(gsl_matrix* R, gsl_matrix* Q, gsl_matrix* B);

int main() {
  /*        1        */
  int n = 3;

  gsl_matrix* A = gsl_matrix_alloc (n,n);
  gsl_matrix* R = gsl_matrix_alloc (n,n);

  for (int i = 0; i<n;i++){
    for (int j = 0; j<n; j++){
      gsl_matrix_set (A, i, j, RND);
    }
  }
  printf("Starting matrix A:\n");
  print_matrix(A);

  gsl_matrix* Q = gsl_matrix_alloc (n,n);
  gsl_matrix_memcpy(Q, A); //Making a copy, Q is equal A
  qr_gs_decomp(Q, R); //Now Q contains what is actually Q

  printf("Resulting matrix Q:\n");
  print_matrix(Q);
  printf("Resulting matrix R:\n");
  print_matrix(R);

  gsl_matrix* B = gsl_matrix_alloc (n,n);
  inverse(R, Q, B);
  printf("Inverse matrix of A:\n");
  print_matrix(B);

  gsl_matrix* I = gsl_matrix_alloc (n,n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, B, 0.0, I);
  printf("AB should be equal 1. AB is calculated to\n");
  print_matrix(I);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, B, A, 0.0, I);
  printf("BA should be equal 1. BA is calculated to\n");
  print_matrix(I);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(Q);
  gsl_matrix_free(B);
  gsl_matrix_free(I);

  return 0;
}
