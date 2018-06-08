#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#define RND (double)rand()/RAND_MAX

void print_matrix(gsl_matrix * A);
void qr_gs_decomp(gsl_matrix * A, gsl_matrix * R);
void backsub (gsl_matrix* R, gsl_matrix* Q, gsl_vector* x, gsl_vector* b);

int main() {
  /*        1        */
  int n = 3;
  int m = 3;

  gsl_matrix* A = gsl_matrix_alloc (n,m);
  gsl_matrix* R = gsl_matrix_alloc (m,m);

  for (int i = 0; i<n;i++){
    for (int j = 0; j<m; j++){
      gsl_matrix_set (A, i, j, RND);
    }
  }
  printf("Starting matrix A:\n");
  print_matrix(A);

  gsl_matrix* Q = gsl_matrix_alloc (n,m);
  gsl_matrix_memcpy(Q, A); //Making a copy, Q is equal A
  qr_gs_decomp(Q, R); //Now Q contains what is actually Q

  printf("Resulting matrix Q:\n");
  print_matrix(Q);
  printf("Resulting matrix R:\n");
  print_matrix(R);

  printf("QTQ should be 1\nCalculated QTQ:\n");
  gsl_matrix* QTQ = gsl_matrix_alloc (m,m);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Q, Q, 0.0, QTQ);
  print_matrix(QTQ);

  printf("QR should be equal A\nCalculated QR-A:\n");
  gsl_matrix* QR = gsl_matrix_alloc (n,m);
  gsl_matrix_memcpy(QR, A); //Making a copy, Q is equal A
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q, R, -1.0, QR);
  print_matrix(QR);

  /*        2        */
  gsl_vector* b = gsl_vector_alloc(m);
  for (int i = 0; i < m; i++) {
    gsl_vector_set(b, i, RND);
  }
  printf("Vector b:\n");
  gsl_vector_fprintf(stdout, b, "%g");

  gsl_vector* x = gsl_vector_alloc(m);
  backsub(R, Q, x, b); //Only works for square matrix

  printf("\nSolution x:\n");
  gsl_vector_fprintf(stdout, x, "%g");

  printf("\nShould be Ax=b, calculated Ax to:\n");
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, b);
  gsl_vector_fprintf(stdout, b, "%g");

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(Q);
  gsl_matrix_free(QTQ);
  gsl_matrix_free(QR);
  gsl_vector_free(b);
  gsl_vector_free(x);

  return 0;
}
