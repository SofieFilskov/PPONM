#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <math.h>

#define RND (double)rand()/RAND_MAX

void print_matrix(gsl_matrix * A);
int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V);

int main(int argc, char const *argv[]) {
  int n;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  else {
    n = 5;
  }

  gsl_matrix* A = gsl_matrix_alloc (n,n);
  gsl_matrix* A_copy = gsl_matrix_alloc (n,n);

  for (int i = 0; i<n;i++){
    for (int j = 0; j<n; j++){
      double random = RND;
      gsl_matrix_set (A, i, j, random);
      gsl_matrix_set (A, j, i, random);
    }
  }
  printf("Starting matrix A:\n");
  print_matrix(A);
  gsl_matrix_memcpy(A_copy, A); //Making a copy, Q is equal A

  gsl_matrix* V = gsl_matrix_alloc (n,n);
  gsl_vector* e = gsl_vector_alloc(n);

  int sweeps = jacobi(A_copy, e, V);

  printf("Orthogonal matrix V with eigenvectors:\n");
  print_matrix(V);

  printf("Vector e with eigenvalues:\n");
  gsl_vector_fprintf(stdout, e, "%g");

  printf("\nSweeps is equal %i\n\n", sweeps);

  gsl_matrix* D1 = gsl_matrix_alloc (n,n);
  gsl_matrix* D2 = gsl_matrix_alloc (n,n);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, V, A, 0.0, D1);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, D1, V, 0.0, D2);

  for (int i = 0; i<n;i++){
    for (int j = 0; j<n; j++){
      double Dij = gsl_matrix_get(D2, i, j);
      if (fabs(Dij) < 0.00001)
        gsl_matrix_set (D2, i, j, 0);
    }
  }

  printf("The matrix VTAV = D is calculated to:\n");
  print_matrix(D2);

  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_vector_free(e);
  gsl_matrix_free(D1);
  gsl_matrix_free(D2);
  gsl_matrix_free(A_copy);

  return 0;
}
