#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <math.h>

#define RND (double)rand()/RAND_MAX

void print_matrix(gsl_matrix * A);
int jacobi_cyclic(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row);

int main(int argc, char const *argv[]) {
  int n;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  else {
    n = 5;
  }

  gsl_matrix* A = gsl_matrix_alloc (n,n);

  for (int i = 0; i<n;i++){
    for (int j = 0; j<n; j++){
      double random = RND;
      gsl_matrix_set (A, i, j, random);
      gsl_matrix_set (A, j, i, random);
    }
  }
  printf("Starting matrix A:\n");
  print_matrix(A);

  gsl_matrix* V = gsl_matrix_alloc (n,n);
  gsl_vector* e = gsl_vector_alloc(n);

  int sweeps = 0;
  for (int i = 0; i < n; i++) {
    sweeps += jacobi_cyclic(A, e, V, i);
    printf("Row %i diagonal elements eliminated:\n", i);
    print_matrix(A);
    printf("The corresponding eigenvalue: %g\n\n", gsl_vector_get(e, i));
  }

  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_vector_free(e);

  return 0;
}
