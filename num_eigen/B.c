#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <math.h>

#define RND (double)rand()/RAND_MAX

void print_matrix(gsl_matrix * A);
int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int* rot);
int jacobi_B(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row, int* rot);
int jacobi_B_high(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row, int* rot);

int main(int argc, char const *argv[]) {
  int n;
  if (argc > 1) {
    n = atoi(argv[1]);
  }
  else {
    n = 5;
  }

  printf("Fully diagonalizing matrix, lowest eigenvalues first\n");

  gsl_matrix* A = gsl_matrix_alloc (n,n);
  gsl_matrix* A_copy = gsl_matrix_alloc (n,n);
  gsl_matrix* A_copy2 = gsl_matrix_alloc (n,n);
  gsl_matrix* A_copy3 = gsl_matrix_alloc (n,n);

  for (int i = 0; i<n;i++){
    for (int j = 0; j<n; j++){
      double random = RND;
      gsl_matrix_set (A, i, j, random);
      gsl_matrix_set (A, j, i, random);
    }
  }
  printf("Starting matrix A:\n");
  print_matrix(A);
  gsl_matrix_memcpy(A_copy, A);
  gsl_matrix_memcpy(A_copy2, A);
  gsl_matrix_memcpy(A_copy3, A);

  gsl_matrix* V = gsl_matrix_alloc (n,n);
  gsl_matrix_set_identity(V);

  gsl_vector* e = gsl_vector_alloc(n);
  for(int i = 0; i < n; i++){
		gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
	}

  int rot = 0;
  int sweeps = 0;
  for (int k = 0; k < n; k++) {
    sweeps += jacobi_B(A, e, V, k, &rot);

    for (int j = k; j<n; j++){
      double Akj = gsl_matrix_get(A, k, j);
      if (Akj < 1e-5){
        gsl_matrix_set(A, k, j, 0);
      }
    }

    printf("Row %i diagonal elements eliminated:\n", k);
    print_matrix(A);
    printf("The corresponding eigenvalue: %g\n\n", gsl_vector_get(e, k));
  }

  printf("Sweeps: %i, Rotations: %i\n\n", sweeps, rot);

  /* Largest eigenvalues */
  gsl_matrix_set_identity(V);
  for(int i = 0; i < n; i++){
		gsl_vector_set(e, i, gsl_matrix_get(A_copy, i, i));
	}
  printf("Finding highest eigenvalues first\n");
  sweeps = 0;
  rot = 0;
  for (int k = 0; k < n; k++) {
    sweeps += jacobi_B_high(A_copy, e, V, k, &rot);

    for (int j = k; j<n; j++){
      double Akj = gsl_matrix_get(A_copy, k, j);
      if (Akj < 1e-5){
        gsl_matrix_set(A_copy, k, j, 0);
      }
    }

    printf("%g\n", gsl_vector_get(e, k));
  }
  printf("Sweeps: %i, Rotations: %i\n\n", sweeps, rot);

  /* Only the first from small to large */
  gsl_matrix_set_identity(V);
  for(int i = 0; i < n; i++){
		gsl_vector_set(e, i, gsl_matrix_get(A_copy3, i, i));
	}
  printf("Finding only the first lowest eigenvalue\n");
  sweeps = 0;
  rot = 0;
  for (int k = 0; k < 1; k++) {
    sweeps += jacobi_B(A_copy3, e, V, k, &rot);

    for (int j = k; j<n; j++){
      double Akj = gsl_matrix_get(A_copy3, k, j);
      if (Akj < 1e-5){
        gsl_matrix_set(A_copy3, k, j, 0);
      }
    }

    printf("%g\n", gsl_vector_get(e, k));
  }
  printf("Sweeps: %i, Rotations: %i\n\n", sweeps, rot);

  /*   Using cyclic method   */

  gsl_matrix_set_identity(V);
  for(int i = 0; i < n; i++){
		gsl_vector_set(e, i, gsl_matrix_get(A_copy2, i, i));
	}
  printf("Cyclic method, eigenvalues are:\n");
  sweeps = 0;
  rot = 0;

  sweeps = jacobi(A_copy2, e, V, &rot);

  for (int j = 0; j<n; j++){
    for (int i = 0; i < n; i++){
      double Aji = gsl_matrix_get(A_copy2, j, i);
      if (Aji < 1e-5){
        gsl_matrix_set(A_copy2, j, i, 0);
      }
    }
    printf("%g\n", gsl_vector_get(e, j));
  }

  printf("Sweeps: %i, Rotations: %i\n\n", sweeps, rot);



  gsl_matrix_free(A);
  gsl_matrix_free(A_copy);
  gsl_matrix_free(A_copy2);
  gsl_matrix_free(A_copy3);
  gsl_matrix_free(V);
  gsl_vector_free(e);

  return 0;
}
