#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#define RND (double)rand()/RAND_MAX

int main(){
  int n = 3;
  int i;
  int j;

  //Making the matrix
  gsl_matrix* M = gsl_matrix_alloc (n,n);

  for (i = 0; i<n;i++)
    for (j = 0; j<n; j++)
      gsl_matrix_set (M, i, j, RND);

  for (i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      printf("M(%i, %i)\t=\t%g\n", i, j, gsl_matrix_get (M, i, j));

  //Making a copy
  gsl_matrix* A = gsl_matrix_alloc (n,n);
  gsl_matrix_memcpy(A, M);

  //Making the vector
  gsl_vector*  v = gsl_vector_alloc(n);

  for (i = 0; i < n; i++)
    gsl_vector_set(v, i, RND);

  for (i = 0; i <n; i++)
    printf("v(%i)\t=\t%g\n", i, gsl_vector_get(v,i));

  //Making an empty vector
  gsl_vector*  x = gsl_vector_alloc(n);

  //Solving the linear system
  gsl_linalg_HH_solve(M, v, x);

  for (i = 0; i <n; i++)
    printf("x(%i)\t=\t%g\n", i, gsl_vector_get(x,i));

  //Making an empty vector
  gsl_vector* b = gsl_vector_alloc(n);

  //Checking that the solution is right
  gsl_blas_dgemv (CblasNoTrans, 1, A, x, 0, b);

  printf("b should be equal to v\n");
  for (i = 0; i <n; i++)
    printf("b(%i) = %g\t v(%i) = %g\n", i, gsl_vector_get(b,i), i, \
      gsl_vector_get(v,i));

  gsl_matrix_free (M);
  gsl_vector_free(v);
  gsl_vector_free(x);
  gsl_matrix_free (A);
  gsl_vector_free(b);
  return 0;
}
