#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <math.h>

#define RND (double)rand()/RAND_MAX

void print_matrix(gsl_matrix * A);
int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, gsl_matrix* U, int* rot);

int main(int argc, char const *argv[]) {
  int n = 2;

  gsl_matrix* A = gsl_matrix_alloc (n,n);
  gsl_matrix* A2 = gsl_matrix_alloc (n,n);
  gsl_matrix* ATA = gsl_matrix_alloc (n,n);
  gsl_matrix* AAT = gsl_matrix_alloc (n,n);
  gsl_matrix* D = gsl_matrix_alloc (n,n);
  gsl_matrix* V = gsl_matrix_alloc (n,n);
  gsl_matrix* U = gsl_matrix_alloc (n,n);
  gsl_matrix* VTATA = gsl_matrix_alloc (n,n);
  gsl_matrix* VTATAV = gsl_matrix_alloc (n,n);
  gsl_matrix* UTAAT = gsl_matrix_alloc (n,n);
  gsl_matrix* UTAATU = gsl_matrix_alloc (n,n);
  gsl_matrix* DVT = gsl_matrix_alloc (n,n);
  gsl_matrix* VVT = gsl_matrix_alloc (n,n);
  gsl_matrix* UUT = gsl_matrix_alloc (n,n);
  gsl_vector* e = gsl_vector_alloc(n);

  for (int i = 0; i<n; i++){
    for (int j = 0; j<n; j++){
      double random = RND;
      gsl_matrix_set (A, i, j, random);
    }
  }
  gsl_matrix_memcpy(D, A);
  printf("Starting matrix A:\n");
  print_matrix(A);

  int rot = 0;
  int sweeps = jacobi(D, e, V, U, &rot);

  printf("D from jacobi:\n");
  for (int i = 0; i<n; i++){
    gsl_matrix_set(D, i, i, gsl_vector_get(e, i));
  }
  /* Making sure that D is positive */
  for (int i = 0; i<n; i++){
    double Dii = gsl_matrix_get(D, i, i);
    if (Dii < 0){
      gsl_matrix_set(D, i, i, Dii*-1);
      for (int j = 0; j<n; j++){
        double Uij = gsl_matrix_get(U, j, i);
        gsl_matrix_set(U, j, i, Uij*-1);
      }
    }
  }
  print_matrix(D);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, D, V, 0, DVT);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, DVT, 0, A2);
  printf("Calculated A = UDVT:\n");
  print_matrix(A2);

  printf("Calculated VTATAV = eigenvalues for ATA:\n");
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, ATA);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, V, ATA, 0, VTATA);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, VTATA, V, 0, VTATAV);
  print_matrix(VTATAV);

  printf("Calculated UTAATU = eigenvalues for AAT:\n");
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, A, A, 0, AAT);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, U, AAT, 0, UTAAT);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, UTAAT, U, 0, UTAATU);
  print_matrix(UTAATU);

  printf("Eigenvalues from D = sqrt(lambda):\n");
  for (int i = 0; i<n; i++){
    printf("%g\n", pow(gsl_vector_get(e, i),2));
  }
  printf("\n");

  printf("Checking that VVT is equal I:\n");
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, V, V, 0, VVT);
  print_matrix(VVT);
  printf("Checking that UUT is equal I:\n");
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, U, U, 0, UUT);
  print_matrix(UUT);

  printf("Sweeps: %i, rotations: %i\n\n", sweeps, rot);


  gsl_matrix_free(A);
  gsl_matrix_free(A2);
  gsl_matrix_free(D);
  gsl_matrix_free(V);
  gsl_matrix_free(U);
  gsl_matrix_free(ATA);
  gsl_matrix_free(AAT);
  gsl_matrix_free(VTATA);
  gsl_matrix_free(VTATAV);
  gsl_matrix_free(UTAAT);
  gsl_matrix_free(UTAATU);
  gsl_matrix_free(DVT);
  gsl_matrix_free(VVT);
  gsl_matrix_free(UUT);
  gsl_vector_free(e);


  return 0;
}
