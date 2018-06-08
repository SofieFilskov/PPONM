#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R) {
  int m = A->size2;
  for (int i = 0; i<m; i++){
    gsl_vector_view colAi = gsl_matrix_column(A, i);
    double normAi = gsl_blas_dnrm2(&colAi.vector);
    gsl_matrix_set(R, i, i, normAi);
    gsl_vector_scale(&colAi.vector, 1.0/normAi); //Now colAi = colQi = a_i/R_ii
    for (int j = i+1; j<m; j++){
      gsl_vector_view colAj = gsl_matrix_column(A, j);
      double dotprod = 0;
      gsl_blas_ddot(&colAi.vector, &colAj.vector, &dotprod);
      gsl_matrix_set(R, i, j, dotprod);
      gsl_blas_daxpy(-dotprod, &colAi.vector, &colAj.vector);
    }
  }
}

void backsub (gsl_matrix* R, gsl_matrix* Q, gsl_vector* x, gsl_vector* b){
  gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x);
  for (int i = x->size-1; i >= 0; i--){
    double s = gsl_vector_get(x,i);
    for (int k = i+1; k < x->size; k++){
      s -= gsl_matrix_get(R, i, k)*gsl_vector_get(x,k);
    }
    gsl_vector_set(x, i, s/gsl_matrix_get(R, i, i));
  }
}

void inverse(gsl_matrix* R, gsl_matrix* Q, gsl_matrix* B){
  int n = Q->size1;
  gsl_vector* e = gsl_vector_calloc(n);
  gsl_vector* x = gsl_vector_alloc(n);
  for (int i = 0; i < n; i++){
    gsl_vector_set(e, i, 1.0);
    backsub(R, Q, x, e);
    gsl_matrix_set_col(B, i, x);
    gsl_vector_set(e, i, 0.0);
  }
  gsl_vector_free(e);
  gsl_vector_free(x);
}
