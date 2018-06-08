#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>

void qr_gs_decomp(gsl_matrix * A, gsl_matrix * R);
void backsub(gsl_matrix * R, gsl_matrix * Q, gsl_vector * x, gsl_vector * b);


int qnewton(double F(gsl_vector* x, gsl_vector* df),
            gsl_vector* xstart,
            double eps) {
  int n = xstart->size;
  int steps = 0;
  int maxsteps = 1000;
  double alpha = 1e-4;
  double lambda;
  double f, f_xnew;
  double dotprod, sTy, uTy;

  gsl_vector* df = gsl_vector_alloc(n);
  gsl_vector* df_new = gsl_vector_alloc(n);
  gsl_matrix* B = gsl_matrix_alloc(n,n);
  gsl_vector* Bdf = gsl_vector_alloc(n);
  gsl_vector* xnew = gsl_vector_alloc(n);
  gsl_vector* dx = gsl_vector_alloc(n);
  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector* u = gsl_vector_alloc(n);

  gsl_matrix_set_identity(B);

  f = F(xstart, df);

  do {
    steps++;
    gsl_blas_dgemv(CblasNoTrans, -1.0, B, df, 0, dx); //dx = B*df
    if (gsl_blas_dnrm2(dx)<eps*gsl_blas_dnrm2(xstart)) {
      break;
    }
    if (gsl_blas_dnrm2(df)<eps) {
      break;
    }
    lambda = 1.0;
    do {
      gsl_vector_memcpy(xnew,xstart);
      gsl_vector_add(xnew, dx);
      f_xnew = F(xnew, df_new);
      gsl_blas_ddot(dx, df, &dotprod);

      if (f_xnew < f+alpha*dotprod) {
        //gsl_vector_memcpy(xstart, xnew);
        break;
      }

      lambda *= 0.5;
      gsl_vector_scale(dx, lambda);
    }
    while(lambda > 1.0/564.0);

    gsl_vector_memcpy(y, df_new);
    gsl_vector_sub(y, df); /*y = grad f(x+s) - grad f(x) */
    gsl_vector_memcpy(u, dx);
    gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u); /* u=s-*B*y */
    gsl_blas_ddot(dx,y,&sTy); /*sTy = s*y */

    if (fabs(sTy) > eps) {
      gsl_blas_ddot(u,y,&uTy);
			double gamma=uTy/2/sTy;
			gsl_blas_daxpy(-gamma,dx,u); /* u=u-gamma*s */
			gsl_blas_dger(1.0/sTy,u,dx,B); /*update B*/
			gsl_blas_dger(1.0/sTy,dx,u,B);
    }

    gsl_vector_memcpy(xstart, xnew);
    f = F(xstart, df);
  }
  while(gsl_blas_dnrm2(df) > eps && steps < maxsteps);


  gsl_vector_free(df);
  gsl_vector_free(Bdf);
  gsl_vector_free(y);
  gsl_vector_free(df_new);
  gsl_vector_free(xnew);
  gsl_matrix_free(B);
  gsl_vector_free(dx);

  return steps;
}
