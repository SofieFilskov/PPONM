#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <math.h>

void rkstep12(
	double t,                                  /* the current value of the variable */
	double h,                                  /* the step to be taken */
	gsl_vector* yt,                                /* the current value y(t) of the sought function */
	void f(double t, gsl_vector* y, gsl_vector* dydt), /* the right-hand-side, dydt = f(t,y) */
	gsl_vector* yth,                               /* output: y(t+h) */
	gsl_vector* err                                /* output: error estimate dy */
  ) {

  int i;
  int n = yt->size;
  gsl_vector* k0 = gsl_vector_alloc(n);
  gsl_vector* yt_plus = gsl_vector_alloc(n);
  gsl_vector* k12 = gsl_vector_alloc(n);

  /* Making y0 + 1/2 h k0 to be used later*/
  f(t, yt, k0); /* k0 = f(t0, y0) */
  for (i = 0; i < n; i++) {
    gsl_vector_set(yt_plus, i, gsl_vector_get(yt, i) + gsl_vector_get(k0,i)*h*0.5);
  }

  /* 2nd order midpoint */
  f(t+h*0.5, yt_plus, k12); /* k1/2 = f(t0+1/2h, y0+1/2*h*k0) */
  for(i = 0; i < n; i++) {
    gsl_vector_set(yth, i, gsl_vector_get(yt,i) + gsl_vector_get(k12,i) * h); /* y(t+h) = y(t) + k1/2 * h */
  }

  /* error estimation */
  for(i = 0; i < n; i++) {
    gsl_vector_set(err,i,(gsl_vector_get(k0,i) - gsl_vector_get(k12,i))*h/2);
  }

gsl_vector_free(k0);
gsl_vector_free(yt_plus);
gsl_vector_free(k12);

}

int driver(
	double* t,               /* the current value of the variable */
	double b,                /* the end-point of the integration */
	double* h,               /* the current step-size */
	gsl_vector* yt,          /* the current y(t) */
	double acc,              /* absolute accuracy goal */
	double eps,              /* relative accuracy goal */
  void rkstep12(
  	double t, double h, gsl_vector* yt,
  	void f(double t, gsl_vector* y, gsl_vector* dydt),
  	gsl_vector* yth, gsl_vector* err),
	void f(double t, gsl_vector* y, gsl_vector* dydt) /* right-hand-side */
) {
  int i, k = 0;
  int n = yt->size;
  double normy, loc_err, tau;
  double x = *t;
  double a = x;
  double hint = *h; /*h_internal*/
  gsl_vector* yth = gsl_vector_alloc(n);
  gsl_vector* err = gsl_vector_alloc(n);

  while(x < b){
    if (x+hint > b) {
      hint = b-x;
    }
    rkstep12(x, hint, yt, f, yth, err);

    normy = gsl_blas_dnrm2(yth);
    loc_err = gsl_blas_dnrm2(err);
    tau = (normy*eps + acc)*sqrt(hint/(b-a));
    if(loc_err < tau) {
      k++;
      x = x + hint;
      gsl_vector_memcpy(yt, yth);

      printf("%g %g %g\n",x, gsl_vector_get(yth, 0), gsl_vector_get(yth, 1));
      
    }
    if (loc_err > 0) {
      hint *= pow(tau/loc_err, 0.25)*0.95;
    }
    else {
      hint *= 2.0;
    }
  }
  *h = hint;

  gsl_vector_free(yth);
  gsl_vector_free(err);

  return k+1; /*Number of steps we have taken*/
}
