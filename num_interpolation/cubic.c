#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef struct {int n; double * x, *y, *b, *c;} cspline;

cspline * cspline_alloc(int n, double * x, double * y);
double cspline_eval(cspline * s, double z);
double cspline_integ(cspline * s, double z);
double cspline_deriv(cspline * s, double z);
void cspline_free(cspline * s);

int main() {
  int n = 11;
  int fine_i = 10;
  double fine_d = 10.0;

  double x[n];
  double y[n];

  double xfine[n*fine_i];
  double yspline[n*fine_i];

  double yinteg[n];
  double yintegspline[n*fine_i];

  double yderiv[n];
  double yderivspline[n*fine_i];

  for (int i = 0; i < n; i++) {
    x[i] = i;
    y[i] = x[i]*x[i];
    yinteg[i] = x[i]*x[i]*x[i]/3.0;
    yderiv[i] = 2.0*x[i];
    printf("%g %g %g %g\n", x[i], y[i], yinteg[i], yderiv[i]);
  }

  printf("\n\n");

  cspline * s = cspline_alloc(n, x, y);

  for (int i = 0; i < n*fine_i; i++) {
    xfine[i] = i/fine_d;

    if (xfine[i] <= x[n-1]){
      yspline[i] = cspline_eval(s, xfine[i]);
      yintegspline[i] = cspline_integ(s, xfine[i]);
      yderivspline[i] = cspline_deriv(s, xfine[i]);
      printf("%g %g %g %g\n", xfine[i], yspline[i], \
      yintegspline[i], yderivspline[i]);
    }

  }

  printf("\n\n");
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, x, y, n);
  gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
  gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline2, x, y, n);
  gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();
  gsl_spline *spline3 = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline3, x, y, n);

  for (double xi = x[0]; xi < x[n-1]; xi+= 0.1){
    double y1 = gsl_spline_eval (spline, xi, acc);
    double y2 = gsl_spline_eval_integ (spline2, 0, xi, acc2);
    double y3 = gsl_spline_eval_deriv (spline3, xi, acc3);
    printf ("%g %g %g %g\n", xi, y1, y2, y3);
  }

  cspline_free(s);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  gsl_spline_free (spline2);
  gsl_interp_accel_free (acc2);
  gsl_spline_free (spline3);
  gsl_interp_accel_free (acc3);

  return 0;
}
