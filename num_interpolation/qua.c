#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef struct {int n; double * x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double * x, double * y);
double qspline_eval(qspline * s, double z);
double qspline_integ(qspline * s, double z);
double qspline_deriv(qspline * s, double z);
void qspline_free(qspline * s);

int main() {
  int n = 10;
  int fine_i = 10;
  double fine_d = 10.0;

  double x[n];
  double y[n];

  double xfine[n*fine_i];
  double yfine[n*fine_i];
  double yspline[n*fine_i];

  double yinteg[n];
  double yintegfine[n*fine_i];
  double yintegspline[n*fine_i];

  double yderiv[n];
  double yderivfine[n*fine_i];
  double yderivspline[n*fine_i];

  for (int i = 0; i < n; i++) {
    x[i] = i;
    y[i] = x[i]*x[i];
    yinteg[i] = x[i]*x[i]*x[i]/3.0;
    yderiv[i] = 2.0*x[i];
    printf("%g %g %g %g\n", x[i], y[i], yinteg[i], yderiv[i]);
  }

  printf("\n\n");

  qspline * s = qspline_alloc(n, x, y);

  for (int i = 0; i < n*fine_i; i++) {
    xfine[i] = i/fine_d;
    yfine[i] = xfine[i]*xfine[i];
    yintegfine[i] = xfine[i]*xfine[i]*xfine[i]/3.0;
    yderivfine[i] = 2.0*xfine[i];

    if (xfine[i] <= x[n-1]){
      yspline[i] = qspline_eval(s, xfine[i]);
      yintegspline[i] = qspline_integ(s, xfine[i]);
      yderivspline[i] = qspline_deriv(s, xfine[i]);
      printf("%g %g %g %g %g %g %g\n", xfine[i], yfine[i], yspline[i], \
      yintegfine[i], yintegspline[i], yderivfine[i], yderivspline[i]);
    }

  }

  qspline_free(s);

  return 0;
}
