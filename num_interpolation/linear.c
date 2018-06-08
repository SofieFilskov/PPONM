#include<assert.h>
#include <stdio.h>

double linterp(int n, double* x, double* y, double z){
  assert ((n > 1) & (z >= x[0]) & (z<= x[n-1]));
  int i = 0;
  int j = n-1;
  while (j-i > 1){
    int m = (i+j)/2;
    if (z > x[m])
      i = m;
    else
      j = m;
  }
  return y[i] + (y[i+1] - y[i]) / (x[i+1] - x[i]) * (z - x[i]);
}

double linterp_integ(int n, double* x, double* y, double z){
  assert ((n > 1) & (z >= x[0]) & (z<= x[n-1]));
  int i = 0;
  int j = n-1;
  while (j-i > 1){
    int m = (i+j)/2;
    if (z > x[m])
      i = m;
    else
      j = m;
  }
  double integ = 0;
  for (int k = 0; k < i; k++){
    integ += y[k]*(x[k+1] - x[k]) + 0.5*(y[k+1] - y[k]) / (x[k+1] - x[k])\
      * (x[k+1] - x[k])*(x[k+1] - x[k]);
  }
  integ += y[i]*(z - x[i]) + 0.5*(y[i+1] - y[i]) / (x[i+1] - x[i])\
    * (z - x[i])*(z - x[i]);
  return integ;
}

int main() {
  int n = 5;
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

  for (int i = 0; i < n; i++) {
    x[i] = i;
    y[i] = x[i]*x[i];
    yinteg[i] = x[i]*x[i]*x[i]/3;
    printf("%g %g %g\n", x[i], y[i], yinteg[i]);
  }

  printf("\n\n");

  for (int i = 0; i < n*fine_i; i++) {
    xfine[i] = i/fine_d;
    yfine[i] = xfine[i]*xfine[i];
    yintegfine[i] = xfine[i]*xfine[i]*xfine[i]/3;

    if (xfine[i] <= x[n-1]){
      yspline[i] = linterp(n, x, y, xfine[i]);
      yintegspline[i] = linterp_integ(n, x, y, xfine[i]);
      printf("%g %g %g %g %g\n", xfine[i], yfine[i], yspline[i], yintegfine[i], yintegspline[i]);
    }
  }

  return 0;
}
