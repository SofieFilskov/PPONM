#include<assert.h>
#include<stdlib.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline* qspline_alloc(int n,double* x, double* y) {
  qspline *s = (qspline*) malloc(sizeof(qspline));
  s->b = (double*) malloc((n-1)*sizeof(double));
  s->c = (double*) malloc((n-1)*sizeof(double));
  s->x = (double*) malloc(n*sizeof(double));
  s->y = (double*) malloc(n*sizeof(double));
  s->n = n;
  for (int i = 0; i<n; i++) {
    s->x[i] = x[i];
    s->y[i] = y[i];
  }
  int i;
  double p[n-1], h[n-1];
  for (i=0; i<n-1; i++) {
    h[i] = x[i+1]-x[i];
    p[i] = (y[i+1]-y[i])/h[i];
  }
  s->c[0] = 0;
  for (i=0; i<n-2; i++) {
    s->c[i+1] = (p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
  }
  s->c[n-2]/=2;
  for (i=n-3; i>=0; i--) {
    s->c[i] = (p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
  }
  for (i=0; i<n-1; i++) {
    s->b[i] = p[i]-s->c[i]*h[i];
  }
  return s;
}

double qspline_eval(qspline *s, double z) {
  assert(z>=s->x[0] && z<= s-> x[s->n-1]);
  int i, j=s->n-1;
  while (j-i>1){
    int m=(i+j)/2;
    if (z>s->x[m])
      i = m;
    else
      j = m;
  }
  double h = z-s->x[i];
  return s->y[i]+h*(s->b[i]+h*s->c[i]);
}

double qspline_integ(qspline *s, double z) {
  assert(z>=s->x[0] && z<= s-> x[s->n-1]);
  int i = 0, j=s->n-1;
  while (j-i>1){
    int m=(i+j)/2;
    if (z>s->x[m]){
      i = m;
    }
    else {
      j = m;
    }
  }
  double integ = 0;
  for (int k = 0; k < i; k++){
    double h = s->x[k+1]-s->x[k];
    integ += s->y[k]*h + 0.5*s->b[k]*h*h + s->c[k]*h*h*h/3.0;
  }
  double h = z-s->x[i];
  integ += s->y[i]*h + 0.5*s->b[i]*h*h + s->c[i]*h*h*h/3.0;
  return integ;
}

double qspline_deriv(qspline *s, double z) {
  assert(s->n>1 && z>=s->x[0] && z<= s->x[s->n-1]);
  int i = 0, j=s->n-1;
  while (j-i>1){
    int m=(i+j)/2;
    if (z>s->x[m]){
      i = m;
    }
    else {
      j = m;
    }
  }
  return s->b[i]+2*(z-s->x[i])*s->c[i];
}

void qspline_free(qspline *s){
  free (s->x);
  free (s->y);
  free (s->b);
  free (s->c);
  free (s);
}
