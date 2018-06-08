#include<assert.h>
#include<stdlib.h>
#include<stdio.h>

typedef struct {int n; double *x, *y, *b, *c, *d;} cspline;

cspline* cspline_alloc(int n,double* x, double* y) {
  cspline *s = (cspline*) malloc(sizeof(cspline));
  s->b = (double*) malloc(n*sizeof(double));
  s->c = (double*) malloc((n-1)*sizeof(double));
  s->d = (double*) malloc((n-1)*sizeof(double));
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
    assert(h[i] > 0);
    p[i] = (y[i+1]-y[i])/h[i];
  }
  double D[n], Q[n-1], B[n];
  D[0] = 2;
  Q[0] = 1;
  B[0] = 3*p[0];
  for(i = 0; i<n-2;i++){
    D[i+1] = 2*h[i]/h[i+1]+2;
    D[n-1] = 2;
    Q[i+1] = h[i]/h[i+1];
    B[i+1] = 3*(p[i]+p[i+1]*h[i]/h[i+1]);
  }
  B[n-1] = 3*p[n-2];

  for(i=1; i<n;i++){
    D[i]-=Q[i-1]/D[i-1];
    B[i]-=B[i-1]/D[i-1];
  }
  s->b[n-1] = B[n-1]/D[n-1];
  for(i = n-2; i>=0; i--) {
    s->b[i] = (B[i]-Q[i]*s->b[i+1])/D[i];
  }
  for (i = 0; i<n-1;i++){
    s->c[i] = (-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
    s->d[i] = (s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
  }
  return s;
}

double cspline_eval(cspline *s, double z) {
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
  return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}

double cspline_integ(cspline *s, double z) {
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
    integ += s->y[k]*h + 0.5*s->b[k]*h*h + s->c[k]*h*h*h/3.0 + s->d[k]*h*h*h*h/4.0;
  }
  double h = z-s->x[i];
  integ += s->y[i]*h + 0.5*s->b[i]*h*h + s->c[i]*h*h*h/3.0 + s->d[i]*h*h*h*h/4.0;
  return integ;
}

double cspline_deriv(cspline *s, double z) {
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
  return s->b[i]+2*(z-s->x[i])*s->c[i]+3*(z-s->x[i])*(z-s->x[i])*s->d[i];
}

void cspline_free(cspline *s){
  free (s->x);
  free (s->y);
  free (s->b);
  free (s->c);
  free (s->d);
  free (s);
}
