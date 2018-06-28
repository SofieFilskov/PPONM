#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>
#include <string.h>

typedef struct filename {char name[30];} fname;

fname making_name(char * name){
  fname info;
  strcpy(info.name, name);
  return info;
}

int driver(double* t, double b, double* h, gsl_vector* yt, double acc, double eps,
  void rkstep12(double t, double h, gsl_vector* yt,
    void f(double t, gsl_vector* y, gsl_vector* dydt),
  	gsl_vector* yth, gsl_vector* err),
	void f(double t, gsl_vector* y, gsl_vector* dydt),
  fname* name);

void sinus(double t, gsl_vector* y, gsl_vector* dydt) {
  gsl_vector_set(dydt,0, gsl_vector_get(y,1));
  gsl_vector_set(dydt,1, -1*gsl_vector_get(y,0));
}

void function(double t, gsl_vector* y, gsl_vector* dydt) {
  /* x^2 - 5x +7 = f(x) */
  gsl_vector_set(dydt, 0, gsl_vector_get(y,1));
  gsl_vector_set(dydt, 1, 2);
}

void rkstep12(
  double t, double h, gsl_vector* yt,
  void f(double t, gsl_vector* y, gsl_vector* dydt),
  gsl_vector* yth, gsl_vector* err);

int main(int argc, char const *argv[]) {
  /* Solving differential equation for sine function */
  int n = 2;
  double t = 0;
  double b = 6.28;
  double h = 0.02;
  double acc = 0.0005;
  double eps = 0.001;

  gsl_vector* yt = gsl_vector_alloc(n);
  /* start conditions */

  gsl_vector_set(yt, 0, 0); /*sin(0) = 0*/
  gsl_vector_set(yt, 1, 1); /*cos(0) = 1*/

  fname sin_name = making_name("sinus.out.txt");
  int k = driver(&t, b, &h, yt, acc, eps, rkstep12, sinus, &sin_name);

  printf("Solving sine from 0 to %g\n", b);
  printf("Final value, y=%g dy/dx=%g\n", gsl_vector_get(yt,0),gsl_vector_get(yt,1));
  printf("Math values: sin(b)=%g cos(b)=%g\n",sin(b),cos(b));
  printf("Steps taken: %i\n", k);

  printf("\n\n");

  gsl_vector_free(yt);

  /* Solving differential equation for 2. degree polynomial */

  n = 2;
  t = 0;
  b = 7;
  h = 0.02;
  acc = 0.0005;
  eps = 0.001;

  gsl_vector* yt2 = gsl_vector_alloc(n);
  /* start conditions */
  gsl_vector_set(yt2, 0, 7);
  gsl_vector_set(yt2, 1, -5);
  fname func_name = making_name("function.out.txt");
  int k2 = driver(&t, b, &h, yt2, acc, eps, rkstep12, function, &func_name);
  printf("Solving function from 0 to %g\n", b);
  printf("Final value, y=%g\n", gsl_vector_get(yt2,0));
  printf("Math values: x^2 - 5 x + 7 =%g\n",b*b-5*b+7);
  printf("Steps taken: %i\n", k2);

  gsl_vector_free(yt2);

  return 0;
}
