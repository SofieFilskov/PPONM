#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multimin.h>

typedef struct {int n; double (*f)(double, double); gsl_vector* data;} ann2D;

ann2D* ann2D_alloc(int n, double(*f)(double, double)) {
  ann2D* network = malloc(sizeof(ann2D));
  network->n = n;
  network->f = f;
  network->data = gsl_vector_alloc(5*n);
  return network;
}

void ann2D_free(ann2D* network){
  gsl_vector_free(network -> data);
  free(network);
}

double ann2D_feed_forward(ann2D* network, double x, double y){
  double s = 0;
  for (int i = 0; i<network -> n; i++){
    double a = gsl_vector_get(network -> data, 5*i+0);
    double b = gsl_vector_get(network -> data, 5*i+1);
    double c = gsl_vector_get(network -> data, 5*i+2);
    double d = gsl_vector_get(network -> data, 5*i+3);
    double w = gsl_vector_get(network -> data, 5*i+4);
    s += network -> f((x-a)/b, (y-c)/d)*w;
  }
  return s;
}


void ann2D_train(ann2D* network, gsl_vector* xlist, gsl_vector* ylist, gsl_matrix* zlist){
  int n = network -> data -> size;
  gsl_vector* p = gsl_vector_alloc(n);
  gsl_vector* startstep = gsl_vector_alloc(n);

  //The function we want to minimize
  double deltap(const gsl_vector* p, void* params){
    gsl_vector_memcpy(network -> data, p);
    double s = 0;
    for (int i = 0; i < xlist -> size; i++){
      for (int j = 0; j < ylist -> size; j++){
        double x = gsl_vector_get(xlist, i);
        double y = gsl_vector_get(ylist, i);
        double f = gsl_matrix_get(zlist, i, j);
        double z = ann2D_feed_forward(network, x, y);
        s += pow(f-z,2);
      }
    }
    return s;
  }

  gsl_vector_memcpy(p, network -> data);

  gsl_multimin_function F = {.f = deltap, .n = n, .params = NULL};
  gsl_vector_set_all(startstep, 0.1);

  gsl_multimin_fminimizer * minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,n);
  gsl_multimin_fminimizer_set(minimizer, &F, p, startstep);

  int iters = 0;
  int flag;
  do{
    iters++;
    flag = gsl_multimin_fminimizer_iterate(minimizer);
    if (flag) {
      break;
    }
    if (minimizer -> size <1e-4) {
      break;
    }
  } while(iters < 1e6);

  gsl_vector_memcpy(network -> data, minimizer -> x);

  gsl_vector_free(p);
  gsl_vector_free(startstep);
	gsl_multimin_fminimizer_free(minimizer);
}
