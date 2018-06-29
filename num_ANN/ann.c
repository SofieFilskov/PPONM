#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_multimin.h>

typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;

ann* ann_alloc(int n, double(*f)(double)) {
  ann* network = malloc(sizeof(ann));
  network->n = n;
  network->f = f;
  network->data = gsl_vector_alloc(3*n);
  return network;
}

void ann_free(ann* network){
  gsl_vector_free(network -> data);
  free(network);
}

double ann_feed_forward(ann* network, double x){
  double s = 0;
  for (int i = 0; i<network -> n; i++){
    double a = gsl_vector_get(network -> data, 3*i+0); //data[0, 3, 6]
    double b = gsl_vector_get(network -> data, 3*i+1); //data[1, 4, 7]
    double w = gsl_vector_get(network -> data, 3*i+2); //data[2, 5, 8]
    s += network -> f((x-a)/b)*w;
  }
  return s;
}


void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist){
  int n = network -> data -> size;
  gsl_vector* p = gsl_vector_alloc(n);
  gsl_vector* startstep = gsl_vector_alloc(n);

  //The function we want to minimize
  double deltap(const gsl_vector* p, void* params){
    gsl_vector_memcpy(network -> data, p);
    double s = 0;
    for (int i = 0; i < xlist -> size; i++){
      double x = gsl_vector_get(xlist, i);
      double f = gsl_vector_get(ylist, i);
      double y = ann_feed_forward(network, x);
      s += pow(f-y,2);
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
