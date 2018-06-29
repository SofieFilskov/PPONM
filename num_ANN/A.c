#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>

typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;
ann* ann_alloc(int n, double(*f)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist);

double activation_func(double x){
  return x*exp(-x*x);
}

double function_to_fit(double x){
  return sin(x);
}

int main(int argc, char const *argv[]) {
  int n = 10; //number of hidden neurons
  ann* network = ann_alloc(n, activation_func);
  double a = 0; //start position of function to fit
  double b = 2*M_PI; //end position of function to fit
  int num_points = 20;

  gsl_vector* x = gsl_vector_alloc(num_points);
  gsl_vector* y = gsl_vector_alloc(num_points);

  //Making all our x's and y's. Like lin-space
  for (int i = 0; i < num_points; i++){
    double xi = a+(b-a)*i/(num_points-1);
    double yi = function_to_fit(xi);
    gsl_vector_set(x, i, xi);
    gsl_vector_set(y, i, yi);

    printf("%g %g\n", xi, yi);
  }

  for (int i = 0; i < n; i++){
    gsl_vector_set(network -> data, 3*i+0, a+(b-a)*i/(n-1)); //a_i
    gsl_vector_set(network -> data, 3*i+1, b); //b_i
    gsl_vector_set(network -> data, 3*i+2, 1); //weighted to be the same
  }

  ann_train(network, x, y);

  printf("\n\n");

  //Want to plot the interpolated function
  double step = (b-a)/100;
  for (double x_interp = a; x_interp < b; x_interp += step){
    double fit = ann_feed_forward(network, x_interp);
    printf("%g %g\n", x_interp, fit);
  }

  ann_free(network);
  return 0;
}
