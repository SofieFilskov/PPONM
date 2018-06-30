#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

typedef struct {int n; double (*f)(double); gsl_vector* data;} ann2D;
ann2D* ann2D_alloc(int n, double(*f)(double, double));
void ann2D_free(ann2D* network);
double ann2D_feed_forward(ann2D* network, double x, double y);
void ann2D_train(ann2D* network, gsl_vector* xlist, gsl_vector* ylist, gsl_matrix* zlist);

double activation_func(double x, double y){
  return x*exp(-x*x)*y*exp(-y*y);
}

double function_to_fit(double x, double y){
  return sin(sqrt(x*x + y*y));
}

int main(int argc, char const *argv[]) {
  int n = 10; //number of hidden neurons
  ann2D* network = ann2D_alloc(n, activation_func);
  double ax = 0; //start position of function to fit
  double ay = 0;
  double bx = 2*M_PI; //end position of function to fit
  double by = 2*M_PI;
  int num_points_x = 20;
  int num_points_y = 20;

  gsl_vector* x = gsl_vector_alloc(num_points_x);
  gsl_vector* y = gsl_vector_alloc(num_points_y);
  gsl_matrix* Z = gsl_matrix_alloc(num_points_x, num_points_y);

  //Making all our x's and y's. Like lin-space
  for (int i = 0; i < num_points_x; i++){
    double xi = ax+(bx-ax)*i/(num_points_x-1);
    gsl_vector_set(x, i, xi);
  }

  for (int j = 0; j < num_points_y; j++){
    double yj = ay+(by-ay)*j/(num_points_y-1);
    gsl_vector_set(y, j, yj);
    for (int i = 0; i < num_points_x; i++){
      double zij = function_to_fit(gsl_vector_get(x, i), yj);
      gsl_matrix_set(Z, i, j, zij);
    }
  }

  for (int i = 0; i < n; i++){
    gsl_vector_set(network -> data, 5*i+0, ax+(bx-ax)*i/(n-1));
    gsl_vector_set(network -> data, 5*i+1, bx);
    gsl_vector_set(network -> data, 5*i+2, ay+(by-ay)*i/(n-1));
    gsl_vector_set(network -> data, 5*i+3, by);
    gsl_vector_set(network -> data, 5*i+4, 1); //weighted to be the same
  }

  ann2D_train(network, x, y, Z);

  //Want to plot the interpolated function
  double stepx = (bx-ax)/100;
  double stepy = (by-ay)/100;
  for (double x_interp = ax; x_interp < bx; x_interp += stepx){
    for(double y_interp = ay; y_interp < by; y_interp += stepy){
      double fit = ann2D_feed_forward(network, x_interp, y_interp);
      printf("%g %g %g %g\n", x_interp, y_interp, fit, function_to_fit(x_interp, y_interp));
    }
  }

  ann2D_free(network);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(Z);

  return 0;
}
