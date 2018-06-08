#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>

typedef struct {int n; double *t,*y,*e;} experimental_data;

int least_square(const gsl_vector * current_point, void *params){
  //f(t) = A*exp(-t/T)+B
  double A = gsl_vector_get(current_point,0);
  double T = gsl_vector_get(current_point,1);
  double B = gsl_vector_get(current_point,2);

  experimental_data *par = (experimental_data*) params;
	int     n = par->n;
	double *t = par->t;
	double *y = par->y;
	double *e = par->e;

  double sum = 0;
  double f(double t) {return A*exp(-t/T)+B;} //making function f(t)
  for (int i = 0; i<n; i++){
    sum += pow((f(t[i])-y[i])/e[i],2);
  }

  return sum; //as it is the sum we try to minimize
}

int main(){
  double t[]= {0.47,1.41,2.36,3.30,4.24,5.18,6.13,7.07,8.01,8.95};
  double y[]= {5.49,4.08,3.54,2.61,2.09,1.91,1.55,1.47,1.45,1.25};
  double e[]= {0.26,0.12,0.27,0.10,0.15,0.11,0.13,0.07,0.15,0.09};
  int n = sizeof(t)/sizeof(t[0]);

  experimental_data data;
  data.n=n;
  data.t=t;
  data.y=y;
  data.e=e;

  const gsl_multimin_fminimizer_type *type;
  gsl_multimin_fminimizer *s;

  const int dim=3;

  gsl_multimin_function F; //Giving the function
  F.f = least_square;
  F.n = dim;
  F.params = (void*) &data;

  type = gsl_multimin_fminimizer_nmsimplex2;
  s = gsl_multimin_fminimizer_alloc(type, dim);

  size_t iter = 0;
  int status;

  /* Starting point */
  gsl_vector* start = gsl_vector_alloc (dim);
  gsl_vector_set (start, 0, 2.0);
  gsl_vector_set (start, 1, 2.0);
  gsl_vector_set (start, 2, 2.0);

  /* Set initial step sizes to 0.1 */
  gsl_vector* step = gsl_vector_alloc (dim);
  gsl_vector_set_all (step, 0.1);

  gsl_multimin_fminimizer_set (s, &F, start, step);

  do{
    iter++;

    status = gsl_multimin_fminimizer_iterate(s);
    if (status)
      break;

    status = gsl_multimin_test_size(s->size,1e-4);
  }
  while (status == GSL_CONTINUE && iter < 1000);

  double A=gsl_vector_get(s->x,0);
  double T=gsl_vector_get(s->x,1);
  double B=gsl_vector_get(s->x,2);

  double f(double t) {return A*exp(-(t)/T) + B;}

  for(double x=t[0];x<t[n-1];x+=0.05){
    printf("%g %g\n",x,f(x));
  }
  printf("\n\n");
  printf("A=%g T=%g B=%g iter=%i\n",A,T,B,iter);

  gsl_vector_free(start);
  gsl_vector_free(step);
  gsl_multimin_fminimizer_free(s);

  return 0;
}
