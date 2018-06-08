#include<math.h>
#include<stdio.h>
#include<stdlib.h>
double (*f)(double);
void print_h_of_1 (double (*h)(double)){
	printf("h of 1 = %g\n", h(1));
	}

double* foo(){static  double a[] = {0,0,0}; return a; }

int main(){

f = &sin;
printf("%g\n", (*f)(1) );
print_h_of_1 (&sin);

double (*g[3])(double) = {sin,cos,tan};
	g[2]=exp;
        for(int i=0;i<3;i++)printf("%g\n",g[i](1));

struct vector {double x,y,z;};
struct vector v = {1,2,3};
struct vector x = {1,2,3};
struct vector u = {.x=1,.y=2,.z=3};
struct vector w = {.z=3,.y=2,.x=1};
struct vector y;
typedef struct vector vector;
vector z;

printf("w.x = %g, v.x = %g\n", w.x, v.x);

  double* a=foo(); a[2]=1;
  double* b=foo(); b[2]=2;
  printf("a[2] = %g\n",a[2]);

struct fe {int n; double (*fe)(double); };
double (*fe)(double) = sin;
double yx = (*fe)(1);
printf("yx = %g\n", yx);

return 0;
}
