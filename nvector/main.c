#include"nvector.h"
#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX

int main(){

int n = 5;

printf("\nmain: testing nvector_alloc \n");
nvector *v = nvector_alloc(n);
if (v == NULL) printf("test failed\n");
else printf("test passed\n");

printf("\nmain: testing nvector_set and nvector_get \n");
double value = RND;
int i = n/2;
nvector_set(v, i, value);
double v_i = nvector_get(v,i);
printf("the set value = %g \nthe extracted v_i = %g\n", value, v_i);
//here we can also use a nvector_equal, but then we need to make it.

printf("\nmain: testing nvector_add\n");
nvector *a = nvector_alloc(n);
nvector *b = nvector_alloc(n);
nvector *c = nvector_alloc(n);
for (int i = 0; i < n; i++) {
	double x = RND, y = RND;
	nvector_set(a, i, x);
	nvector_set(b, i, y);
	nvector_set(c, i, x + y);
	}
nvector_add(a, b);
printf("a+b should = %g\n", (*c).data[1]);
printf("a+b actually = %g\n", (*a).data[1]);

printf("\nmain: testing nvector_dot_product\n");
double d = nvector_dot_product(a,b);
printf("a * b = %g\n", d);

printf("\nmain: testing nvector_equal\n");
if (nvector_equal(c,a)) {printf("test passed\n");}
else {printf("test failed\n");}

nvector_free(a);
nvector_free(b);
nvector_free(c);
nvector_free(v);

return 0;
}
