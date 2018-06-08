#include<stdio.h>
#include"nvector.h"

nvector* nvector_alloc(int n){
	nvector* v = malloc(sizeof(nvector));
	(*v).size = n;
	(*v).data = malloc(n*sizeof(double));
	if( v == NULL) printf(stderr, "error in nvector_alloc\n");
	return v;
	}

void nvector_free(nvector* v) {free(v->data); free(v);}

void nvector_set(nvector* v, int i, double value) { (*v).data[i] = value;}

double nvector_get(nvector* v, int i) {return (*v).data[i];}

void nvector_add(nvector* a, nvector* b) {
	int i = (*a).size;
	for(int j = 0; j < i; j++) {
		(*a).data[j] = (*a).data[j]+(*b).data[j];
		}
	}

double nvector_dot_product (nvector* u, nvector* v) {
	int i = (*u).size;
	double sum = 0;
        for(int j = 0; j < i; j++) {
                sum = sum + (*u).data[j]*(*v).data[j]; 
		}
	return sum;
	}

int nvector_equal (nvector* a, nvector* b) {
	double tau = 0.1;
	int sum = 0;
	int i = (*a).size;
	for(int j = 0; j<i; j++) {
		double c = (*a).data[j];
		double d = (*b).data[j];
		if (fabs(c-d)<tau) {sum++;}
		}
	if (sum == i) {return 1;}
	else {return 0;}
	}

//Also dot product, print, add, equal - try to play with it.
