#ifndef HAVE_NVECTOR_H

typedef struct {int size; double* data;} nvector;

nvector*	nvector_alloc		(int n);
void		nvector_free 		(nvector* v);
void		nvector_set		(nvector* v, int i, double value);
double 		nvector_get		(nvector* v, int i);
double 		nvector_dot_product	(nvector* u, nvector* v);
void		nvector_add		(nvector* a, nvector* b);
int 		nvector_equal		(nvector* a, nvector* b);

#define HAVE_NVECTOR_H
#endif
