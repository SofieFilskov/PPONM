#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

int jacobi_B(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row, int* rot){
	int changed, sweeps = 0, n=A->size1;

	int p,q;
	for (p = row; p < row+1; p++){
		do {
			changed = 0; sweeps++;
			for (q = p+1; q<n; q++){
				*rot += 1;
				double app = gsl_vector_get(e,p);
				double aqq = gsl_vector_get(e,q);
				double apq = gsl_matrix_get(A, p, q);
				double phi = 0.5*atan2(2.0*apq, aqq-app);
				double c = cos(phi), s = sin(phi);
				double app1 = c*c*app-2.0*s*c*apq+s*s*aqq;
				double aqq1 = s*s*app+2.0*s*c*apq+c*c*aqq;
				if(app1 != app || aqq1 != aqq){
					changed = 1;
					gsl_vector_set(e, p, app1);
					gsl_vector_set(e, q, aqq1);
					gsl_matrix_set(A, p, q, 0.0);
					for (int i = 0; i<p; i++){
						double aip = gsl_matrix_get(A, i, p);
						double aiq = gsl_matrix_get(A, i, q);
						gsl_matrix_set(A, i, p, c*aip-s*aiq);
						gsl_matrix_set(A, i, q, c*aiq+s*aip);
					}
					for(int i = p+1; i<q; i++){
						double api = gsl_matrix_get(A, p, i);
						double aiq = gsl_matrix_get(A, i, q);
						gsl_matrix_set(A, p, i, c*api-s*aiq);
						gsl_matrix_set(A, i, q, c*aiq+s*api);
					}
					for(int i = q+1; i<n; i++){
						double api = gsl_matrix_get(A, p, i);
						double aqi = gsl_matrix_get(A, q, i);
						gsl_matrix_set(A, p, i, c*api-s*aqi);
						gsl_matrix_set(A, q, i, c*aqi+s*api);
					}
					for(int i = 0; i<n; i++){
						double vip = gsl_matrix_get(V, i, p);
						double viq = gsl_matrix_get(V, i, q);
						gsl_matrix_set(V, i, p, c*vip-s*viq);
						gsl_matrix_set(V, i, q, c*viq+s*vip);
					}
				}
			}
		}
		while(changed!=0);
	}
	return sweeps;
}

int jacobi_B_high(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, int row, int* rot){
	int changed, sweeps = 0, n=A->size1;

	int p,q;
	for (p = row; p < row+1; p++){
		do {
			changed = 0; sweeps++;
			for (q = p+1; q<n; q++){
				*rot += 1;
				double app = gsl_vector_get(e,p);
				double aqq = gsl_vector_get(e,q);
				double apq = gsl_matrix_get(A, p, q);
				double phi = 0.5*atan2(2.0*apq, aqq-app)+M_PI*0.5;
				double c = cos(phi), s = sin(phi);
				double app1 = c*c*app-2.0*s*c*apq+s*s*aqq;
				double aqq1 = s*s*app+2.0*s*c*apq+c*c*aqq;
				if(app1 != app || aqq1 != aqq){
					changed = 1;
					gsl_vector_set(e, p, app1);
					gsl_vector_set(e, q, aqq1);
					gsl_matrix_set(A, p, q, 0.0);
					for (int i = 0; i<p; i++){
						double aip = gsl_matrix_get(A, i, p);
						double aiq = gsl_matrix_get(A, i, q);
						gsl_matrix_set(A, i, p, c*aip-s*aiq);
						gsl_matrix_set(A, i, q, c*aiq+s*aip);
					}
					for(int i = p+1; i<q; i++){
						double api = gsl_matrix_get(A, p, i);
						double aiq = gsl_matrix_get(A, i, q);
						gsl_matrix_set(A, p, i, c*api-s*aiq);
						gsl_matrix_set(A, i, q, c*aiq+s*api);
					}
					for(int i = q+1; i<n; i++){
						double api = gsl_matrix_get(A, p, i);
						double aqi = gsl_matrix_get(A, q, i);
						gsl_matrix_set(A, p, i, c*api-s*aqi);
						gsl_matrix_set(A, q, i, c*aqi+s*api);
					}
					for(int i = 0; i<n; i++){
						double vip = gsl_matrix_get(V, i, p);
						double viq = gsl_matrix_get(V, i, q);
						gsl_matrix_set(V, i, p, c*vip-s*viq);
						gsl_matrix_set(V, i, q, c*viq+s*vip);
					}
				}
			}
		}
		while(changed!=0);
	}
	return sweeps;
}
