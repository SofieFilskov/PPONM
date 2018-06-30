#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

int jacobi(gsl_matrix* A, gsl_vector* e, gsl_matrix* V, gsl_matrix* U, int* rot){
	int changed, sweeps = 0, n=A->size1;
	for(int i = 0; i < n; i++){
		gsl_vector_set(e, i, gsl_matrix_get(A, i, i));
	}
	gsl_matrix_set_identity(V);
	gsl_matrix_set_identity(U);
	do {
		changed = 0; sweeps++; int p,q;
		for (p=0; p<n; p++){
			for (q = p+1; q<n; q++){

				double app = gsl_vector_get(e,p);
				double aqq = gsl_vector_get(e,q);
				double apq = gsl_matrix_get(A, p, q);
				double aqp = gsl_matrix_get(A, q, p);

				/*We do the Givens*/
				double theta = atan2(aqp-apq, app+aqq);
				double ct = cos(theta), st = sin(theta);
				double a = ct*app+st*aqp;
				double b = -st*app+ct*aqp;
				double c = -apq*st+aqq*ct;

				app = a;
				aqq = c;
				apq = b;
				aqp = b;

				/*This is Jacobi*/
				double phi = 0.5*atan2(2*apq, aqq-app);
				double cp = cos(phi);
				double sp = sin(phi);
				double app1 = cp*cp*app-2*sp*cp*apq+sp*sp*aqq;
				double aqq1 = sp*sp*app+2*sp*cp*apq+cp*cp*aqq;
				if(app1 != app || aqq1 != aqq){
					/*Givens*/
					for (int i = 0; i<n; i++){
						double aip = gsl_matrix_get(A, p, i);
						double aiq = gsl_matrix_get(A, q, i);
						gsl_matrix_set(A, p, i, ct*aip+st*aiq);
						gsl_matrix_set(A, q, i, ct*aiq-st*aip);
					}

					/*Jacobi*/
					*rot += 1;
					changed = 1;
					gsl_vector_set(e, p, app1);
					gsl_vector_set(e, q, aqq1);
					gsl_matrix_set(A, p, q, 0.0);
					gsl_matrix_set(A, q, p, 0.0);
					for (int i = 0; i<p; i++){
						double aip = gsl_matrix_get(A, i, p);
						double aiq = gsl_matrix_get(A, i, q);
						gsl_matrix_set(A, i, p, cp*aip-sp*aiq);
						gsl_matrix_set(A, i, q, cp*aiq+sp*aip);
					}
					for(int i = p+1; i<q; i++){
						double api = gsl_matrix_get(A, p, i);
						double aiq = gsl_matrix_get(A, i, q);
						gsl_matrix_set(A, p, i, cp*api-sp*aiq);
						gsl_matrix_set(A, i, q, cp*aiq+sp*api);
					}
					for(int i = q+1; i<n; i++){
						double api = gsl_matrix_get(A, p, i);
						double aqi = gsl_matrix_get(A, q, i);
						gsl_matrix_set(A, p, i, cp*api-sp*aqi);
						gsl_matrix_set(A, q, i, cp*aqi+sp*api);
					}
					for(int i = 0; i<n; i++){
						double vip = gsl_matrix_get(V, i, p);
						double viq = gsl_matrix_get(V, i, q);
						gsl_matrix_set(V, i, p, cp*vip-sp*viq);
						gsl_matrix_set(V, i, q, cp*viq+sp*vip);
					}
					for(int i = 0; i<n; i++){
						double uip = gsl_matrix_get(U, i, p);
						double uiq = gsl_matrix_get(U, i, q);
						gsl_matrix_set(U, i, p, uip*(ct*cp+st*sp)+uiq*(st*cp-ct*sp));
						gsl_matrix_set(U, i, q, uip*(ct*sp-st*cp)+uiq*(st*sp+ct*cp));
					}
				}
			}
		}
	}
	while(changed!=0);
	return sweeps;
}
