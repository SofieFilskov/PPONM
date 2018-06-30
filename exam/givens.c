#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

void givens(gsl_matrix* A, gsl_matrix* M){
	int n=A->size1;
	gsl_matrix_set_identity(M);
	int p,q;
	for (p=0; p<n; p++){
		for (q = p+1; q<n; q++){

			double app = gsl_matrix_get(A, p, p);
			double aqq = gsl_matrix_get(A, q, q);
			double apq = gsl_matrix_get(A, p, q);
			double aqp = gsl_matrix_get(A, q, p);

			/*We do the Givens*/
			double theta = atan2(aqp-apq, app+aqq);
			double ct = cos(theta), st = sin(theta);
			double a = ct*app+st*aqp;
			double b = -st*app+ct*aqp;
			double c = -apq*st+aqq*ct;

      gsl_matrix_set(A, p, p, a);
      gsl_matrix_set(A, p, q, b);
      gsl_matrix_set(A, q, p, b);
      gsl_matrix_set(A, q, q, c);

      gsl_matrix_set(M, p, q, theta);

		}
	}
}
