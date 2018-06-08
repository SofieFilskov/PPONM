#include<complex.h>
#include<math.h>
#include<stdio.h>
int main()
{	
	printf("Exercise 1\n");
	int x = 5;
	printf("x = %i, gamma(x) = %g\n", x, tgamma(x));
	int y = 4*3*2;
	printf("Should be %i\n", y);

	double z = 0.5;
	printf("z = %g, bessel(z) = %g\n", z, j1(z));

	double complex a = csqrt(-2);
	printf("sqrt(-2) = %g + %g I\n", creal(a), cimag(a));

	double complex b = cexp(I * 1.0);
	printf("exp(i) = %g + %g I\n", creal(b), cimag(b));

	double complex c = cexp(I * M_PI);
	printf("exp(i * pi) = %g + %g I\n", creal(c), cimag(c));

	double complex d = cpow(I, M_E);
	double e = M_E;
	printf("e = %g\n", e);
	printf("I^e = %g + %g I\n", creal(d), cimag(d));
	
	printf("\nExercise 2\n");
	double dtal = 0.1111111111111111111111111111;
	long double ldtal = 0.1111111111111111111111111111L;
	float ftal = 0.1111111111111111111111111111;

	printf("double = %0.25lg \nlong double = %0.25Lg\n", dtal, ldtal);
	printf("float = %0.25f\n", ftal);

return 0;
}
