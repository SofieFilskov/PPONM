#include<gsl/gsl_sf_airy.h>
#include<stdio.h>

int main() {
  gsl_mode_t prec = GSL_PREC_DOUBLE;
  for (double x = -2*3.1415; x < 2*3.1415; x += 0.010) {
    printf("%g\t%g\t%g\n", x, gsl_sf_airy_Ai(x, prec), gsl_sf_airy_Bi(x, prec));
  }

  return 0;
}
