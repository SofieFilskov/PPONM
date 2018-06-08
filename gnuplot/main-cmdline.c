#include<math.h>
#include<stdlib.h>
#include<stdio.h>

int main(int argc, char** argv) {

for(int i = 1; i<argc; i++) {
	double x = atof(argv[i]);
	printf("%lg \t %lg\n", x, sin(x));
	}

/*write
./main-cmdline 1 2 3 4 5 6 > test.cmd.out.txt
in commad line */

return 0;
}
