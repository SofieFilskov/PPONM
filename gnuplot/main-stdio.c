#include<stdio.h>
#include<math.h>

int main(){

double x;
while( scanf("%lg", &x) != EOF ) {
	printf("%lg \t %lg\n", x, cos(x));
	//printf("%lg\n",x); 
	}

return 0;
}


/* write
echo 1 2 3 4 5 | ./main-stdio > test.io.out.txt
in command line in terminal to run*/
