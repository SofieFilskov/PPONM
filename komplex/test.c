#include <stdio.h>
//void f(int* p){*p = 4;}
/*double func(){
        double x[4];
        x[0] = 0;
        return x[0];}
*/
int k=2; /* file scope */
void fun(){printf("k=%i\n",k);}

int main(){

/*
int i = 1;
f(&i);
printf("i=%i\n",i);
*/

//func();
//printf("x[0] = %g\n", x[0]);

int k=1; /* function scope */
	{
		int k=0; /* block scope */
		printf("k=%i\n",k);
	}
	printf("k=%i\n",k);
	fun();

return 0;
}



