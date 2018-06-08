#include<stdio.h>
#include"komplex.h"
//#define TINY 1e-6

int main(){
	komplex a = {1, 2};
	komplex b = {3, 4};

	printf("testing komplex_add\n");
	komplex r = komplex_add(a, b);
	komplex R = {4, 6};
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("a+b should be = ", R);
	komplex_print("a+b is calculated to = ", r);

	printf("testing komplex_sub\n");
	komplex s = komplex_sub(b, a);
	komplex S = {2, 2};
	komplex_print("b-a should be = ", S);
	komplex_print("b-a is calculated to = ", s);

}

