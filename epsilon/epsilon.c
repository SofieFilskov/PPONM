#include<limits.h>
#include<float.h>
#include<math.h>
#include<stdio.h>
#include "equal.c"
int equal(double a, double b, double tau, double epsilon);

int main()
{
//	Number 1	//
//MAXIMUM
printf("1. FINDING MAXIMUM\n");
int i = 1;
while(i+1>i)
	{i++;}
printf("Using while my max is %i\n", i);
printf("The number max is %i\n", INT_MAX);

int j = 1;
for(j; j+1>j; j++){}
printf("Using for my max is %i\n", j);
printf("The number max is %i\n", INT_MAX);

int h = 1;
do
	{h++;}
	while(h+1>h);
printf("Using do my max is %i\n", h);
printf("The number max is %i\n", INT_MAX);


//MINIMUM
printf("\n1. FINDING MINIMUM\n");
int x = 1;
while(x>x-1)
        {x--;}
printf("Using while my min is %i\n", x);
printf("The number min is %i\n", INT_MIN);

int y = 1;
for(y; y>y-1; y--){}
printf("Using for my min is %i\n", y);
printf("The number min is %i\n", INT_MIN);

int z = 1;
do
        {z--;}
        while(z>z-1);
printf("Using do my min is %i\n", z);
printf("The number min is %i\n", INT_MIN);


//MACHINE EPSILON
printf("\n1. MACHINE EPSILON\n");

double a = 1;
long double b = 1;
float c = 1;

while(1+a!=1) //this means while 1+a not equal to 1
	{a/=2;} //this means a = a/2
a*=2; //this means a = a*2
printf("***Using while***\n");
printf("For double my epsilon is = %g\n", a);
printf("Double epsilon = %g\n", DBL_EPSILON);

while(1+b!=1)
        {b/=2;}
b*=2;
printf("For long double my epsilon is = %Lg\n", b);
printf("Long double epsilon = %Lg\n", LDBL_EPSILON);

while(1+c!=1)
        {c/=2;}
c*=2;
printf("For float my epsilon is  = %g\n", c);
printf("Float epsilon = %g\n", FLT_EPSILON);

//Now we do it with for
double d = 1;
long double f = 1;
float g = 1;
printf("\n***Using for***\n");
for(d; 1+d!=1; d/=2){}
d*=2;
printf("For double my epsilon is = %g\n", d);
printf("Double epsilon = %g\n", DBL_EPSILON);

for(f; 1+f!=1; f/=2){}
f*=2;
printf("For long double my epsilon is = %Lg\n", f);
printf("Long double epsilon = %Lg\n", LDBL_EPSILON);

for(g; 1+g!=1; g/=2){}
g*=2;
printf("For float my epsilon is  = %g\n", g);
printf("Float epsilon = %g\n", FLT_EPSILON);

//Now with do
double k = 1;
long double l = 1;
float m = 1;
printf("\n***Using do***\n");
do
{k/=2;}
while(1+k!=1);
k*=2;
printf("For double my epsilon is = %g\n", k);
printf("Double epsilon = %g\n", DBL_EPSILON);

do
{l/=2;}
while(1+l!=1);
l*=2;
printf("For double my epsilon is = %g\n", l);
printf("Double epsilon = %g\n", DBL_EPSILON);

do
{m/=2;}
while(1+m!=1);
m*=2;
printf("For double my epsilon is = %g\n", m);
printf("Double epsilon = %g\n", DBL_EPSILON);


//	Number 2	//
printf("\n2. ITERATION\n");
int max = INT_MAX/3;
float sum_up_float = 0.0;

for(int l = 1; l<= max; l++)
	{
	sum_up_float += 1.0f/l;
	}
printf("sum_up_float = %g\n", sum_up_float);

float sum_down_float = 0.0;
int s = max;
while(s>=1){
	sum_down_float += 1.0f/s;
	s--;
	}
printf("sum_down_float = %g\n", sum_down_float);

printf("The difference must have something to do with rounding.\n");
printf("It does converge, as it cannot handle so many decimals\n");

//now with double
double sum_up_dbl = 0.0;

for(int n = 1; n<= max; n++)
        {
        sum_up_dbl += 1.0f/n;
        }
printf("sum_up_double = %g\n", sum_up_dbl);

int o = max;
double sum_down_dbl = 0.0;
while(o>=1){
        sum_down_dbl += 1.0f/o;
        o--;
        }
printf("sum_down_double = %g\n", sum_down_dbl);
printf("Double has more decimals and is therefore more precise\n");


//	EQUAL	//
printf("\n3. EQUAL\n");

double p = 12.0;
double q = 12.0;
double tau = 9.0;
double epsilon = 0.0003;
equal(p, q, tau, epsilon);


return 0;
}
