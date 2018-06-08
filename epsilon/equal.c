int equal(double a, double b, double tau, double epsilon){

if (fabs(a-b) < tau){
        printf("Smaller than tau, 1\n");
        return 1;
        }
else if (fabs(a-b)/(fabs(a)+fabs(b)) < epsilon/2){
        printf("Smaller than epsilon, 1\n");
        return 1;
        }
else {
        printf("Nope... 0\n");
        return 0;
        }
}
