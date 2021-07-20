/* Testing for Fast_Acc_Sum */
#include "acc_sum.h"
#include <stdio.h>
/* Main Function */
int main(int argc,char **argv)
{
	double q[59] = {-92,1,3,2,9,4,4,10,11,8,8,8,18,22,8,8,17,21,16,
					16,16,16,16,16,16,16,16,16,16,16,16,16,43,42,86,32,32,32,32,
					32,32,32,32,7,32,32,32,32,32,32,32,32,32,32,32,32,32,32,50};
	double p[62];
	p[0] = pow(2,51) - 2060;
	p[60] = -(832+115712*eps)*eps;
	p[61] = -1/(4*eps) - 46.5;
	for(int i=1; i<60; i++)
	{
		p[i] = 1 - eps*q[i-1];
	}
	printf("naive sum = %.16f\n",sum(p,62));
	double s = fast_acc_sum(p,62);
	printf("sum = %.16f\n",s);
	
	double se = -2047.5 - (2058 + 115712*eps)*eps;
	printf("true sum = %.16f\n",se);
}