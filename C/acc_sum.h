#ifndef ACC_SUM
#define ACC_SUM
#include "eft.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
/* Global Constants */
double eps = DBL_EPSILON/2;
double eta = DBL_TRUE_MIN;
/* Unit in First Place */
double ufp(const double p)
{
	double q = p/(2*eps) + p;
	return fabs(q - (1 - eps)*q);
}
/* Two Sum */
void two_sum(const double a,const double b,struct eft* res)
{
	res->fl_res = a + b;
	double t = res->fl_res - a;
	res->fl_err = (a - (res->fl_res - t)) + (b - t);
}
/* Error Free Array Extraction */
double extract(double* p,double sigma,const unsigned int n)
{
	struct eft res;
	for(int i=0; i<n; ++i)
	{
		two_sum(sigma,p[i],&res);
		sigma = res.fl_res;
		p[i] = res.fl_err;
	}
	return sigma;
}
/* Sum */
double sum(const double* p,const unsigned int n)
{
	double s = p[0];
	for(int i=1; i<n; ++i)
	{
		s += p[i];
	}
	return s;
}
/* Absolute Sum */
double abs_sum(const double* p,const unsigned int n)
{
	double s = fabs(p[0]);
	for(int i=1; i<n; ++i)
	{
		s += fabs(p[i]);
	}
	return s;
}
/* Fast Accurate Summation */
double fast_acc_sum(double* p,const unsigned int n)
{
	// variables
	double phi, sigma, sigma_new, T, t, t_new, tau, u;
	// goto label
	start:
	T = abs_sum(p,n)/fma(-n,eps,1);
	if(T<=eta/eps)
	{
		return sum(p,n);	// no rounding error
	}
	t = 0;
	do {
		sigma = (2*T)/(1 - fma(3,n,1)*eps);
		sigma_new = extract(p,sigma,n);
		tau = sigma_new - sigma;
		t_new = t;
		t = t_new + tau;
		if(t==0)
		{
			goto start;		// intermediate sum is zero, recursively apply fast_acc_sum to array of lower order parts
		}
		u = ufp(sigma);
		phi = ((2*n*(n+2)*eps)*u)/fma(-5,eps,1);
		T = fmin((fma(4,eps,1.5)*(n*eps))*sigma,(2*n*eps)*u);
	} while((fabs(t) < phi) && (4*T > eta/eps));
	tau = (t_new - t) + tau;
	// return
	return t + (tau + sum(p,n));
}
/* Fast Complex Accurate Summation */
double complex fast_cmplx_acc_sum(double complex* p,const unsigned int n)
{
	// variables
	double* realp = (double*)malloc(n*sizeof(double));
	double* imagp = (double*)malloc(n*sizeof(double));
	for(int i=0; i<n; i++)
	{
		realp[i] = creal(p[i]);
		imagp[i] = cimag(p[i]);
	}
	// summations
	double real = fast_acc_sum(realp,n);
	double imag = fast_acc_sum(imagp,n);
	// free
	free(realp);
	free(imagp);
	// return
	return real+I*imag;
}
#endif