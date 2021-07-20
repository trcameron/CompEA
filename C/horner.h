#ifndef HORNER
#define HORNER
#include "eft.h"
#include <float.h>
#include <math.h>
#include <mpfr.h>
#include <mpc.h>
/* Global Constants */
double EPS = DBL_EPSILON/2;
double ETA = DBL_TRUE_MIN;
/* Unit in First Place */
double ufp(const double p)
{
	double q = p/DBL_EPSILON + p;
	return fma(q,EPS-1,q);
}
/* Two Sum */
void two_sum(const double a,const double b,struct eft* res)
{
	res->fl_res = a + b;
	double t = res->fl_res - a;
	res->fl_err = (a - (res->fl_res - t)) + (b - t);
}
/* Two Product */
void two_prod(const double a,const double b,struct eft* res)
{
	res->fl_res = a*b;
	res->fl_err = fma(a,b,-res->fl_res);
}
/* Two Sum Cmplx */
void two_sum_cmplx(const double complex a,const double complex b,struct eft* eft_arr,struct eft_cmplx_sum* res)
{
	//struct eft real_eft, imag_eft;
	two_sum(creal(a),creal(b),eft_arr);
	two_sum(cimag(a),cimag(b),eft_arr+1);
	res->fl_res = eft_arr[0].fl_res + I*eft_arr[1].fl_res;
	res->fl_err = eft_arr[0].fl_err + I*eft_arr[1].fl_err;
}
/* Two Product Cmplx */
void two_prod_cmplx(const double complex a,const double complex b,struct eft* eft_arr,struct eft_cmplx_prod* res)
{
	//struct eft real_prod1, real_prod2, real_prod3, real_prod4, real_sum1, real_sum2;
	two_prod(creal(a),creal(b),eft_arr);
	two_prod(cimag(a),cimag(b),eft_arr+1);
	two_prod(creal(a),cimag(b),eft_arr+2);
	two_prod(cimag(a),creal(b),eft_arr+3);
	two_sum(eft_arr[0].fl_res,-eft_arr[1].fl_res,eft_arr+4);
	two_sum(eft_arr[2].fl_res,eft_arr[3].fl_res,eft_arr+5);
	res->fl_res = eft_arr[4].fl_res + I*eft_arr[5].fl_res;
	res->fl_err1 = eft_arr[0].fl_err + I*eft_arr[2].fl_err;
	res->fl_err2 = -eft_arr[1].fl_err + I*eft_arr[3].fl_err;
	res->fl_err3 = eft_arr[4].fl_err + I*eft_arr[5].fl_err;
}
/* Error Free Array Extraction */
double extract(double* p,double sigma)
{
	struct eft res;
	for(int i=0; i<4; i++)
	{
		two_sum(sigma,p[i],&res);
		sigma = res.fl_res;
		p[i] = res.fl_err;
	}
	return sigma;
}
/* Sum */
double sum(const double* p)
{
	double s = p[0];
	for(int i=1; i<4; i++)
	{
		s += p[i];
	}
	return s;
}
/* Absolute Sum */
double abs_sum(const double* p)
{
	double s = fabs(p[0]);
	for(int i=1; i<4; i++)
	{
		s += fabs(p[i]);
	}
	return s;
}
/* Fast Accurate Summation */
double fast_acc_sum(double* p)
{
	// variables
	double phi, sigma, sigma_new, T, t, t_new, tau, u;
	// goto label
	start:
	T = abs_sum(p)/fma(-4,EPS,1);
	if(T<=ETA/EPS)
	{
		return sum(p);	// no rounding error
	}
	t = 0;
	do {
		sigma = (2*T)/(1 - 13*EPS);
		sigma_new = extract(p,sigma);
		tau = sigma_new - sigma;
		t_new = t;
		t = t_new + tau;
		if(t==0)
		{
			goto start;		// intermediate sum is zero, recursively apply fast_acc_sum to array of lower order parts
		}
		u = ufp(sigma);
		phi = ((48*EPS)*u)/fma(-5,EPS,1);
		T = fmin((fma(4,EPS,1.5)*(4*EPS))*sigma,(8*EPS)*u);
	} while((fabs(t) < phi) && (4*T > ETA/EPS));
	tau = (t_new - t) + tau;
	// return
	return t + (tau + sum(p));
}
/* Fast Complex Accurate Summation */
double complex fast_cmplx_acc_sum(double complex* p)
{
	// variables
	double realp[4] = {creal(p[0]),creal(p[1]),creal(p[2]),creal(p[3])};
	double imagp[4] = {cimag(p[0]),cimag(p[1]),cimag(p[2]),cimag(p[3])};
	// return
	return fast_acc_sum(realp)+I*fast_acc_sum(imagp);
}
/* Horner Method with Double Real Arithmetic */
void horner_dble(const double* poly,const double x,const unsigned int deg,double* h)
{
	// Horner's method
	*h = poly[deg];
	for(int i=deg-1; i>=0; i--)
	{
		*h = fma(*h,x,poly[i]);
	}
}
/* Reversal Horner Method with Double Real Arithmetic */
void rhorner_dble(const double* poly,const double x,const unsigned int deg,double* h)
{
	// Reversal Horner's method
	*h = poly[0];
	for(int i=1; i<=deg; i++)
	{
		*h = fma(*h,x,poly[i]);
	}
}
/* Horner Method with Complex Arithmetic */
void horner_cmplx(const double complex *poly,const double complex x,const unsigned int deg,double complex* h,double complex* hd)
{
	// Horner's method
	*h = poly[deg]; *hd = 0;
	for(int i=deg-1; i>=0; i--)
	{
		*hd = *hd*x + *h;
		*h = *h*x + poly[i];
	}
}
/* Reversal Horner Method with Complex Arithmetic */
void rhorner_cmplx(const double complex *poly,const double complex x,const unsigned int deg,double complex* h,double complex* hd)
{
	// Reversal Horner's method
	*h = poly[0]; *hd = 0;
	for(int i=1; i<=deg; i++)
	{
		*hd = *hd*x + *h;
		*h = *h*x + poly[i];
	}
}
/* Horner's Method with Complex Compensated Arithmetic */
void horner_comp_cmplx(const double complex* poly,const double complex x,const unsigned int deg,double complex* h,double complex* hd,double* eb)
{
	// local variables
	struct eft eft_arr[6];
	struct eft_cmplx_sum tsc;
	struct eft_cmplx_prod tpc;
	double complex e = 0, ed = 0;
	double complex p[4];
	double ap[4];
	// Horner's method
	*h = poly[deg]; *hd = 0; *eb = 0;
	for(int i=deg-1; i>=0; i--)
	{
		// product and sum for derivative evaluation
		two_prod_cmplx(*hd,x,eft_arr,&tpc);
		two_sum_cmplx(tpc.fl_res,*h,eft_arr,&tsc);
		// update hd and ed
		*hd = tsc.fl_res;
		//ed = ed*x + e + (tpc.fl_err1 + tpc.fl_err2 + tpc.fl_err3 + tsc.fl_err);
		p[0] = tpc.fl_err1; p[1] = tpc.fl_err2; p[2] = tpc.fl_err3; p[3] = tsc.fl_err; 
		ed = ed*x + e + fast_cmplx_acc_sum(p);
		// product and sum for polynomial evaluation
		two_prod_cmplx(*h,x,eft_arr,&tpc);
		two_sum_cmplx(tpc.fl_res,poly[i],eft_arr,&tsc);
		// update h and e
		*h = tsc.fl_res;
		//e = e*x + (tpc.fl_err1 + tpc.fl_err2 + tpc.fl_err3 + tsc.fl_err);
		p[0] = tpc.fl_err1; p[1] = tpc.fl_err2; p[2] = tpc.fl_err3; p[3] = tsc.fl_err;
		e = e*x + fast_cmplx_acc_sum(p);
		// update error bound
		ap[0] = cabs(tpc.fl_err1); ap[1] = cabs(tpc.fl_err2); ap[2] = cabs(tpc.fl_err3); ap[3] = cabs(tsc.fl_err);
		*eb = *eb*cabs(x) + fast_acc_sum(ap);
	}
	// add error back into result
	*h += e;
	*hd += ed;
}
/* Reversal Horner's Method with Complex Compensated Arithmetic */
void rhorner_comp_cmplx(const double complex* poly,const double complex x,const unsigned int deg,double complex* h,double complex* hd,double* eb)
{
	// local variables
	struct eft eft_arr[6];
	struct eft_cmplx_sum tsc;
	struct eft_cmplx_prod tpc;
	double complex e = 0, ed = 0;
	double complex p[4];
	double ap[4];
	// Horner's method
	*h = poly[0]; *hd = 0; *eb = 0;
	for(int i=1; i<=deg; i++)
	{
		// product and sum for derivative evaluation
		two_prod_cmplx(*hd,x,eft_arr,&tpc);
		two_sum_cmplx(tpc.fl_res,*h,eft_arr,&tsc);
		// update hd and ed
		*hd = tsc.fl_res;
		//ed = ed*x + e + (tpc.fl_err1 + tpc.fl_err2 + tpc.fl_err3 + tsc.fl_err);
		p[0] = tpc.fl_err1; p[1] = tpc.fl_err2; p[2] = tpc.fl_err3; p[3] = tsc.fl_err; 
		ed = ed*x + e + fast_cmplx_acc_sum(p);
		// product and sum for polynomial evaluation
		two_prod_cmplx(*h,x,eft_arr,&tpc);
		two_sum_cmplx(tpc.fl_res,poly[i],eft_arr,&tsc);
		// update h and e
		*h = tsc.fl_res;
		//e = e*x + (tpc.fl_err1 + tpc.fl_err2 + tpc.fl_err3 + tsc.fl_err);
		p[0] = tpc.fl_err1; p[1] = tpc.fl_err2; p[2] = tpc.fl_err3; p[3] = tsc.fl_err;
		e = e*x + fast_cmplx_acc_sum(p);
		// update error bound
		ap[0] = cabs(tpc.fl_err1); ap[1] = cabs(tpc.fl_err2); ap[2] = cabs(tpc.fl_err3); ap[3] = cabs(tsc.fl_err);
		*eb = *eb*cabs(x) + fast_acc_sum(ap);
	}
	// add error back into result
	*h += e;
	*hd += ed;
}
/* Horner's Method with Real Quadruple Precision */
void horner_quad_real(const mpfr_t* poly_quad,const mpfr_t x_quad,const unsigned int deg,mpfr_t* h_quad)
{
	/* initialize and set mpfr variables */
	mpfr_t mul;
	mpfr_prec_t prec = mpfr_get_prec(*h_quad);
	mpfr_init2(mul,prec);
	mpfr_set(*h_quad,poly_quad[deg],MPFR_RNDN);
	/* evaluate polynomial in quadruple precision */
	for(int i=deg-1; i>=0; i--)
	{
		mpfr_mul(mul,*h_quad,x_quad,MPFR_RNDN);
		mpfr_add(*h_quad,mul,poly_quad[i],MPFR_RNDN);
	}		
}
/* Reversal Horner's Method with Real Quadruple Precision */
void rhorner_quad_real(const mpfr_t* poly_quad,const mpfr_t x_quad,const unsigned int deg,mpfr_t* h_quad)
{
	/* initialize and set mpfr variables */
	mpfr_t mul;
	mpfr_prec_t prec = mpfr_get_prec(*h_quad);
	mpfr_init2(mul,prec);
	mpfr_set(*h_quad,poly_quad[0],MPFR_RNDN);
	/* evaluate polynomial in quadruple precision */
	for(int i=1; i<=deg; i++)
	{
		mpfr_mul(mul,*h_quad,x_quad,MPFR_RNDN);
		mpfr_add(*h_quad,mul,poly_quad[i],MPFR_RNDN);
	}		
} 
/* Horner's Method with Complex Quadruple Precision */
void horner_quad_cmplx(const mpc_t* poly_quad,const mpc_t x_quad,const unsigned int deg,mpc_t* h_quad,mpc_t* hd_quad)
{
	/* initialize and set mpc variables */
	mpc_t mul;
	mpfr_prec_t prec = mpc_get_prec(*h_quad);
	mpc_init2(mul,prec);
	mpc_set(*h_quad,poly_quad[deg],MPC_RNDNN); mpc_set_dc(*hd_quad,0,MPC_RNDNN);
	/* evaluate polynomial in quadruple precision */
	for(int i=deg-1; i>=0; i--)
	{
		// product and sum for derivative evaluation
		mpc_mul(mul,*hd_quad,x_quad,MPC_RNDNN);
		mpc_add(*hd_quad,mul,*h_quad,MPC_RNDNN);
		// product and sum for polynomial evaluation
		mpc_mul(mul,*h_quad,x_quad,MPC_RNDNN);
		mpc_add(*h_quad,mul,poly_quad[i],MPC_RNDNN);
	}
}
/* Reversal Horner's Method with Complex Quadruple Precision */
void rhorner_quad_cmplx(const mpc_t* poly_quad,const mpc_t x_quad,const unsigned int deg,mpc_t* h_quad,mpc_t* hd_quad)
{
	/* initialize and set mpc variables */
	mpc_t mul;
	mpfr_prec_t prec = mpc_get_prec(*h_quad);
	mpc_init2(mul,prec);
	mpc_set(*h_quad,poly_quad[0],MPC_RNDNN); mpc_set_dc(*hd_quad,0,MPC_RNDNN);
	/* evaluate polynomial in quadruple precision */
	for(int i=1; i<=deg; i++)
	{
		// product and sum for derivative evaluation
		mpc_mul(mul,*hd_quad,x_quad,MPC_RNDNN);
		mpc_add(*hd_quad,mul,*h_quad,MPC_RNDNN);
		// product and sum for polynomial evaluation
		mpc_mul(mul,*h_quad,x_quad,MPC_RNDNN);
		mpc_add(*h_quad,mul,poly_quad[i],MPC_RNDNN);
	}
}
#endif