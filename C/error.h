#ifndef ERROR
#define ERROR
#include "horner.h"
/* Global Constants */
const unsigned int PREC = 131072;
/* maximum value */
double max_value(const double* values,const unsigned int n)
{
	double max = values[0];
	for(int i=1; i<n; i++)
	{
		if(max < values[i])
		{
			max = values[i];
		}
	}
	return max;
}
/* roots to polynomial */
void roots_to_poly(double complex* poly,const double complex* roots,const double complex lead_coeff,const unsigned int deg)
{
	/* initialize multi-precision variables */
	mpfr_t err_mp, abs_mp;
	mpfr_init2(err_mp,PREC); mpfr_init2(abs_mp,PREC);
	mpc_t poly_mp[deg+1], roots_mp[deg], mul_mp, sub_mp;
	mpc_init2(mul_mp,PREC); mpc_init2(sub_mp,PREC);
	for(int i=0; i<deg; i++)
	{
		mpc_init2(poly_mp[i],PREC);
		mpc_set_dc(poly_mp[i],0,MPC_RNDNN);
		mpc_init2(roots_mp[i],PREC);
		mpc_set_dc(roots_mp[i],roots[i],MPC_RNDNN);		
	}
	mpc_init2(poly_mp[deg],PREC);
	mpc_set_dc(poly_mp[deg],lead_coeff,MPC_RNDNN);
	/* expand polynomial */
	for(int i=0; i<deg; i++)
	{
		for(int j=i; j>=0; j--)
		{
			if(j==0)
			{
				mpc_sub(poly_mp[j],poly_mp[j],roots_mp[i],MPC_RNDNN);
			}
			else
			{
				mpc_mul(mul_mp,roots_mp[i],poly_mp[j-1],MPC_RNDNN);
				mpc_sub(poly_mp[j],poly_mp[j],mul_mp,MPC_RNDNN);	
			}
		}
	}
	/* multiply by leading coefficient */
	for(int i=0; i<deg; i++)
	{
		mpc_mul(poly_mp[i],poly_mp[i],poly_mp[deg],MPC_RNDNN);
	}
	/* round to double precision */
	for(int i=0; i<deg; i++)
	{
		poly[i] = mpc_get_dc(poly_mp[deg-i-1],MPC_RNDNN);
	}
	poly[deg] = lead_coeff;
}
/* backward error */
double back_err(const double complex* poly,const double complex* roots,const unsigned int deg)
{
	/* initialize multi-precision variables */
	mpfr_t err_mp, abs_mp, max_mp;
	mpfr_init2(err_mp,PREC); mpfr_init2(abs_mp,PREC); mpfr_init2(max_mp,PREC);
	mpc_t poly_mp[deg+1], roots_mp[deg], mul_mp, sub_mp;
	mpc_init2(mul_mp,PREC); mpc_init2(sub_mp,PREC);
	for(int i=0; i<deg; i++)
	{
		mpc_init2(poly_mp[i],PREC);
		mpc_set_dc(poly_mp[i],0,MPC_RNDNN);
		mpc_init2(roots_mp[i],PREC);
		mpc_set_dc(roots_mp[i],roots[i],MPC_RNDNN);		
	}
	mpc_init2(poly_mp[deg],PREC);
	mpc_set_dc(poly_mp[deg],poly[deg],MPC_RNDNN);
	/* expand polynomial */
	for(int i=0; i<deg; i++)
	{
		for(int j=i; j>=0; j--)
		{
			if(j==0)
			{
				mpc_sub(poly_mp[j],poly_mp[j],roots_mp[i],MPC_RNDNN);
			}
			else
			{
				mpc_mul(mul_mp,roots_mp[i],poly_mp[j-1],MPC_RNDNN);
				mpc_sub(poly_mp[j],poly_mp[j],mul_mp,MPC_RNDNN);	
			}
		}
	}
	/* multiply by leading coefficient */
	for(int i=0; i<deg; i++)
	{
		mpc_mul(poly_mp[i],poly_mp[i],poly_mp[deg],MPC_RNDNN);
	}
	/* compute backward error */
	mpfr_set_d(err_mp,0,MPC_RNDNN);
	for(int i=0; i<deg; i++)
	{
		mpc_set_dc(mul_mp,poly[i],MPC_RNDNN);
		mpc_sub(sub_mp,poly_mp[deg - i - 1],mul_mp,MPC_RNDNN);
		mpc_abs(abs_mp,sub_mp,MPFR_RNDN);
		if(mpfr_cmp(abs_mp,err_mp)>0)
		{
			mpfr_set(err_mp,abs_mp,MPFR_RNDN);
		}
	}
	mpc_abs(max_mp,poly_mp[0],MPFR_RNDN);
	for(int i=1; i<=deg; i++)
	{
		mpc_abs(abs_mp,poly_mp[i],MPFR_RNDN);
		if(mpfr_cmp(abs_mp,max_mp)>0)
		{
			mpfr_set(max_mp,abs_mp,MPFR_RNDN);
		}
	}
	mpfr_div(err_mp,err_mp,max_mp,MPFR_RNDN);
	// return largest value between err_mp and EPS
	double p[2] = {mpfr_get_d(err_mp,MPFR_RNDN),EPS};
	return max_value(p,2);
}
/* forward error */
double forw_err(const double complex* roots,const double complex* exact_roots,const unsigned int deg)
{
	// spectral variation 1
	double sv1 = 0;
	bool* used = (bool*)malloc(deg*sizeof(bool));
	for(int i=0; i<deg; i++)
	{
		used[i] = false;
	}
	for(int i=0; i<deg; i++)
	{
		double min = DBL_MAX;
		int ind;
		for(int j=0; j<deg; j++)
		{
			if(!used[j])
			{
				double temp = cabs(roots[i]-exact_roots[j]);
				if(temp < min)
				{
					min = temp;
					ind = j;
				}
			}
		}
		used[ind] = true;
		if(min > sv1)
		{
			sv1 = min;
		}
	}
	// spectral variation 2
	double sv2 = 0;
	for(int i=0; i<deg; i++)
	{
		used[i] = false;
	}
	for(int i=0; i<deg; i++)
	{
		double min = DBL_MAX;
		int ind;
		for(int j=0; j<deg; j++)
		{
			if(!used[j])
			{
				double temp = cabs(exact_roots[i] - roots[j]);
				if(temp < min)
				{
					min = temp;
					ind = j;
				}
			}
		}
		used[ind] = true;
		if(min > sv2)
		{
			sv2 = min;
		}
	}
	// largest exact root
	double max = cabs(exact_roots[0]);
	for(int i=1; i<deg; i++)
	{
		double temp = cabs(exact_roots[i]);
		if(temp > max)
		{
			max = temp;
		}
	}
	// return largest value of sv1/max, sv2/max, and EPS
	double p[3] = {sv1/max,sv2/max,EPS};
	return max_value(p,3);
}
/* maximum condition number */
double max_cond(const double complex* poly,const double complex* roots,const unsigned int deg)
{
	/* initialize multi-precision variables */
	mpfr_t alpha_mp[deg+1], b_mp, x_mp, y_mp;
	mpfr_init2(b_mp,PREC); mpfr_init2(x_mp,PREC); mpfr_init2(y_mp,PREC);
	mpc_t poly_mp[deg+1], roots_mp[deg], h_mp, hd_mp;
	mpc_init2(h_mp,PREC); mpc_init2(hd_mp,PREC);
	for(int i=0; i<deg; i++)
	{
		mpc_init2(poly_mp[i],PREC);
		mpc_set_dc(poly_mp[i],poly[i],MPC_RNDNN);
		mpfr_init2(alpha_mp[i],PREC);
		mpc_abs(alpha_mp[i],poly_mp[i],MPFR_RNDN);
		mpc_init2(roots_mp[i],PREC);
		mpc_set_dc(roots_mp[i],roots[i],MPC_RNDNN);		
	}
	mpc_init2(poly_mp[deg],PREC);
	mpc_set_dc(poly_mp[deg],poly[deg],MPC_RNDNN);
	mpfr_init2(alpha_mp[deg],PREC);
	mpc_abs(alpha_mp[deg],poly_mp[deg],MPFR_RNDN);
	/* compute maximum condition number */
	double max = 0;
	for(int i=0; i<deg; i++)
	{
		mpc_abs(x_mp,roots_mp[i],MPFR_RNDN);
		horner_quad_real(alpha_mp,x_mp,deg,&b_mp);
		horner_quad_cmplx(poly_mp,roots_mp[i],deg,&h_mp,&hd_mp);
		mpc_abs(y_mp,hd_mp,MPFR_RNDN);
		mpfr_mul(x_mp,x_mp,y_mp,MPFR_RNDN);
		mpfr_div(x_mp,b_mp,x_mp,MPFR_RNDN);
		double temp = mpfr_get_d(x_mp,MPFR_RNDN);
		if(temp > max)
		{
			max = temp;
		}
	}			
	// return maximum condition number
	return max;
}
#endif