#ifndef EHRLICH_ABERTH
#define EHRLICH_ABERTH
#include "init_est.h"
#include "horner.h"
#include <stdbool.h>
#include <stdio.h>
/* point convergence structure */
typedef struct
{
	bool x, y;
}	point_conv;
/* ehrlich aberth correction term */
double complex correction(const double complex* roots,const double complex h,const double complex hd,const unsigned int deg,const unsigned int j)
{
	double complex corr = 0;
	for(int i=0; i<j; i++)
	{
		corr += 1/(roots[j] - roots[i]);
	}
	for(int i=j+1; i<deg; i++)
	{
		corr += 1/(roots[j] - roots[i]);
	}
	return h/(hd - h*corr);
}
/* ehrlich aberth correction term in quadruple precision */
void correction_quad(const mpc_t* roots_quad,const mpc_t h_quad,const mpc_t hd_quad,const unsigned int deg,const unsigned int j,mpc_t* corr_quad)
{
	mpc_t temp;
	mpc_init2(temp,113);
	mpc_set_d(*corr_quad,0,MPC_RNDNN);
	for(int i=0; i<j; i++)
	{
		mpc_sub(temp,roots_quad[j],roots_quad[i],MPC_RNDNN);
		mpc_ui_div(temp,1,temp,MPC_RNDNN);
		mpc_add(*corr_quad,*corr_quad,temp,MPC_RNDNN);
	}
	for(int i=j+1; i<deg; i++)
	{
		mpc_sub(temp,roots_quad[j],roots_quad[i],MPC_RNDNN);
		mpc_ui_div(temp,1,temp,MPC_RNDNN);
		mpc_add(*corr_quad,*corr_quad,temp,MPC_RNDNN);
	}
	mpc_mul(*corr_quad,h_quad,*corr_quad,MPC_RNDNN);
	mpc_sub(*corr_quad,hd_quad,*corr_quad,MPC_RNDNN);
	mpc_div(*corr_quad,h_quad,*corr_quad,MPC_RNDNN);
}
/* reverse ehrlich aberth correction term */
double complex rcorrection(const double complex* roots,const double complex h,const double complex hd,const unsigned int deg,const unsigned int j)
{
	double complex corr = 0;
	for(int i=0; i<j; i++)
	{
		corr += 1/(roots[j] - roots[i]);
	}
	for(int i=j+1; i<deg; i++)
	{
		corr += 1/(roots[j] - roots[i]);
	}
	return cpow(roots[j],2)*h/(deg*roots[j]*h - hd - cpow(roots[j],2)*h*corr);
}
/* reverse ehrlich aberth correction term in quadruple precision */
void rcorrection_quad(const mpc_t* roots_quad,const mpc_t h_quad,const mpc_t hd_quad,const unsigned int deg,const unsigned int j,mpc_t* corr_quad)
{
	mpc_t temp;
	mpc_init2(temp,113);
	mpc_set_d(*corr_quad,0,MPC_RNDNN);
	for(int i=0; i<j; i++)
	{
		mpc_sub(temp,roots_quad[j],roots_quad[i],MPC_RNDNN);
		mpc_ui_div(temp,1,temp,MPC_RNDNN);
		mpc_add(*corr_quad,*corr_quad,temp,MPC_RNDNN);
	}
	for(int i=j+1; i<deg; i++)
	{
		mpc_sub(temp,roots_quad[j],roots_quad[i],MPC_RNDNN);
		mpc_ui_div(temp,1,temp,MPC_RNDNN);
		mpc_add(*corr_quad,*corr_quad,temp,MPC_RNDNN);
	}
	mpc_mul(*corr_quad,h_quad,*corr_quad,MPC_RNDNN);
	mpc_sqr(temp,roots_quad[j],MPC_RNDNN);
	mpc_mul(*corr_quad,temp,*corr_quad,MPC_RNDNN);
	mpc_add(*corr_quad,hd_quad,*corr_quad,MPC_RNDNN);
	mpc_mul(temp,roots_quad[j],h_quad,MPC_RNDNN);
	mpc_mul_ui(temp,temp,deg,MPC_RNDNN);
	mpc_sub(*corr_quad,temp,*corr_quad,MPC_RNDNN);
	mpc_sqr(temp,roots_quad[j],MPC_RNDNN);
	mpc_mul(temp,temp,h_quad,MPC_RNDNN);
	mpc_div(*corr_quad,temp,*corr_quad,MPC_RNDNN);
}
/* ehrlich_aberth*/
void ehrlich_aberth(const double complex* poly,double complex* roots,const unsigned int deg,const unsigned int itmax)
{
	// local variables
	bool s;
	double b;
	double complex h,hd;
	// local arrays
	double alpha[deg+1];
	bool conv[deg];
	// initial estimates
	for(int i=0; i<deg; i++)
	{
		alpha[i] = cabs(poly[i]);
		conv[i] = false;
	}
	alpha[deg] = cabs(poly[deg]);
	init_est(alpha,deg,roots);
	// update initial estimates
	for(int i=0; i<=deg; i++)
	{
		alpha[i] = alpha[i]*fma(3.8,i,1);
	}
	for(int i=0; i<itmax; i++)
	{
		for(int j=0; j<deg; j++)
		{
			if(conv[j]==0)
			{
				if(cabs(roots[j])>1)
				{
					rhorner_dble(alpha,1/cabs(roots[j]),deg,&b);
					rhorner_cmplx(poly,1/roots[j],deg,&h,&hd);
					if(cabs(h)>EPS*b)
					{
						roots[j] = roots[j] - rcorrection(roots,h,hd,deg,j);
					}
					else
					{
						conv[j] = true;
					}
					
				}
				else
				{
					horner_dble(alpha,cabs(roots[j]),deg,&b);
					horner_cmplx(poly,roots[j],deg,&h,&hd);	
					if(cabs(h)>EPS*b)
					{
						roots[j] = roots[j] - correction(roots,h,hd,deg,j);	
					}
					else
					{
						conv[j] = true;
					}
				}
			}
		}
		s = conv[0];
		for(int j=1; j<deg; j++)
		{
			s = s && conv[j];
		}
		if(s)
		{
			break;
		}
	}
	if(!s)
	{
		printf("not all roots converged\n");
	}
}
/* ehrlich_aberth with compensated arithmetic*/
void ehrlich_aberth_comp(const double complex* poly,double complex* roots,const unsigned int deg,const unsigned int itmax)
{
	// local variables
	int s;
	double b;
	double complex h,hd;
	// local arrays
	double alpha[deg+1];
	point_conv conv[deg];
	// initial estimates
	for(int i=0; i<deg; i++)
	{
		alpha[i] = cabs(poly[i]);
		conv[i].x = false;
		conv[i].y = false;
	}
	alpha[deg] = cabs(poly[deg]);
	init_est(alpha,deg,roots);
	// update initial estimates
	for(int i=0; i<=deg; i++)
	{
		alpha[i] = alpha[i]*fma(3.8,i,1);
	}
	for(int i=0; i<itmax; i++)
	{
		for(int j=0; j<deg; j++)
		{
			if(!conv[j].x)
			{
				if(cabs(roots[j])>1)
				{
					rhorner_dble(alpha,1/cabs(roots[j]),deg,&b);
					rhorner_cmplx(poly,1/roots[j],deg,&h,&hd);
					if(cabs(h)>EPS*b)
					{
						roots[j] = roots[j] - rcorrection(roots,h,hd,deg,j);
					}
					else
					{
						conv[j].x = true;
					}
					
				}
				else
				{
					horner_dble(alpha,cabs(roots[j]),deg,&b);
					horner_cmplx(poly,roots[j],deg,&h,&hd);	
					if(cabs(h)>EPS*b)
					{
						roots[j] = roots[j] - correction(roots,h,hd,deg,j);	
					}
					else
					{
						conv[j].x = true;	
					}
				}
			}
			else if(!conv[j].y)
			{
				if(cabs(roots[j])>1)
				{
					rhorner_comp_cmplx(poly,1/roots[j],deg,&h,&hd,&b);
					double errBound = EPS*cabs(h) + (gamma_const(4*deg+2)*b + 2*pow(EPS,2)*cabs(h));
					if(cabs(h) > 4*errBound)
					{
						double complex corr = rcorrection(roots,h,hd,deg,j);
						if(cabs(corr) > 4*EPS*cabs(roots[j]))
						{
							roots[j] = roots[j] - corr;
						}
						else
						{
							conv[j].y = true;
						}
					}
					else
					{
						conv[j].y = true;
					}
					
				}
				else
				{
					horner_comp_cmplx(poly,roots[j],deg,&h,&hd,&b);
					double errBound = EPS*cabs(h) + (gamma_const(4*deg+2)*b + 2*pow(EPS,2)*cabs(h));
					if(cabs(h) > 4*errBound)
					{
						double complex corr = correction(roots,h,hd,deg,j);
						if(cabs(corr) > 4*EPS)
						{
							roots[j] = roots[j] - corr;	
						}
						else
						{
							conv[j].y = true;
						}
					}
					else
					{
						conv[j].y = true;
					}
				}	
			}
		}
		s = conv[0].y;
		for(int j=1; j<deg; j++)
		{
			s = s && conv[j].y;
		}
		if(s)
		{
			break;
		}
	}
	if(!s)
	{
		printf("not all roots comp converged\n");
	}
}
/* ehrlich_aberth in quadruple precision */
void ehrlich_aberth_quad(const mpc_t* poly_quad,mpc_t* roots_quad,const unsigned int deg,const unsigned int itmax)
{
	// local variables
	bool s;
	mpfr_t b_quad, eps_quad, x_quad;
	mpc_t corr_quad, h_quad, hd_quad;
	mpfr_init2(b_quad,113); mpfr_init2(eps_quad,113); mpfr_init2(x_quad,113);
	mpc_init2(corr_quad,113); mpc_init2(h_quad,113); mpc_init2(hd_quad,113);
	mpfr_set_d(eps_quad,pow(2,-113),MPFR_RNDN);
	// local arrays
	mpfr_t alpha_quad[deg+1];
	bool conv[deg];
	// initial estimates
	for(int i=0; i<deg; i++)
	{
		mpfr_init2(alpha_quad[i],113);
		mpc_abs(alpha_quad[i],poly_quad[i],MPFR_RNDN);
		conv[i] = false;
	}
	mpfr_init2(alpha_quad[deg],113);
	mpc_abs(alpha_quad[deg],poly_quad[deg],MPFR_RNDN);
	init_est_quad(alpha_quad,deg,roots_quad);
	// update initial estiamtes
	for(int i=0; i<=deg; i++)
	{
		mpfr_mul_d(alpha_quad[i],alpha_quad[i],fma(3.8,i,1),MPFR_RNDN);
	}
	for(int i=0; i<itmax; i++)
	{
		for(int j=0; j<deg; j++)
		{
			if(!conv[j])
			{
				mpc_abs(x_quad,roots_quad[j],MPFR_RNDN);
				if(mpfr_cmp_ui(x_quad,1)>0)
				{
					mpfr_ui_div(x_quad,1,x_quad,MPFR_RNDN);
					rhorner_mp_real(alpha_quad,x_quad,deg,&b_quad);
					mpc_ui_div(corr_quad,1,roots_quad[j],MPC_RNDNN);
					rhorner_mp_cmplx(poly_quad,corr_quad,deg,&h_quad,&hd_quad);
					mpfr_mul(b_quad,b_quad,eps_quad,MPFR_RNDN);
					mpc_abs(x_quad,h_quad,MPFR_RNDN);
					if(mpfr_cmp(x_quad,b_quad)>0)
					{
						rcorrection_quad(roots_quad,h_quad,hd_quad,deg,j,&corr_quad);
						mpc_sub(roots_quad[j],roots_quad[j],corr_quad,MPC_RNDNN);	
					}
					else
					{
						conv[j] = true;
					}
				}
				else
				{
					horner_mp_real(alpha_quad,x_quad,deg,&b_quad);
					horner_mp_cmplx(poly_quad,roots_quad[j],deg,&h_quad,&hd_quad);
					mpfr_mul(b_quad,b_quad,eps_quad,MPFR_RNDN);
					mpc_abs(x_quad,h_quad,MPFR_RNDN);
					if(mpfr_cmp(x_quad,b_quad)>0)
					{
						correction_quad(roots_quad,h_quad,hd_quad,deg,j,&corr_quad);
						mpc_sub(roots_quad[j],roots_quad[j],corr_quad,MPC_RNDNN);	
					}
					else
					{
						conv[j] = true;
					}
				}
			}
		}
		s = conv[0];
		for(int j=1; j<deg; j++)
		{
			s = s && conv[j];
		}
		if(s)
		{
			break;
		}
	}
	if(!s)
	{
		printf("not all roots quad converged\n");
	}
	// clear mpfr variables
	mpfr_clears(b_quad, eps_quad, x_quad, (mpfr_ptr) 0);
	for(int i=0; i<=deg; i++)
	{
		mpfr_clear(alpha_quad[i]);
	}
	// clear mpc variables
	mpc_clear(corr_quad); mpc_clear(h_quad); mpc_clear(hd_quad);
}
#endif