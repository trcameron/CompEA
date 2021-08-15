/* Horner Test */
#include "horner.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/* Main Function */
int main(int argc,char **argv)
{
	/* horner error bound test */
	int deg = 5;
	/* initialize polynomials */
	double alpha[deg+1];
	double complex poly[deg+1];
	mpfr_t alpha_quad[deg+1];
	mpc_t poly_quad[deg+1];
	mpz_t mp_bin;
	mpz_init(mp_bin);
	for(int k=0; k<=deg; k++)
	{
		mpz_bin_uiui(mp_bin,deg,k);
		poly[k] = cpow(-1-I,deg-k)*mpz_get_ui(mp_bin);
		alpha[k] = cabs(poly[k]);
		mpc_init2(poly_quad[k],113);
		mpc_set_dc(poly_quad[k],poly[k],MPC_RNDNN);
		mpfr_init2(alpha_quad[k],113);
		mpfr_set_d(alpha_quad[k],alpha[k],MPFR_RNDN);
	}
	/* initialize writing file */
	FILE* f;
	f = fopen("data_files/horner_errbound_test.dat","w+");
	fprintf(f,"x, a priori error bound, running error bound, forward error\n");
	/* initialize multi-precision variables */
	mpfr_t eb_quad, x_quad;
	mpfr_init2(eb_quad,113); mpfr_init2(x_quad,113);
	mpc_t err_quad, h_quad, hd_quad, xc_quad;
	mpc_init2(err_quad,113); mpc_init2(h_quad,113); mpc_init2(hd_quad,113); mpc_init2(xc_quad,113);
	/* initialize horner_comp variables */
	double eb;
	double complex h, hd;
	/* initialize testing variables */
	const unsigned int num = 2E+3;
	const double dx = 1E-5;
	double comp_et = 0, quad_et = 0;
	double complex x = 0.99 + I;
	clock_t begin, end;
	/* run test */
	for(int i=0; i<=num; i++)
	{
		// write creal(x) value
		fprintf(f,"%.5e, ",creal(x));
		// quadruple precision
		mpc_set_dc(xc_quad,x,MPC_RNDNN);
		begin = clock();
		horner_mp_cmplx(poly_quad,xc_quad,deg,&h_quad,&hd_quad);
		end = clock();
		quad_et += (double)(end - begin) / CLOCKS_PER_SEC;
		// a priori error bound
		mpfr_set_d(x_quad,cabs(x),MPFR_RNDN);
		horner_mp_real(alpha_quad,x_quad,deg,&eb_quad);
		mpc_abs(x_quad,h_quad,MPFR_RNDN);
		fprintf(f,"%.5e, ",EPS*mpfr_get_d(x_quad,MPFR_RNDN)+pow(gamma_const(2*deg),2)*mpfr_get_d(eb_quad,MPFR_RNDN));
		// compensated arithmetic
		begin = clock();
		horner_comp_cmplx(poly,x,deg,&h,&hd,&eb);
		end = clock();
		comp_et += (double)(end - begin) / CLOCKS_PER_SEC;
		// running error bound
		fprintf(f,"%.5e, ",EPS*cabs(h) + (gamma_const(4*deg+2)*eb + 2*pow(EPS,2)*cabs(h)));
		// forward error bound
		mpc_set_dc(hd_quad,-h,MPC_RNDNN);
		mpc_add(err_quad,h_quad,hd_quad,MPC_RNDNN);
		fprintf(f,"%.5e\n",cabs(mpc_get_dc(err_quad,MPC_RNDNN)));
		// update x
		x += dx;
	}
	/* close file */
	fclose(f);
	/* clear mpz variables */
	mpz_clear(mp_bin);
	/* clear mpfr variables */
	mpfr_clears(eb_quad,x_quad,(mpfr_ptr) 0);
	for(int i=0; i<=deg; i++)
	{
		mpfr_clear(alpha_quad[i]);
	}
	/* clear mpc variables */
	mpc_clear(err_quad);
	for(int i=0; i<=deg; i++)
	{
		mpc_clear(poly_quad[i]);
	}
	/* print elapsed times */
	printf("comp horner elapsed time = %.4e\n",comp_et);
	printf("quad horner elapsed time = %.4e\n",quad_et);
	printf("quad to comp ratio = %.4e\n",quad_et/comp_et);
	
	/* horner time test */
	f = fopen("data_files/horner_time_test.dat","w+");
	fprintf(f,"deg, double_time, comp_time, quad_time\n");
	/* initialize testing variables */
	const unsigned int min_deg = 100, max_deg = 2500, itnum = 100;
	for(deg = min_deg; deg <= max_deg; deg+=10)
	{
		// write deg value
		fprintf(f,"%d, ",deg);
		// initialize random polynomial
		double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
		// initialize quad-precision polynomial
		mpc_t* poly_quad = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
		// seed random variable
		srand(time(NULL));
		// initialize timing variables
		double dble_et = 0;
		comp_et = 0; quad_et = 0;
		for(int i=0; i<itnum; i++)
		{
			// store random polynomial
			for(int i=0; i<=deg; i++)
			{
				double real = (double)2*rand()/RAND_MAX - 1;
				double imag = (double)2*rand()/RAND_MAX - 1;
				poly[i] = real + I*imag;
			}
			// store random value
			x = ((double)2*rand()/RAND_MAX - 1) + I*((double)2*rand()/RAND_MAX - 1);
			// store quad-precision polynomial
			for(int i=0; i<=deg; i++)
			{
				mpc_init2(poly_quad[i],113);
				mpc_set_dc(poly_quad[i],poly[i],MPC_RNDNN);
			}
			// double precision
			begin = clock();
			horner_cmplx(poly,x,deg,&h,&hd);
			end = clock();
			dble_et += (double)(end - begin) / CLOCKS_PER_SEC;
			// compensated arithmetic
			begin = clock();
			horner_comp_cmplx(poly,x,deg,&h,&hd,&eb);
			end = clock();
			comp_et += (double)(end - begin) / CLOCKS_PER_SEC;
			// quadruple precision
			mpc_set_dc(xc_quad,x,MPC_RNDNN);
			begin = clock();
			horner_mp_cmplx(poly_quad,xc_quad,deg,&h_quad,&hd_quad);
			end = clock();
			quad_et += (double)(end - begin) / CLOCKS_PER_SEC;
		}
		// free memory
		free(poly); free(poly_quad);
		// write avg elapsed times
		fprintf(f,"%.5e, ",dble_et/itnum);
		fprintf(f,"%.5e, ",comp_et/itnum);
		fprintf(f,"%.5e\n",quad_et/itnum);
	}
	/* close file */
	fclose(f);
	/* clear mpc variables */
	 mpc_clear(h_quad); mpc_clear(hd_quad); mpc_clear(xc_quad);
	/* return */
	return 0;
}