/* Special Polynomial Test */
#include "ehrlich_aberth.h"
#include "error.h"
#include "specpoly.h"
#include <stdio.h>
#include <stdlib.h>
/* Main Function */
int main(int argc,char **argv)
{
	/* testing variables */
	FILE* f;
	const unsigned int itmax = 50;
	unsigned int deg_min = 5, deg_max = 80;
	/* Chebyshev Polynomial */
	f = fopen("data_files/chebyshev_test.dat","w+");
	fprintf(f,"cond, ea_err, ea_bound, comp_err, comp_bound\n");
	for(int deg = deg_min; deg <= deg_max; deg++)
	{
		// initialize storage for error
		double ea_err, comp_err;
		// initialize storage for bounds
		double ea_bound, comp_bound;
		// initialize storage for condition number
		double cond;
		// initialize Chebyshev polynomial based on recursive formula and known roots
		double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
		double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
		double complex* exact_roots = (double complex*)malloc(deg*sizeof(double complex));
		// compute chebyshev polynomial coefficients
		double complex* a = (double complex*)malloc((deg+1)*sizeof(double complex));
		double complex* b = (double complex*)malloc((deg+1)*sizeof(double complex));
		for(int i=0; i<=deg; i++)
		{
			a[i] = 0; b[i] = 0;
		}
		a[0] = 1; b[1] = 1;
		for(int i=2; i<=deg; i++)
		{
			poly[0] = -a[0];
			for(int j=1; j<=i; j++)
			{
				poly[j] = 2*b[j-1] - a[j];
			}
			for(int j=0; j<=i; j++)
			{
				a[j] = b[j];
				b[j] = poly[j];
			}
		}
		free(a); free(b);
		// compute chevyshev polynomial roots
		mpfr_t pie, root;
		mpfr_init2(pie,128); mpfr_init2(root,128);
		mpfr_const_pi(pie,MPFR_RNDN);
		for(int i=0; i<deg; i++)
		{
			mpfr_mul_ui(root,pie,2*i+1,MPFR_RNDN);
			mpfr_div_ui(root,root,2*deg,MPFR_RNDN);
			mpfr_cos(root,root,MPFR_RNDN);
			exact_roots[i] = mpfr_get_d(root,MPFR_RNDN);
		}
		// free mpfr variables
		mpfr_clears(pie,root,(mpfr_ptr) 0);
		// ehrlich-aberth methods
		ehrlich_aberth(poly,roots,deg,itmax);
		ea_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_comp(poly,roots,deg,itmax);
		comp_err = forw_err(roots,exact_roots,deg);
		// limiting accuracy bounds
		cond = max_cond2(poly,exact_roots,deg);
		ea_bound = gamma_const(2*deg)*cond;
		if(ea_bound > 1)
		{
			ea_bound = 1;
		}
		comp_bound = EPS + pow(gamma_const(2*deg),2)*cond;
		if(comp_bound > 1)
		{
			comp_bound = 1;
		}
		// write to file
		fprintf(f,"%.5e, ",cond);
		fprintf(f,"%.5e, ",ea_err);
		fprintf(f,"%.5e, ",ea_bound);
		fprintf(f,"%.5e, ",comp_err);
		fprintf(f,"%.5e\n",comp_bound);
		// free memory
		free(poly); free(roots); free(exact_roots);
	}
	// close file
	fclose(f);
	/* Wilkinson Polynomial */
	deg_min = 5, deg_max = 20;
	f = fopen("data_files/wilkinson_test.dat","w+");
	fprintf(f,"degree, ea_err, ea_comp_err, ea_quad_err, max_cond\n");
	for(int deg = deg_min; deg <= deg_max; deg++)
	{
		// initialize storage for error
		double ea_err, comp_err, quad_err;
		// initialize Wilkonson polynomial based on known roots
		double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
		double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
		double complex* exact_roots = (double complex*)malloc(deg*sizeof(double complex));
		for(int i=0; i<deg; i++)
		{
			exact_roots[i] = (i+1);
		}
		roots_to_poly(poly,exact_roots,1,deg);
		// initialize quad-precision polynomial and storage for quad-precision roots
		mpc_t* poly_quad = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
		mpc_t* roots_quad = (mpc_t*)malloc(deg*sizeof(mpc_t));
		for(int i=0; i<deg; i++)
		{
			mpc_init2(poly_quad[i],113); mpc_init2(roots_quad[i],113);
			mpc_set_dc(poly_quad[i],poly[i],MPC_RNDNN);
		}
		mpc_init2(poly_quad[deg],113);
		mpc_set_dc(poly_quad[deg],poly[deg],MPC_RNDNN);
		// ehrlich-aberth methods
		ehrlich_aberth(poly,roots,deg,itmax);
		ea_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_comp(poly,roots,deg,itmax);
		comp_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
		for(int i=0; i<deg; i++)
		{
			roots[i] = mpc_get_dc(roots_quad[i],MPC_RNDNN);
		}
		quad_err = forw_err(roots,exact_roots,deg);
		// write to file
		fprintf(f,"%d, ",deg);
		fprintf(f,"%.5e, ",ea_err);
		fprintf(f,"%.5e, ",comp_err);
		fprintf(f,"%.5e, ",quad_err);
		fprintf(f,"%.5e\n",max_cond2(poly,exact_roots,deg));
		// free memory
		free(poly); free(roots); free(exact_roots); free(poly_quad); free(roots_quad);
	}
	// close file
	fclose(f);
	/* Prescribed-Roots Polynomial */
	deg_min = 5, deg_max = 20;
	f = fopen("data_files/presroots_test.dat","w+");
	fprintf(f,"degree, ea_err, ea_comp_err, ea_quad_err, max_cond\n");
	for(int deg = deg_min; deg <= deg_max; deg++)
	{
		// initialize storage for error
		double ea_err, comp_err, quad_err;
		// initialize prescribed-roots poly based on known roots
		double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
		double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
		double complex* exact_roots = (double complex*)malloc(deg*sizeof(double complex));
		for(int i=0; i<deg; i++)
		{
			exact_roots[i] = pow(2,-floor((double)deg/2)+i) - 3;
		}
		roots_to_poly(poly,exact_roots,1,deg);
		// initialize quad-precision polynomial and storage for quad-precision roots
		mpc_t* poly_quad = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
		mpc_t* roots_quad = (mpc_t*)malloc(deg*sizeof(mpc_t));
		for(int i=0; i<deg; i++)
		{
			mpc_init2(poly_quad[i],113); mpc_init2(roots_quad[i],113);
			mpc_set_dc(poly_quad[i],poly[i],MPC_RNDNN);
		}
		mpc_init2(poly_quad[deg],113);
		mpc_set_dc(poly_quad[deg],poly[deg],MPC_RNDNN);
		// ehrlich-aberth methods
		ehrlich_aberth(poly,roots,deg,itmax);
		ea_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_comp(poly,roots,deg,itmax);
		comp_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
		for(int i=0; i<deg; i++)
		{
			roots[i] = mpc_get_dc(roots_quad[i],MPC_RNDNN);
		}
		quad_err = forw_err(roots,exact_roots,deg);
		// write to file
		fprintf(f,"%d, ",deg);
		fprintf(f,"%.5e, ",ea_err);
		fprintf(f,"%.5e, ",comp_err);
		fprintf(f,"%.5e, ",quad_err);
		fprintf(f,"%.5e\n",max_cond2(poly,exact_roots,deg));
		// free memory
		free(poly); free(roots); free(exact_roots); free(poly_quad); free(roots_quad);
	}
	// close file
	fclose(f);
	/* Small-Imaginary Polynomial */
	deg_min = 5, deg_max = 25;
	f = fopen("data_files/smallimag_test.dat","w+");
	fprintf(f,"degree, ea_err, ea_comp_err, ea_quad_err, max_cond\n");
	for(int deg = deg_min; deg <= deg_max; deg++)
	{
		// initialize storage for error
		double ea_err, comp_err, quad_err;
		// initialize small-imag poly based on known roots
		double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
		double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
		double complex* exact_roots = (double complex*)malloc(deg*sizeof(double complex));
		for(int i=0; i<deg; i++)
		{
			exact_roots[i] = (i+1) + I*pow(-1,i+1)*8*EPS;
		}
		roots_to_poly(poly,exact_roots,1,deg);
		// initialize quad-precision polynomial and storage for quad-precision roots
		mpc_t* poly_quad = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
		mpc_t* roots_quad = (mpc_t*)malloc(deg*sizeof(mpc_t));
		for(int i=0; i<deg; i++)
		{
			mpc_init2(poly_quad[i],113); mpc_init2(roots_quad[i],113);
			mpc_set_dc(poly_quad[i],poly[i],MPC_RNDNN);
		}
		mpc_init2(poly_quad[deg],113);
		mpc_set_dc(poly_quad[deg],poly[deg],MPC_RNDNN);
		// ehrlich-aberth methods
		ehrlich_aberth(poly,roots,deg,itmax);
		ea_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_comp(poly,roots,deg,itmax);
		comp_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
		for(int i=0; i<deg; i++)
		{
			roots[i] = mpc_get_dc(roots_quad[i],MPC_RNDNN);
		}
		quad_err = forw_err(roots,exact_roots,deg);
		// write to file
		fprintf(f,"%d, ",deg);
		fprintf(f,"%.5e, ",ea_err);
		fprintf(f,"%.5e, ",comp_err);
		fprintf(f,"%.5e, ",quad_err);
		fprintf(f,"%.5e\n",max_cond2(poly,exact_roots,deg));
		// free memory
		free(poly); free(roots); free(exact_roots); free(poly_quad); free(roots_quad);
	}
	// close file
	fclose(f);
	/* polynomials with multiple and near multipe roots */
	f = fopen("data_files/multroots_test.dat","w+");
	fprintf(f,"poly_num, limit_acc, ea_err, ea_comp_err, ea_quad_err\n");
	for(unsigned int poly_num = 1; poly_num <= 7; poly_num++)
	{
		// initialize storage for error
		double ea_err, comp_err, quad_err;
		// initialize poly and exact_roots
		unsigned int deg;
		double complex* poly;
		double complex* exact_roots;
		spec_poly(&poly,&exact_roots,&deg,poly_num);
		// initialize roots
		double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
		// initialize quad-precision polynomial and storage for quad-precision roots
		mpc_t* poly_quad = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
		mpc_t* roots_quad = (mpc_t*)malloc(deg*sizeof(mpc_t));
		for(int i=0; i<deg; i++)
		{
			mpc_init2(poly_quad[i],113); mpc_init2(roots_quad[i],113);
			mpc_set_dc(poly_quad[i],poly[i],MPC_RNDNN);
		}
		mpc_init2(poly_quad[deg],113);
		mpc_set_dc(poly_quad[deg],poly[deg],MPC_RNDNN);
		// ehrlich-aberth methods
		ehrlich_aberth(poly,roots,deg,itmax);
		ea_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_comp(poly,roots,deg,itmax);
		comp_err = forw_err(roots,exact_roots,deg);
		ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
		for(int i=0; i<deg; i++)
		{
			roots[i] = mpc_get_dc(roots_quad[i],MPC_RNDNN);
		}
		quad_err = forw_err(roots,exact_roots,deg);
		// write to file
		fprintf(f,"%d, ",poly_num);
		fprintf(f,"%.5e, ",EPS + pow(gamma_const(2*deg),2)*max_cond2(poly,exact_roots,deg));
		fprintf(f,"%.5e, ",ea_err);
		fprintf(f,"%.5e, ",comp_err);
		fprintf(f,"%.5e\n",quad_err);
		// free memory
		free(poly); free(roots); free(exact_roots); free(poly_quad); free(roots_quad);
	}
	// close file
	fclose(f);
	// return
	return 0;
}