/* Special Polynomial Test */
#include "ehrlich_aberth.h"
#include "error.h"
#include <stdio.h>
#include <stdlib.h>
/* Main Function */
int main(int argc,char **argv)
{
	/* testing variables */
	FILE* f;
	const unsigned int itmax = 30;
	unsigned int deg_min = 5, deg_max = 20;
	/* Wilkinson Polynomial */
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
		fprintf(f,"%.5e\n",max_cond(poly,exact_roots,deg));
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
		// initialize Wilkonson polynomial based on known roots
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
		fprintf(f,"%.5e\n",max_cond(poly,exact_roots,deg));
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
		// initialize Wilkonson polynomial based on known roots
		double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
		double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
		double complex* exact_roots = (double complex*)malloc(deg*sizeof(double complex));
		for(int i=0; i<deg; i++)
		{
			exact_roots[i] = (i+1) + I*pow(-1,i+1)*10*EPS;
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
		fprintf(f,"%.5e\n",max_cond(poly,exact_roots,deg));
		// free memory
		free(poly); free(roots); free(exact_roots); free(poly_quad); free(roots_quad);
	}
	// close file
	fclose(f);
	// return
	return 0;
}