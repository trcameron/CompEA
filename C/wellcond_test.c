/* Well-Conditioned Test */
#include "ehrlich_aberth.h"
#include "error.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/* Main Function */
int main(int argc,char **argv)
{
	/* testing variables */
	const unsigned deg_min = 10, deg_max = 1280, itmax = 50, itnum = 10; 
	double ea_et, comp_et, quad_et;
	clock_t ea_begin, ea_end, comp_begin, comp_end, quad_begin, quad_end;
	/* random coefficients test */
	FILE* f;
	f = fopen("data_files/wellcond_rand_test.dat","w+");
	fprintf(f,"degree, ea_time, ea_err, ea_comp_time, ea_comp_err, ea_quad_time, ea_quad_err\n");
	ea_et = 0; comp_et = 0; quad_et = 0;
	for(unsigned int deg = deg_min; deg <= deg_max; deg *= 2)
	{
		// initialize storage for error
		double* ea_err = (double*)malloc(itnum*sizeof(double));
		double* comp_err = (double*)malloc(itnum*sizeof(double));
		double* quad_err = (double*)malloc(itnum*sizeof(double));
		for(unsigned int it = 0; it < itnum; it++)
		{
			// initialize random polynomial and storage for roots
			double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
			double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
			srand(time(NULL));
			for(int i=0; i<=deg; i++)
			{
				double real = (double)2*rand()/RAND_MAX - 1;
				double imag = (double)2*rand()/RAND_MAX - 1;
				poly[i] = real + I*imag;
			}
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
			ea_begin = clock();
			ehrlich_aberth(poly,roots,deg,itmax);
			ea_end = clock();
			ea_et += (double)(ea_end - ea_begin) / CLOCKS_PER_SEC;
			ea_err[it] = back_err(poly,roots,deg);
			comp_begin = clock();
			ehrlich_aberth_comp(poly,roots,deg,itmax);
			comp_end = clock();
			comp_et += (double)(comp_end - comp_begin) / CLOCKS_PER_SEC;
			comp_err[it] = back_err(poly,roots,deg);
			quad_begin = clock();
			ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
			quad_end = clock();
			quad_et += (double)(quad_end - quad_begin) / CLOCKS_PER_SEC;
			for(int i=0; i<deg; i++)
			{
				roots[i] = mpc_get_dc(roots_quad[i],MPC_RNDNN);
			}
			quad_err[it] = back_err(poly,roots,deg);
			// free memory
			free(poly); free(roots); free(poly_quad); free(roots_quad);	
		}
		// write to file
		fprintf(f,"%d, ",deg);
		fprintf(f,"%.5e, ",ea_et);
		fprintf(f,"%.5e, ",max_value(ea_err,itnum));
		fprintf(f,"%.5e, ",comp_et);
		fprintf(f,"%.5e, ",max_value(comp_err,itnum));
		fprintf(f,"%.5e, ",quad_et);
		fprintf(f,"%.5e\n",max_value(quad_err,itnum));
		// free memory
		free(ea_err); free(comp_err); free(quad_err);
	}
	// close file
	fclose(f);
	/* natural coefficients test */
	f = fopen("data_files/wellcond_nat_test.dat","w+");
	fprintf(f,"degree, ea_time, ea_err, ea_comp_time, ea_comp_err, ea_quad_time, ea_quad_err\n");
	ea_et = 0; comp_et = 0; quad_et = 0;
	for(unsigned int deg = deg_min; deg <= deg_max; deg *= 2)
	{
		// initialize storage for error
		double* ea_err = (double*)malloc(itnum*sizeof(double));
		double* comp_err = (double*)malloc(itnum*sizeof(double));
		double* quad_err = (double*)malloc(itnum*sizeof(double));
		for(unsigned int it = 0; it < itnum; it++)
		{
			// initialize random polynomial and storage for roots
			double complex* poly = (double complex*)malloc((deg+1)*sizeof(double complex));
			double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
			srand(time(NULL));
			for(int i=0; i<=deg; i++)
			{
				poly[i] = (i+1);
			}
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
			ea_begin = clock();
			ehrlich_aberth(poly,roots,deg,itmax);
			ea_end = clock();
			ea_et += (double)(ea_end - ea_begin) / CLOCKS_PER_SEC;
			ea_err[it] = back_err(poly,roots,deg);
			comp_begin = clock();
			ehrlich_aberth_comp(poly,roots,deg,itmax);
			comp_end = clock();
			comp_et += (double)(comp_end - comp_begin) / CLOCKS_PER_SEC;
			comp_err[it] = back_err(poly,roots,deg);
			quad_begin = clock();
			ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
			quad_end = clock();
			quad_et += (double)(quad_end - quad_begin) / CLOCKS_PER_SEC;
			for(int i=0; i<deg; i++)
			{
				roots[i] = mpc_get_dc(roots_quad[i],MPC_RNDNN);
			}
			quad_err[it] = back_err(poly,roots,deg);
			// free memory
			free(poly); free(roots); free(poly_quad); free(roots_quad);	
		}
		// write to file
		fprintf(f,"%d, ",deg);
		fprintf(f,"%.5e, ",ea_et);
		fprintf(f,"%.5e, ",max_value(ea_err,itnum));
		fprintf(f,"%.5e, ",comp_et);
		fprintf(f,"%.5e, ",max_value(comp_err,itnum));
		fprintf(f,"%.5e, ",quad_et);
		fprintf(f,"%.5e\n",max_value(quad_err,itnum));
		// free memory
		free(ea_err); free(comp_err); free(quad_err);
	}
	// close file
	fclose(f);
	// return
	return 0;
}