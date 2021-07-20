/* Testing for Ehrlich Aberth */
#include "ehrlich_aberth.h"
#include "error.h"
/* Main Function */
int main(int argc,char **argv)
{
	int deg = 4, itmax = 30;
	double complex poly[5] = {1,5,2,7,1};
	double complex* roots = (double complex*)malloc(deg*sizeof(double complex));
	ehrlich_aberth_comp(poly,roots,deg,itmax);
	for(int i=0; i<deg; i++)
	{
		printf("roots[i] = %f + i%f\n",creal(roots[i]),cimag(roots[i]));
	}
	
	mpc_t* poly_quad = (mpc_t*)malloc((deg+1)*sizeof(mpc_t));
	mpc_t* roots_quad = (mpc_t*)malloc(deg*sizeof(mpc_t));
	for(int i=0; i<deg; i++)
	{
		mpc_init2(poly_quad[i],113); mpc_init2(roots_quad[i],113);
		mpc_set_d(poly_quad[i],poly[i],MPC_RNDNN);
	}
	mpc_init2(poly_quad[deg],113);
	mpc_set_d(poly_quad[deg],poly[deg],MPC_RNDNN);
	ehrlich_aberth_quad(poly_quad,roots_quad,deg,itmax);
	
	mpfr_t abs_diff;
	mpc_t diff;
	mpfr_init2(abs_diff,113);
	mpc_init2(diff,113);
	for(int i=0; i<deg; i++)
	{
		mpc_set_dc(diff,roots[i],MPC_RNDNN);
		mpc_sub(diff,roots_quad[i],diff,MPC_RNDNN);
		mpc_abs(abs_diff,diff,MPFR_RNDN);
		printf("err[i] = %.4e\n",mpfr_get_d(abs_diff,MPFR_RNDN));
	}
	
	printf("comp backward error = %.4e\n",back_err(poly,roots,deg));
	
	free(roots);
	free(poly_quad);
	free(roots_quad);
	return 0;
}