#ifndef INIT_EST
#define INIT_EST
#include <complex.h>
#include <float.h>
#include <math.h>
#include <mpfr.h>
#include <mpc.h>
#include <stdlib.h>
/* point structure */
typedef struct
{
	int x;
	double y;
}	point;
/* ccs: Three points are a counter-clockwise turn if ccw > 0, clockwise if
 * ccw < 0, and collinear if ccw = 0 because ccw is a determinant that
 * gives the signed area of the triangle formed by p1, p2 and p3.
 */
static double ccw(point* p1, point* p2, point* p3)
{
	return (p2->x - p1->x)*(p3->y - p1->y) - (p2->y - p1->y)*(p3->x - p1->x);
}
/* convex_hull: Returns a list of points on the upper envelope of the 
 * convex hull in counter-clockwise order.
 */
void convex_hull(point* points,const unsigned int npoints,point* hull,unsigned int* hullsize)
{
	unsigned int k=0;
	for(int i=npoints-1; i>=0; i--)
	{
		while(k>=2 && ccw(&hull[k-2], &hull[k-1], &points[i]) <= 0) --k;
		hull[k++] = points[i];
	}
	*hullsize = k;
}
/* init_est: 
 */
void init_est(const double* alpha,const unsigned int deg,double complex* roots)
{
	// local pointers
	point* points = (point*)malloc((deg+1)*sizeof(point));
	point* hull = (point*)malloc((deg+1)*sizeof(point));
	unsigned int* hullsize = (unsigned int*)malloc(sizeof(unsigned int));
	// local variables
	const double pi2 = 6.28318530717958647693, sigma = 0.7;
	int i, j, k = 0, nzeros;
	double a1, a2, ang, r, th = pi2/deg;
	// create points
	for(i=0; i<=deg; i++)
	{
		if(alpha[i]>0)
		{
			points[i].x = i;
			points[i].y = log(alpha[i]);
 		}
		else
		{
			points[i].x = i;
			points[i].y = -1E+30;
		}
	}
	// compute convex hull
	convex_hull(points,deg+1,hull,hullsize);
	// compute initial estimates
	for(i=*hullsize-2; i>=0; i--)
	{
		nzeros = hull[i].x - hull[i+1].x;
		a1 = pow(alpha[hull[i+1].x],1.0/nzeros);
		a2 = pow(alpha[hull[i].x],1.0/nzeros);
		r = a1/a2;
		ang = pi2/nzeros;
		for(j=0; j<nzeros; j++)
		{
			roots[k+j] = r*(cos(ang*j+th*i+sigma) + I*sin(ang*j+th*i+sigma));
		}
		k += nzeros;
	}
	// free pointers
	free(points);
	free(hull);
	free(hullsize);
}
/* init_est_quad */
void init_est_quad(const mpfr_t* alpha_quad,const unsigned int deg,mpc_t* roots_quad)
{
	// local pointers
	point* points = (point*)malloc((deg+1)*sizeof(point));
	point* hull = (point*)malloc((deg+1)*sizeof(point));
	unsigned int* hullsize = (unsigned int*)malloc(sizeof(unsigned int));
	// local variables
	const double pi2 = 6.28318530717958648, sigma = 0.7;
	int i, j, k = 0, nzeros;
	double a1, a2, ang, r, th = pi2/deg;
	// create points
	for(i=0; i<=deg; i++)
	{
		if(mpfr_cmp_d(alpha_quad[i],0)>0)
		{
			points[i].x = i;
			points[i].y = log(mpfr_get_d(alpha_quad[i],MPFR_RNDN));
 		}
		else
		{
			points[i].x = i;
			points[i].y = -1E+30;
		}
	}
	// compute convex hull
	convex_hull(points,deg+1,hull,hullsize);
	// compute initial estimates
	for(i=*hullsize-2; i>=0; i--)
	{
		nzeros = hull[i].x - hull[i+1].x;
		a1 = pow(mpfr_get_d(alpha_quad[hull[i+1].x],MPFR_RNDN),1.0/nzeros);
		a2 = pow(mpfr_get_d(alpha_quad[hull[i].x],MPFR_RNDN),1.0/nzeros);
		r = a1/a2;
		ang = pi2/nzeros;
		for(j=0; j<nzeros; j++)
		{
			mpc_set_dc(roots_quad[k+j],r*(cos(ang*j+th*i+sigma) + I*sin(ang*j+th*i+sigma)),MPC_RNDNN);
		}
		k += nzeros;
	}
	// free pointers
	free(points);
	free(hull);
	free(hullsize);
}
#endif