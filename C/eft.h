#ifndef EFT
#define EFT
#include <complex.h>
/* EFT Data Structure */
struct eft
{
	double fl_res, fl_err;
};
/* EFT Cmplx Data Structure for Sum */
struct eft_cmplx_sum
{
	double complex fl_res, fl_err;
};
/* EFT Cmplx Data Structure for Product */
struct eft_cmplx_prod
{
	double complex fl_res, fl_err1, fl_err2, fl_err3;
};
#endif