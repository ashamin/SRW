//#ifndef MINRES_5D_OCL_SSOR_H
//#define MINRES_5D_OCL_SSOR_H
//
//#include "minres5d.h"
//#include "MathArea2d.h"
//#include "SSORpar.h"
//#include "omp.h"
//#include "CL\cl.h"
//
//#include <iostream>
//#include "math.h"
//
//class minres5dOclSSOR : public minres5d
//{
//public:
//	minres5dOclSSOR(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit);
//	virtual ~minres5dOclSSOR();
//	double* solve();
//  
//	virtual double exec_time();
//	int it_num();
//
//private:
//	/** Main diagonal of x related part of parallel tridiagonal SSOR preconditioner */
//	double* dx_d;
//	/** Lower diagonal of x */
//	double* dx_l;
//	/** Upper diagonal of x */
//	double* dx_u;
//	  /** Main diagonal of y */
//	double* dy_d;
//	  /** Lower diagonal of y */
//	double* dy_l;
//	  /** Upper diagonal of y */
//	double* dy_u;
//
//	/** OpenCL variables **/
//	cl_platform_id platform_id;
//    cl_device_id device_id;
//    cl_context context;
//    cl_command_queue command_queue;
//    cl_program program;
//    cl_kernel kernel;
//    cl_uint ret_num_devices;
//    cl_uint ret_num_platforms;
//    cl_int ret;
//
//	int ocl_clean();
//	cl_int setArgsToKernel(cl_kernel kernel, cl_mem *pars, int parnum);
//	cl_int reSetArgToKernel(cl_kernel kernel, cl_mem par, int parnum);
//  
//  
//	/** Execution time */
//	double time; 
//  
//};
//
//#endif
