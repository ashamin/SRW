//#include "minres5dOclSSOR.h"
//
//#include "stdlib.h"
//
//#define PTR_FLAG CL_MEM_COPY_HOST_PTR
//
//minres5dOclSSOR::minres5dOclSSOR(MathArea2d* const area, const SSORpar* const precond, const double epsilon, const int maxit)
//{
//	ap = area->getAp();
//    an = area->getAn();
//    as = area->getAs();
//    ae = area->getAe();
//    aw = area->getAw();
//    f = area->getF();
//    
//    dx_d = precond->dx_d;
//    dx_l = precond->dx_l;
//    dx_u = precond->dx_u;
//    dy_d = precond->dy_d;
//    dy_l = precond->dy_l;
//    dy_u = precond->dy_u;
//    
//    this->epsilon = epsilon;
//    this->maxit = maxit;
//    
//    m = area->getN();
//    
//    x = new double[m];
//
//    for (int i = 0; i<m; i++){
//		x[i] = (double)i/m;
//    }
//
//    corr = new double[m];
//    r = new double[m];
//    Aw = new double[m];
//    
//    ixs = area->getI() - 2;
//
//	cl_platform_id platform_id = NULL;
//    cl_device_id device_id = NULL;
//    cl_context context = NULL;
//    cl_command_queue command_queue = NULL;
//    cl_program program = NULL;
//    cl_kernel kernel = NULL;
//}
//
//minres5dOclSSOR::~minres5dOclSSOR(){
//    delete [] corr;
//    delete [] r;
//    delete [] Aw;
//	//ocl_clean();
//}
//
//double minres5dOclSSOR::exec_time()
//{
//    return time;
//}
//
//int minres5dOclSSOR::it_num()
//{
//    return maxit;
//}
//
//cl_int minres5dOclSSOR::setArgsToKernel(cl_kernel kernel, cl_mem *pars, int parnum){
//	cl_int ret;
//	for (int i = 0; i<parnum; i++)
//		ret |= clSetKernelArg(kernel, i, sizeof(cl_mem), (void *)&pars[i]);
//	return ret;
//}
//
//cl_int minres5dOclSSOR::reSetArgToKernel(cl_kernel kernel, cl_mem par, int parnum){
//	return clSetKernelArg(kernel, parnum, sizeof(cl_mem), (void *)&par);
//}
//
//double* minres5dOclSSOR::solve()
//{  
//    // due to OpenMP does not allow class members in OpenMP clauses
//    double* ap = this->ap;
//    double* as = this->as;
//    double* an = this->an;
//    double* ae= this->ae;
//    double* aw = this->aw;
//    
//    double* dx_d = this->dx_d;
//    double* dx_l = this->dx_l;
//    double* dx_u = this->dx_u;
//    double* dy_d = this->dy_d;
//    double* dy_l = this->dy_l;
//    double* dy_u = this->dy_u;
//    
//    double* x = this->x;
//    double* corr = this->corr;
//    double* r = this->r;
//    double* f = this->f;
//    double* Aw = this->Aw;
//    
//    double epsilon = this->epsilon;
//    int ixs = this->ixs;
//    int m = this->m; 
//    
//    double* tmp_v = new double[m];
//    
//    double max_it_local = maxit;
//    maxit = 0;
//    // tau parameter
//    double tau = 1;
//    // norm of residual vector r
//    double rnorm = 0;
//    double tmp = 0;
//    // (Aw, r) and (Aw, Aw) dot products
//    double dp_Aw_r = 0, dp_Aw_Aw = 0;
//
//    // matrix always is square so we can split this way
//    int n = sqrt(m), i, j, k = 0;
//    
//    int THREAD_NUM = omp_get_max_threads();
//
//    omp_set_dynamic(0);
//    omp_set_num_threads(THREAD_NUM);
//
//
//	using namespace std;
//	cl_mem a, b, c;
//
//	FILE *fp;
//    char *program_buffer;
//    
//	FILE *program_handle = fopen("mincorr_kernels.cl", "rb");
//	fseek(program_handle, 0, SEEK_END);
//	int program_size = ftell(program_handle);
//	rewind(program_handle);
//	program_buffer = (char*)malloc(program_size+1);
//	program_buffer[program_size] = '\0';
//	fread(program_buffer, sizeof(char),	program_size, program_handle);
//	fclose(program_handle);
//
//	ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
//    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, &ret_num_devices);
//
//    context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
//
//    command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
//	
//	cl_mem fBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), f, &ret);
//	cl_mem asBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), as, &ret);
//	cl_mem apBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), ap, &ret);
//	cl_mem anBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), an, &ret);
//	cl_mem aeBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), ae, &ret);
//	cl_mem awBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), aw, &ret);
//	cl_mem xBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), x, &ret);
//	cl_mem ixsBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, sizeof(int), &ixs, &ret);
//	cl_mem rBUFF = clCreateBuffer(context, CL_MEM_WRITE_ONLY, m*sizeof(double), NULL, &ret);
//	
//	cl_mem AwBUFF = clCreateBuffer(context, CL_MEM_WRITE_ONLY, m*sizeof(double), NULL, &ret);
//	cl_mem corrBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), corr, &ret);
//
//	cl_mem dx_lBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), dx_l, &ret);
//	cl_mem dx_dBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), dx_d, &ret);
//	cl_mem dx_uBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), dx_u, &ret);
//	cl_mem dy_lBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), dy_l, &ret);
//	cl_mem dy_dBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), dy_d, &ret);
//	cl_mem dy_uBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), dy_u, &ret);
//	cl_mem nBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, sizeof(int), &n, &ret);
//
//	double *loc_c = new double[m];
//	double *loc_d = new double[m];
//
//	cl_mem loc_cBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), loc_c, &ret);
//	cl_mem loc_dBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), loc_d, &ret);
//
//	cl_mem solveBUFF = clCreateBuffer(context, CL_MEM_WRITE_ONLY, m*sizeof(double), NULL, &ret);
//	// r is because whatever
//	cl_mem tmpBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), r, &ret);
//
//	program = clCreateProgramWithSource(context, 1, (const char **)&program_buffer, (const size_t *)&program_size, &ret);
//    ret     = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
//
//	cl_mem *residual_comp_pars = new cl_mem[9];
//	residual_comp_pars[0] = fBUFF;
//	residual_comp_pars[1] = asBUFF;
//	residual_comp_pars[2] = apBUFF;
//	residual_comp_pars[3] = anBUFF;
//	residual_comp_pars[4] = aeBUFF;
//	residual_comp_pars[5] = awBUFF;
//	residual_comp_pars[6] = xBUFF;
//	residual_comp_pars[7] = ixsBUFF;
//	residual_comp_pars[8] = rBUFF;
//
//	cl_kernel residual_comp_kernel = clCreateKernel(program, "computeResidual", &ret);
//
//	ret = setArgsToKernel(residual_comp_kernel, residual_comp_pars, 9);
//
//	cl_mem *Aw_comp_pars = new cl_mem[8];
//	Aw_comp_pars[0] = asBUFF;
//	Aw_comp_pars[1] = apBUFF;
//	Aw_comp_pars[2] = anBUFF;
//	Aw_comp_pars[3] = aeBUFF;
//	Aw_comp_pars[4] = awBUFF;
//	Aw_comp_pars[5] = corrBUFF;
//	Aw_comp_pars[6] = ixsBUFF;
//	Aw_comp_pars[7] = AwBUFF;
//
//	cl_kernel Aw_comp_kernel = clCreateKernel(program, "computeAw", &ret);
//
//	ret = setArgsToKernel(Aw_comp_kernel, Aw_comp_pars, 8);
//
//	cl_mem *TDMA_pars_X = new cl_mem[8];
//	TDMA_pars_X[0] = dx_lBUFF;
//	TDMA_pars_X[1] = dx_dBUFF;
//	TDMA_pars_X[2] = dx_uBUFF;
//	TDMA_pars_X[3] = solveBUFF;
//	TDMA_pars_X[4] = rBUFF;
//	TDMA_pars_X[5] = nBUFF;
//	TDMA_pars_X[6] = loc_cBUFF;
//	TDMA_pars_X[7] = loc_dBUFF;	
//
//	cl_mem *TDMA_pars_Y = new cl_mem[8];
//	TDMA_pars_Y[0] = dy_lBUFF;
//	TDMA_pars_Y[1] = dy_dBUFF;
//	TDMA_pars_Y[2] = dy_uBUFF;
//	TDMA_pars_Y[3] = solveBUFF;
//	TDMA_pars_Y[4] = NULL;
//	TDMA_pars_Y[5] = nBUFF;
//	TDMA_pars_Y[6] = loc_cBUFF;
//	TDMA_pars_Y[7] = loc_dBUFF;	
//
//	cl_kernel TDMA_kernel = clCreateKernel(program, "TDMA", &ret);
//
//	ret = setArgsToKernel(TDMA_kernel, TDMA_pars_X, 8);
//
//
//
//
//
//	size_t global_item_size = m;
//    
//    time = omp_get_wtime();
//
//    //while(maxit++ < max_it_local){
//    while (1){
//		maxit++;
//		cout << maxit << endl;
//    
//		//computing r and norm(r)
//		//r = f - A*x
//		rnorm = 0;
//
//		clReleaseMemObject(xBUFF);
//		xBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), x, &ret);
//		residual_comp_pars[6] = xBUFF;
//
//		ret = reSetArgToKernel(residual_comp_kernel, xBUFF, 6);
//
//		ret = clEnqueueNDRangeKernel(command_queue, residual_comp_kernel, 1, NULL, 
//								 &global_item_size, NULL, 0, NULL, NULL);
//
//		//ret = clEnqueueBarrier(command_queue);
//
//		ret = clEnqueueReadBuffer(command_queue, rBUFF, CL_TRUE, 0, m*sizeof(double), r, 0, NULL, NULL);
//
//
//		for (i = 0; i<m; i++){
//			tmp = fabs(r[i]);
//			if (rnorm < tmp) rnorm = tmp;
//		}
//
//
//		global_item_size = n;
//
//		clReleaseMemObject(tmpBUFF);
//		tmpBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), r, &ret);
//		TDMA_pars_X[4] = tmpBUFF;
//
//		ret = setArgsToKernel(TDMA_kernel, TDMA_pars_X, 8);
//
//		ret = clEnqueueNDRangeKernel(command_queue, TDMA_kernel, 1, NULL, 
//								 &global_item_size, NULL, 0, NULL, NULL);
//
//		//ret = clEnqueueBarrier(command_queue);
//
//		ret = clEnqueueReadBuffer(command_queue, solveBUFF, CL_TRUE, 0, m*sizeof(double), tmp_v, 0, NULL, NULL);
//
//		clReleaseMemObject(tmpBUFF);
//		tmpBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), tmp_v, &ret);
//		TDMA_pars_Y[4] = tmpBUFF;
//
//		ret = setArgsToKernel(TDMA_kernel, TDMA_pars_Y, 8);
//
//		ret = clEnqueueNDRangeKernel(command_queue, TDMA_kernel, 1, NULL, 
//								 &global_item_size, NULL, 0, NULL, NULL);
//
//		//ret = clEnqueueBarrier(command_queue);
//
//		ret = clEnqueueReadBuffer(command_queue, solveBUFF, CL_TRUE, 0, m*sizeof(double), corr, 0, NULL, NULL);
//
//
//		global_item_size = m;
//
//
//
//		//computing Aw = A*r and dot products (Aw, r) and (Aw, Aw)
//		dp_Aw_r = 0;
//		dp_Aw_Aw = 0;
//
//		clReleaseMemObject(corrBUFF);
//		corrBUFF = clCreateBuffer(context, CL_MEM_READ_ONLY | PTR_FLAG, m*sizeof(double), corr, &ret);
//
//		ret = reSetArgToKernel(Aw_comp_kernel, corrBUFF, 5);
//
//		ret = clEnqueueNDRangeKernel(command_queue, Aw_comp_kernel, 1, NULL, 
//									&global_item_size, NULL, 0, NULL, NULL);
//
//		//ret = clEnqueueBarrier(command_queue);
//
//		ret = clEnqueueReadBuffer(command_queue, AwBUFF, CL_TRUE, 0, m*sizeof(double), Aw, 0, NULL, NULL);
//
//		for (k = 0; k<m; k++){
//			dp_Aw_r += Aw[k]*r[k];
//			dp_Aw_Aw += Aw[k]*Aw[k];
//		}
//
//		//computing tau
//		tau = dp_Aw_r / dp_Aw_Aw;
//
//		//computes new approximation of x
//	#pragma omp parallel for shared(corr, x) \
//	firstprivate(m, tau) private(k) \
//	schedule(static)
//    
//		for (k = 0; k<m; k++)
//			x[k] += corr[k]*tau;
//
//		if (rnorm < epsilon) break;
//		//if (maxit == max_it_local)
//		//  std::cout << "Iteration process obviously won't converge. \\n Try to increase \" maxit \" value" << std::endl;
//    }
//    
//    time = omp_get_wtime() - time;
//    
//    delete [] tmp_v;
//
//    return x;
//
//}
//
//int minres5dOclSSOR::ocl_clean(){
//	ret = clFlush(command_queue);
//    ret = clFinish(command_queue);
//    ret = clReleaseKernel(kernel);
//    ret = clReleaseProgram(program);
//    ret = clReleaseCommandQueue(command_queue);
//    ret = clReleaseContext(context);
//	return 0;
//}
