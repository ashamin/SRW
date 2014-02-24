#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void computeResidual(__global double *f, __global double *as, __global double *ap,
                                __global double *an, __global double *ae, __global double *aw,
                                __global double *x, __global int *ixs, __global double *r)
{
    int k = get_global_id(0);
    r[k] = f[k] - as[k]*x[k-1] - ap[k]*x[k] - an[k]*x[k+1] - ae[k]*x[k+*ixs]
            - aw[k]*x[k-*ixs];
}

__kernel void computeAw(__global double *as, __global double *ap, __global double *an, 
                                __global double *ae, __global double *aw, __global double *corr, 
                                __global int *ixs, __global double *Aw)
{
    int k = get_global_id(0);
    Aw[k] = as[k]*corr[k-1] + ap[k]*corr[k] + an[k]*corr[k+1] + ae[k]*corr[k+*ixs]
            + aw[k]*corr[k-*ixs];
}

__kernel void TDMA(__global double *a, __global double *b, __global double *c, 
                                __global double *x, __global double *d, __global int *n, 
                                __global double *loc_c, __global double *loc_d)
{
    int k = get_global_id(0)**n;
    
	loc_c[k] = c[k] / b[k];
	loc_d[k] = d[k] / b[k];

	for (int i = k+1; i<k+*n; i++){
		double tmp = b[i] - loc_c[i-1]*a[i];
		loc_c[i] = c[i] / tmp;
		loc_d[i] = (d[i] - loc_d[i-1]*a[i]) / tmp;
	}

	x[k+*n-1] = loc_d[k+*n-1];
	for (int i = k+*n-2; i>=k; i--)
		x[i] = loc_d[i] - loc_c[i]*x[i+1];
    
}
