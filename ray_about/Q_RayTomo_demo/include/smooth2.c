


void tripd(float *d, float *e, float *b, int n)
/*****************************************************************************
Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Input:
d[]	the diagonal of A 
e[]	the superdiagonal of A
b[]	the rhs of Ax=b

Output:
b[]	b[] is overwritten with the solution to Ax=b

*****************************************************************************
Notes:

Given an n-by-n symmetric, tridiagonal, positive definite matrix A and
 n-vector b, following algorithm overwrites b with the solution to Ax = b

*****************************************************************************
Author: Zhenyue Liu, Colorado School of Mines, 1993.
*****************************************************************************/
{
	int k; 
	float temp;
	
	/* decomposition */
	for(k=1; k<n; ++k){
           temp = e[k-1];
           e[k-1] = temp/d[k-1];
           d[k] -= temp*e[k-1];
	}

	/* substitution	*/
        for(k=1; k<n; ++k)  b[k] -= e[k-1]*b[k-1];
	
        b[n-1] /=d[n-1];
        for(k=n-1; k>0; --k)  b[k-1] = b[k-1]/d[k-1] - e[k-1]*b[k]; 
	
 }

 int smooth2 (float **vel, float r1, float r2, float *tripd_d, float *tripd_e, float *tripd_f, int n1, int n2)
 {
 	int ix, iz;

	for(iz=0; iz<n1; ++iz)
	{
		for(ix=0; ix<n2-1; ++ix)
		{
			tripd_d[ix] = 1.0+r2*2.0;
			tripd_e[ix] = -r2;
			tripd_f[ix] = vel[ix][iz];
		}
		tripd_d[0] -= r2;
		tripd_d[n2-1] = 1.0+r2;
		tripd_f[n2-1] = vel[n2-1][iz];
		tripd(tripd_d,tripd_e,tripd_f,n2);
		for(ix=0; ix<n2; ++ix)
		vel[ix][iz] = tripd_f[ix];
	}

	for(ix=0; ix<n2; ++ix)
	{
		for(iz=0; iz<n1-2; ++iz)
		{
			tripd_d[iz] = 1.0+r1*2.0;
			tripd_e[iz] = -r1;
			tripd_f[iz] = vel[ix][iz+1];
		}
		tripd_f[0] += r1*vel[ix][0];
		tripd_d[n1-2] = 1.0+r1;
		tripd_f[n1-2] = vel[ix][n1-1];
		tripd(tripd_d,tripd_e,tripd_f,n1-1);
		for(iz=0; iz<n1-1; ++iz)
		vel[ix][iz+1] = tripd_f[iz];
	}



	return 0;
 }
