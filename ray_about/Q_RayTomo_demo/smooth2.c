#include<stdlib.h>
#include<stdio.h>

#include "include/alloc.h"
#include "include/alloc.c"

void tripd(float *d, float *e, float *b, int n)
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

int smooth2(int n1, int n2, float r1, float r2, float **v0)
{

	int nmax;	/* max of n1 and n2 */
	int ix, iz;	/* counters */
	int win[4];
	float **v;	/* array of output velocities */
	float **w;	/* intermediate array */
	float *d, *e;	/* input arrays for subroutine tripd */
	float *f;	/* intermediate array */
	float rw;	/* smoothing parameter for window */

	/* scale the smoothing parameter */
	r1 = r1*r1*0.25;
	r2 = r2*r2*0.25;

	/* allocate space */
	nmax = (n1<n2)?n2:n1;
	v = alloc2float(n1,n2);
	w = alloc2float(n1,n2);
	d = alloc1float(nmax);
	e = alloc1float(nmax);
	f = alloc1float(nmax);

	/* save the original velocity */
        for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			v[ix][iz]=v0[ix][iz];

	/* get parameters for window function */

	win[0] = 0;
	win[1] = n1;
	win[2] = 0;
	win[3] = n2;

	rw = 0.;
	rw = rw*rw*0.25;
 
	/* define the window function */
	for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			w[ix][iz] = 0;	
	for(ix=win[2]; ix<win[3]; ++ix)
	 	for(iz=win[0]; iz<win[1]; ++iz)
			w[ix][iz] = 1;	

	if(win[0]>0 || win[1]<n1 || win[2]>0 || win[3]<n2){
	/*	smooth the window function */
         	for(iz=0; iz<n1; ++iz){
	 		for(ix=0; ix<n2; ++ix){
				d[ix] = 1.0+2.0*rw;
				e[ix] = -rw;
				f[ix] = w[ix][iz];
			}
        		d[0] -= rw;
         		d[n2-1] -= rw;
         		tripd(d,e,f,n2);
	 		for(ix=0; ix<n2; ++ix)
				w[ix][iz] = f[ix];
		}
         	for(ix=0; ix<n2; ++ix){
	 		for(iz=0; iz<n1; ++iz){
				d[iz] = 1.0+2.0*rw;
				e[iz] = -rw;
				f[iz] = w[ix][iz];
		}
        		d[0] -= rw;
         		d[n1-1] -= rw;
         		tripd(d,e,f,n1);
	 		for(iz=0; iz<n1; ++iz)
				w[ix][iz] = f[iz];
		}
	}

	/*      solving for the smoothing velocity */
        for(iz=0; iz<n1; ++iz){
	 	for(ix=0; ix<n2-1; ++ix){
			d[ix] = 1.0+r2*(w[ix][iz]+w[ix+1][iz]);
			e[ix] = -r2*w[ix+1][iz];
			f[ix] = v[ix][iz];
		}
        	d[0] -= r2*w[0][iz];
         	d[n2-1] = 1.0+r2*w[n2-1][iz];
		f[n2-1] = v[n2-1][iz];
         	tripd(d,e,f,n2);
	 	for(ix=0; ix<n2; ++ix)
			v[ix][iz] = f[ix];
	}
         for(ix=0; ix<n2; ++ix){
	 	for(iz=0; iz<n1-2; ++iz){
			d[iz] = 1.0+r1*(w[ix][iz+1]+w[ix][iz+2]);
			e[iz] = -r1*w[ix][iz+2];
			f[iz] = v[ix][iz+1];
		}
		f[0] += r1*w[ix][1]*v[ix][0];
         	d[n1-2] = 1.0+r1*w[ix][n1-1];
		f[n1-2] = v[ix][n1-1];
         	tripd(d,e,f,n1-1);
	 	for(iz=0; iz<n1-1; ++iz)
			v[ix][iz+1] = f[iz];
	}

    for(ix=0; ix<n2; ++ix)
	for(iz=0; iz<n1; ++iz)
		v0[ix][iz]=v[ix][iz];
 
	free2float(v);
	free2float(w);
	free1float(d);
	free1float(e);
	free1float(f);
	
	return(0);
}

void smooth2f_(int *NVX, int *NVZ, float *R1, float *R2, float *vel)
{
	int n1, n2;
	float r1, r2;
	int ix, iz;
	
	n1=*NVZ;
	n2=*NVX;
	r1=*R1;
	r2=*R2;

	float	**v0=alloc2float(n1,n2);

	int	igrid=0;
	for(ix=0; ix<n2; ++ix)
	for(iz=0; iz<n1; ++iz)
	{
		v0[ix][iz]=vel[igrid++];
	}

	smooth2(n1, n2, r1, r2, v0);

	igrid=0;
	for(ix=0; ix<n2; ++ix)
	for(iz=0; iz<n1; ++iz)
	{
		vel[igrid++]=v0[ix][iz];
	}

	free2float(v0);

}
