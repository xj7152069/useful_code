#define DEBUG_MAIN

#ifdef	DEBUG_MAIN
	#include<stdio.h>
	#include<stdlib.h>
	#include<math.h>
	#include<omp.h>
	#define PI (3.141592653589793)
	#include "fballoc.h"
	#include "fballoc.c"
	#include "spline_dy_ddy.c"
	#include "spline_z_dz_ddz.c"
#endif

void	splineInterpolation_1d(
		float	*in,		// in[n1].
		int	n1_in,
		float	d1_in,
		float	f1_in,
		float	*out,		// out[n1].
		int	n1_out,
		float	d1_out,
		float	f1_out
		)
{
	float	*z_node		= alloc1float (n1_in);
	float	*v_node		= alloc1float (n1_in);
	float	*dz_node	= alloc1float (n1_in);
	float	*s_spl_tmp	= alloc1float (n1_in);

	for ( int i1 = 0; i1 < n1_in; ++i1 )
	{
		z_node[i1] = i1 * d1_in + f1_in;
		v_node[i1] = in[i1];
	}

	// Init.
	dz_node[0] = (v_node[1] - v_node[0]) / (z_node[1] - z_node[0]);
	dz_node[n1_in-1] = (v_node[n1_in-1] - v_node[n1_in-2]) / (z_node[n1_in-1] - z_node[n1_in-2]);

	spl_dy (z_node, v_node, n1_in, dz_node, s_spl_tmp);

	for ( int i1 = 0; i1 < n1_out; ++i1 )
	{
		float	time	= f1_out + d1_out*(float)i1;
		spl_z (z_node, v_node, n1_in, dz_node, time, &out[i1], s_spl_tmp);
	}

	// free buffers.
	free(z_node);
	free(v_node);
	free(dz_node);
	free(s_spl_tmp);
}


void	splineInterpolation_2d(
		float	**in,		// in[n2][n1].
		int	n1_in,		// the fast dimension.
		float	d1_in,
		float	f1_in,
		int	n2_in,		// the second dimension.
		float	d2_in,
		float	f2_in,
		float	**out,		// out[n2][n1].
		int	n1_out,
		float	d1_out,
		float	f1_out,
		int	n2_out,
		float	d2_out,
		float	f2_out
		)
{
	// step 1: interpolation along the fast-dim.
	float	**tmp1	= alloc2float(n1_out, n2_in);	// temp-result, smooth along the fast-dim.
	#pragma omp parallel for
	for ( int i2 = 0 ; i2 < n2_in; ++ i2 )
	{
		splineInterpolation_1d( &in[i2][0], n1_in, d1_in, f1_in, &tmp1[i2][0], n1_out, d1_out, f1_out);
	}

	// step 2: interpolation along the 2nd-dim.
	float	**tmp2_in	= alloc2float(n2_in,  n1_out);	// transpose of "tmp1".
	float	**tmp2_out	= alloc2float(n2_out, n1_out);	// temp-result, smooth along the 2nd-dim.
	#pragma omp parallel for
	for ( int i1 = 0 ; i1 < n1_out; ++ i1 )
	{
		for ( int i2 = 0 ; i2 < n2_in; ++ i2 )
			tmp2_in[i1][i2]	= tmp1[i2][i1];		// transpose.

		// interpolation along the 2nd-dim.
		splineInterpolation_1d( &tmp2_in[i1][0], n2_in, d2_in, f2_in, &tmp2_out[i1][0], n2_out, d2_out, f2_out);

		for ( int i2 = 0 ; i2 < n2_out; ++ i2 )
			out[i2][i1]	= tmp2_out[i1][i2];	// transpose.
	}

	free2float(tmp1);
	free2float(tmp2_in);
	free2float(tmp2_out);
}

void	splineInterpolation_3d(
		float	***in,		// in[n3][n2][n1].
		int	n1_in,		// the fast dimension.
		float	d1_in,
		float	f1_in,
		int	n2_in,		// the 2nd-dimension.
		float	d2_in,
		float	f2_in,
		int	n3_in,		// the 3rd-dimension.
		float	d3_in,
		float	f3_in,
		float	***out,		// out[n3][n2][n1].
		int	n1_out,
		float	d1_out,
		float	f1_out,
		int	n2_out,
		float	d2_out,
		float	f2_out,
		int	n3_out,
		float	d3_out,
		float	f3_out
		)
{
	// step 1: 2d-interpolation along the n1 x n2 plane.
	float	***tmp1	= alloc3float(n1_out, n2_out, n3_in);	// temp-result, smooth along the fast-dim.
	for ( int i3 = 0 ; i3 < n3_in; ++ i3 )
	{
		splineInterpolation_2d(in[i3], n1_in, d1_in, f1_in, n2_in, d2_in, f2_in,
			tmp1[i3], n1_out, d1_out, f1_out, n2_out, d2_out, f2_out);
	}

	// step 2: 1d-interpolation along the 3rd-dimension.
	int	n12	= n1_out*n2_out;
	float	**tmp2_in	= alloc2float(n3_in, n12);
	float	**tmp2_out	= alloc2float(n3_out, n12);
	// feed data.
	#pragma omp parallel for
	for ( int i2 = 0 ; i2 < n2_out; ++ i2 )
	for ( int i1 = 0 ; i1 < n1_out; ++ i1 )
	{
		int	i12	= i2*n1_out + i1;
		for ( int i3 = 0 ; i3 < n3_in; ++ i3 )
			tmp2_in[i12][i3]	= tmp1[i3][i2][i1];
	}

	#pragma omp parallel for
	for ( int i12 = 0 ; i12 < n12; ++ i12 )
	{
		splineInterpolation_1d( &tmp2_in[i12][0], n3_in, d3_in, f3_in, &tmp2_out[i12][0], n3_out, d3_out, f3_out);
	}

	// fetch data.
	#pragma omp parallel for
	for ( int i2 = 0 ; i2 < n2_out; ++ i2 )
	for ( int i1 = 0 ; i1 < n1_out; ++ i1 )
	{
		int	i12	= i2*n1_out + i1;
		for ( int i3 = 0 ; i3 < n3_out; ++ i3 )
			out[i3][i2][i1]	= tmp2_out[i12][i3];
	}

	free3float(tmp1);
	free2float(tmp2_in);
	free2float(tmp2_out);
}

#ifdef	DEBUG_MAIN
#define FILE_NAME_MAX_LENGTH 1024
#include "fbrw.h"
#include "fbrw.c"
int	main()
{
/*
	// for Overthrust3D model.
	int	nz	= 251;
	int	nx	= 801;
	int	ny	= 801;
	float	dz	= 25.0;
	float	dx	= 25.0;
	float	dy	= 25.0;
	float	vz0	= 0;
	float	vx0	= 0;
	float	vy0	= 0;

	int	nresamp	= 2;
	int	nz_1	= (int)(nz/nresamp + 1.);
	int	nx_1	= (int)(nx/nresamp + 1.);
	int	ny_1	= (int)(ny/nresamp + 1.);
	float	dz_1	= dz*(float)nresamp;
	float	dx_1	= dx*(float)nresamp;
	float	dy_1	= dy*(float)nresamp;
	float	vz0_1	= 0;
	float	vx0_1	= 0;
	float	vy0_1	= 0;
	printf("Dense-grid:  nx=%d ny=%d nz=%d\n", nx, ny, nz);
	printf("Coarse-grid: nx=%d ny=%d nz=%d\n", nx_1, ny_1, nz_1);

	float	***vtrue_3d	= alloc3float(nz, nx, ny);
	float	***vdense_3d	= alloc3float(nz, nx, ny);
	float	***vcoarse_3d	= alloc3float(nz_1, nx_1, ny_1);

	// define the path for the test model.
	char	path[FILE_NAME_MAX_LENGTH]="/ssd/ssd_3.2T/fb/testdata/Overthrust3D/velocity";
	char	fnvtrue_3d[FILE_NAME_MAX_LENGTH]="";
	char	fnvcoarse_3d[FILE_NAME_MAX_LENGTH]="";
	char	fnvdense_3d[FILE_NAME_MAX_LENGTH]="";
	sprintf(fnvtrue_3d,   "%s/overthrust_nz251_nx801_ny801.true.dat",       path);
	sprintf(fnvcoarse_3d, "%s/overthrust_nz126_nx401_ny401.true.dx50m.dat", path);
	sprintf(fnvdense_3d,  "%s/overthrust_nz251_nx801_ny801.interp3d.dat",   path);

	read_3d_float_rb(vtrue_3d, ny, nx, nz, fnvtrue_3d);

	for ( int iy = 0 ; iy < ny;  iy += nresamp )
	for ( int ix = 0 ; ix < nx;  ix += nresamp )
	for ( int iz = 0 ; iz < nz;  iz += nresamp )
		vcoarse_3d[iy/2][ix/2][iz/2] = vtrue_3d[iy][ix][iz];
	write_3d_float_wb(vcoarse_3d, ny_1, nx_1, nz_1, fnvcoarse_3d);

	// test of 3-D interpolation.
	splineInterpolation_3d(
		vcoarse_3d,
		nz_1, dz_1, vz0_1,
		nx_1, dx_1, vx0_1,
		ny_1, dy_1, vy0_1,
		vdense_3d,
		nz, dz, vz0,
		nx, dx, vx0,
		ny, dy, vy0);
	write_3d_float_wb(vdense_3d, ny, nx, nz, fnvdense_3d);
*/
	int	nz_1	= 3010;
	int	nx_1	= 233;
	int	ny_1	= 337;
	float	dz_1	= 1.0;
	float	dx_1	= 20.0;
	float	dy_1	= 10.0;
	float	vz0_1	= 0;
	float	vx0_1	= 493.0;
	float	vy0_1	= 35.0;

	int	nz= 3010;
	int	nx= 2000;
	int	ny= 800;
	float	dz= 1.0;
	float	dx= 1.0;
	float	dy= 1.0;
	float	vz0= 0;
	float	vx0= 2001;
	float	vy0= 1801;
	printf("Dense-grid:  nx=%d ny=%d nz=%d\n", nx, ny, nz);
	printf("Coarse-grid: nx=%d ny=%d nz=%d\n", nx_1, ny_1, nz_1);
	// define the path for the gradient.
	//char	path[FILE_NAME_MAX_LENGTH]="/ssd/ssd_raid5/fb/LJ2017_cmp/velocity";
	char	path[FILE_NAME_MAX_LENGTH]="/data2_2/fb/LJ2017_cmp/velocity";
	char	fngrad_CoarseGrid[FILE_NAME_MAX_LENGTH]="";
	char	fngrad_DenseGrid[FILE_NAME_MAX_LENGTH]="";
	sprintf(fngrad_CoarseGrid, "%s/LJ_rms_new.dat", path);
	sprintf(fngrad_DenseGrid,  "%s/rmsVel_ncdp2000_nline800_ntau3010.dat",  path);

	float	***grad_coarse	= alloc3float(nz_1, nx_1, ny_1);
	float	***grad_dense	= alloc3float(nz, nx, ny);
	printf("Now, reading file: %s\n", fngrad_CoarseGrid);
	read_3d_float_rb(grad_coarse, ny_1, nx_1, nz_1, fngrad_CoarseGrid);

	// test of 3-D interpolation.
	printf("Now, doing 3-D interpolation. \n");
	splineInterpolation_3d(
		grad_coarse,
		nz_1, dz_1, vz0_1,
		nx_1, dx_1, vx0_1,
		ny_1, dy_1, vy0_1,
		grad_dense,
		nz, dz, vz0,
		nx, dx, vx0,
		ny, dy, vy0);

	printf("Now, writing file: %s\n", fngrad_DenseGrid);
	write_3d_float_wb(grad_dense, ny, nx, nz, fngrad_DenseGrid);

	return	0;

}
#endif
