/***************************************
*     Author:luofei
*     Date:2016-06-26 18:50
*     Filename:calculate_fatray_kernels.c
*     Description:
*
*     Last modified:2016-06-26 18:56
****************************************/
#include"./include/_lf.h"
#include"./include/_lfFunC.h"
#include "include/openCfile.h"
#include "include/calculate_fat_ray_kernels.h"
#include "include/get_vel_differential.c"
#include "include/functions_for_kernel.c"
#include "include/calculate_kernel_fresnel_norm_for_tomo.c"
#include "include/calculate_kernel_ray.c"
							  
int calculate_kernel_fresnel_(float *x_loc_ray,float *z_loc_ray,float *xp_ray,float *zp_ray,float *tau_ray, int *tflag,
							 int *n_p,int *nvx,int *nvz,float *dvx,float *dvz,float *t_period,float *i_sigma_rule,
							 float *sx_coord,float *sz_coord, float *dstep,int *itr,float *elev,float *v,float *ray,float *kernel,int *myid)
{
	int i,ix,iz,nn, n_yn, n_table,n_point;
	float ds,area_grid;
	float *q1, **vel, *p1, *q2, *p2, *deriv_2nd_src, *deriv_2nd_rec, *radius_fresnell_2,
		  **tab_v, **tab_vx, **tab_vxx;
	float **kernel_fresnel_norm, **kernel_ray;
	FILE *fp,*fp1,*fp2;
	int FileStatus,locate,TFlag;
	model_rt  m_rt;
	ray_rt    r_rt;
	ray_lamda r_lamda;

	m_rt.nx 	=*nvx;
	m_rt.nz 	=*nvz;
	m_rt.dx 	=*dvx;
	m_rt.dz     =*dvz;
	m_rt.x_src  =*sx_coord;
	m_rt.z_src  =*sz_coord;
	m_rt.x_min  =0.0;
	m_rt.x_max  = m_rt.dx * (float)(m_rt.nx - 1);
	m_rt.z_min  =0.0;
	m_rt.z_max  = m_rt.dz * (float)(m_rt.nz - 1);

	ds            =*dstep;
	n_table       = 101;
	nn= (NINT(m_rt.dx/MIN(m_rt.dx,m_rt.dz)) * m_rt.nx + NINT(m_rt.dz/MIN(m_rt.dx,m_rt.dz))*m_rt.nz) * 2;	
	area_grid     = m_rt.dx * m_rt.dz;
	n_point       =*n_p;
	TFlag         =*tflag;

	FileStatus = 3;
	locate     = 0;

//	printf("n_point=%d t_period=%f i_sigma_rule=%f\n",n_point,*t_period,*i_sigma_rule);fflush(stdout);
//	printf("dstep=%f nx=%d nz=%d dx=%f dz=%f x_src=%f z_src=%f x_max=%f z_max=%f\n",ds,m_rt.nx,m_rt.nz,m_rt.dx,m_rt.dz,m_rt.x_src,m_rt.z_src,m_rt.x_max,m_rt.z_max);fflush(stdout);

//	printf("nz=%d,nx=%d,nn=%d\n",m_rt.nz, m_rt.nx,nn);fflush(stdout);

	kernel_fresnel_norm = NULL;
	kernel_fresnel_norm = alloc2float(m_rt.nz, m_rt.nx);
	
	if (NULL == kernel_fresnel_norm)
	{
		printf ("Memory Allocation for 'kernel_fresnel_norm' Error.\n");
		return 1;
	}
	memset(&kernel_fresnel_norm[0][0],0,sizeof(float)*m_rt.nz*m_rt.nx);

	kernel_ray = NULL;
	kernel_ray = alloc2float (m_rt.nz, m_rt.nx);
	if (NULL == kernel_ray)
	{
		printf ("Memory Allocation for 'kernel_ray' Error.\n");
		return 1;
	}
	memset(&kernel_ray[0][0],0,sizeof(float)*m_rt.nz*m_rt.nx);

	r_rt.x_loc_ray = NULL;
	r_rt.x_loc_ray = alloc1float (nn);
	if (NULL == r_rt.x_loc_ray)
	{
		printf ("Memory Allocation for 'x_loc_ray' Error.\n");
		return 1;
	}
	memset(&r_rt.x_loc_ray[0],0,sizeof(float)*nn);

	r_rt.z_loc_ray = NULL;
	r_rt.z_loc_ray = alloc1float (nn);
	if (NULL == r_rt.z_loc_ray)
	{
		printf ("Memory Allocation for 'z_loc_ray' Error.\n");
		return 1;
	}
	memset(&r_rt.z_loc_ray[0],0,sizeof(float)*nn);

	r_rt.xp_ray = NULL;
	r_rt.xp_ray = alloc1float (nn);
	if (NULL == r_rt.xp_ray)
	{
		printf ("Memory Allocation for 'xp_ray' Error.\n");
		return 1;
	}
	memset(&r_rt.xp_ray[0],0,sizeof(float)*nn);

	r_rt.zp_ray = NULL;
	r_rt.zp_ray = alloc1float (nn);
	if (NULL == r_rt.zp_ray)
	{
		printf ("Memory Allocation for 'zp_ray' Error.\n");
		return 1;
	}
	memset(&r_rt.zp_ray[0],0,sizeof(float)*nn);

	r_rt.len_ray = NULL;
	r_rt.len_ray = alloc1float (nn);
	if (NULL == r_rt.len_ray)
	{
		printf ("Memory Allocation for 'len_ray' Error.\n");
		return 1;
	}
	memset(&r_rt.len_ray[0],0,sizeof(float)*nn);

	r_rt.tau_ray = NULL;
	r_rt.tau_ray = alloc1float (nn);
	if (NULL == r_rt.tau_ray)
	{
		printf ("Memory Allocation for 'tau_ray' Error.\n");
		return 1;
	}
	memset(&r_rt.tau_ray[0],0,sizeof(float)*nn);

	r_lamda.vel_s = NULL;
	r_lamda.vel_s = alloc1float (nn);
	if (NULL == r_lamda.vel_s)
	{
		printf ("Allocate Memory for 'r_lamda.vel_s' Error.\n");
		return 1;
	}
	memset(&r_lamda.vel_s[0],0,sizeof(float)*nn);

	r_lamda.vel_s_nn = NULL;
	r_lamda.vel_s_nn = alloc1float (nn);
	if (NULL == r_lamda.vel_s_nn)
	{
		printf ("Allocate Memory for 'r_lamda.vel_s_nn' Error.\n");
		return 1;
	}
	memset(&r_lamda.vel_s_nn[0],0,sizeof(float)*nn);

	q1 = NULL;
	q1 = alloc1float (nn);
	if (NULL == q1)
	{
		printf ("Allocate Memory for 'q1' Error.\n");
		return 1;
	}
	memset(&q1[0],0,sizeof(float)*nn);

	p1 = NULL;
	p1 = alloc1float (nn);
	if (NULL == p1)
	{
		printf ("Allocate Memory for 'p1' Error.\n");
		return 1;
	}
	memset(&p1[0],0,sizeof(float)*nn);

	q2 = NULL;
	q2 = alloc1float (nn);
	if (NULL == q2)
	{
		printf ("Allocate Memory for 'q2' Error.\n");
		return 1;
	}
	memset(&q2[0],0,sizeof(float)*nn);

	p2 = NULL;
	p2 = alloc1float (nn);
	if (NULL == p2)
	{
		printf ("Allocate Memory for 'p2' Error.\n");
		return 1;
	}
	memset(&p2[0],0,sizeof(float)*nn);

	deriv_2nd_src = NULL;
	deriv_2nd_src = alloc1float (nn);
	if (NULL == deriv_2nd_src)
	{
		printf ("Allocate Memory for 'deriv_2nd_src' Error.\n");
		return 1;
	}
	memset(&deriv_2nd_src[0],0,sizeof(float)*nn);

	deriv_2nd_rec = NULL;
	deriv_2nd_rec = alloc1float (nn);
	if (NULL == deriv_2nd_rec)
	{
		printf ("Allocate Memory for 'deriv_2nd_rec' Error.\n");
		return 1;
	}
	memset(&deriv_2nd_rec[0],0,sizeof(float)*nn);

	radius_fresnell_2 = NULL;
	radius_fresnell_2 = alloc1float (nn);
	if (NULL == radius_fresnell_2)
	{
		printf ("Allocate Memory for 'radius_fresnell_2' Error.\n");
		return 1;
	}
	memset(&radius_fresnell_2[0],0,sizeof(float)*nn);

	tab_v = NULL;
	tab_v = alloc2float (4, n_table);
	if (NULL == tab_v)
	{
		printf ("Memory Allocation Error.\n");
		return 1;
	}
	memset(&tab_v[0][0],0,sizeof(float)*n_table*4);

	tab_vx = NULL;
	tab_vx = alloc2float (4, n_table);
	if (NULL == tab_vx)
	{
		printf ("Memory Allocation Error.\n");
		return 1;
	}
	memset(&tab_vx[0][0],0,sizeof(float)*n_table*4);

	tab_vxx = NULL;
	tab_vxx = alloc2float (4, n_table);
	if (NULL == tab_vxx)
	{
		printf ("Memory Allocation Error.\n");
		return 1;
	}
	memset(&tab_vxx[0][0],0,sizeof(float)*n_table*4);
	
	vel = NULL;
	vel = alloc2float (m_rt.nz, m_rt.nx);
	if (NULL == vel)
	{
		printf ("Memory Allocation for 'vel' Error.\n");
		return 1;
	}
	memset(&vel[0][0],0,sizeof(float)*m_rt.nz*m_rt.nx);
	int igrid=0;
	for(ix=0; ix<m_rt.nx; ix++){
		for(iz=0; iz<m_rt.nz;iz++){
			vel[ix][iz]=v[igrid++];
		}
	}

	n_yn = extrpolation_coefficient_12O (tab_v, tab_vx, tab_vxx, n_table);

	if (0 != n_yn)
	{
		printf ("Calculating extrpolation coefficient error.\n");
		return 1;
	}

	for(i=0;i<n_point;i++){
		r_rt.x_loc_ray[i]  =x_loc_ray[i];
		r_rt.z_loc_ray[i]  =z_loc_ray[i];
		r_rt.xp_ray[i]     =xp_ray[i];
		r_rt.zp_ray[i]     =zp_ray[i];
		r_rt.tau_ray[i]    =tau_ray[i];
		// printf("%f %f %f %f %f\n",r_rt.x_loc_ray[i],r_rt.z_loc_ray[i],r_rt.xp_ray[i],r_rt.zp_ray[i],r_rt.tau_ray[i]);
	}

	//Fresnel, normalized
	n_yn = calculate_kernel_fresnel_norm (kernel_fresnel_norm, &m_rt,
			&r_rt, vel, q1, p1, q2, p2, deriv_2nd_src, deriv_2nd_rec, radius_fresnell_2,
			r_lamda, ds, n_point, *t_period, *i_sigma_rule, area_grid, tab_v, tab_vx,
			tab_vxx, n_table, elev, TFlag);
	if (0 != n_yn)
	{
		printf ("itr=%d,FUNCTION 'calculate_kernel_fresnel_norm' error.\n",*itr);
		free1float(q1);
		free1float(p1);
		free1float(q2);
		free1float(p2);
		free1float(deriv_2nd_src);
		free1float(deriv_2nd_rec);
		free1float(radius_fresnell_2);
		free1float(r_rt.x_loc_ray);
		free1float(r_rt.z_loc_ray);
		free1float(r_rt.xp_ray);
		free1float(r_rt.zp_ray);
		free1float(r_rt.len_ray);
		free1float(r_rt.tau_ray);
		free1float(r_lamda.vel_s);
		free1float(r_lamda.vel_s_nn);
		free2float(tab_v);
		free2float(tab_vx);
		free2float(tab_vxx);
		free2float(kernel_fresnel_norm);
		free2float(kernel_ray);
		free2float(vel);
		return 1;
	}


	int igrad=0;
	for(iz=0;iz<m_rt.nz;iz++){
		for(ix=0;ix<m_rt.nx;ix++){
			kernel[igrad]=kernel_fresnel_norm[ix][iz];
			igrad=igrad+1;
		}
	}

	free1float(q1);
	free1float(p1);
	free1float(q2);
	free1float(p2);
	free1float(deriv_2nd_src);
	free1float(deriv_2nd_rec);
	free1float(radius_fresnell_2);
	free1float(r_rt.x_loc_ray);
	free1float(r_rt.z_loc_ray);
	free1float(r_rt.xp_ray);
	free1float(r_rt.zp_ray);
	free1float(r_rt.len_ray);
	free1float(r_rt.tau_ray);
	free1float(r_lamda.vel_s);
	free1float(r_lamda.vel_s_nn);
	free2float(tab_v);
	free2float(tab_vx);
	free2float(tab_vxx);
	free2float(kernel_fresnel_norm);
	free2float(kernel_ray);
	free2float(vel);
	return 0;
}
