
/*---------------------------------------------------------------------------------------*/

#ifndef PI
#define PI 3.14159265358979
#endif

int get_qp (float *q, float *p, ray_lamda *r_lamda, int n_point, float ds);
int get_value_along_ray(ray_rt *r_rt, model_rt *m_rt, ray_lamda *r_lamda, float **vel,
	int n_point, float **tab_v, float **tab_vx, float **tab_vxx, int n_table);
int get_perpendicular_foot (float xm, float zm, float *xn, float *zn, int *in, model_rt *m_rt,
	ray_rt *r_rt, int n_point);
int get_deri_2nd_src_rec (float *q1, float *p1, float *q2, float *p2, float *deriv_2nd_src, float *deriv_2nd_rec,
	int n_point);
int get_radius_fresnell_2 (float *deriv_2nd_src, float *deriv_2nd_rec, float t_period,
	float *radius_fresnell_2, int n_point);
int get_kernel_gauss_norm (float *kernel, float radius_fresnell_2, float distance2_from_ray,
	float i_sigma_rule, float dz);
int get_radius_approximate_fresnell_2(ray_rt *r_rt, model_rt *m_rt, ray_lamda *r_lamda, float **vel,
	int n_point, float **tab_v, float **tab_vx, float **tab_vxx, int n_table ,float *radius_fresnell_2,float t_period);

/*---------------------------------------------------------------------------------------*/
/*  Function 'calculate_fresnell_kernel' calculate kernel funtion for tomography.        */
/*  Input:                                                                               */
/*    model_rt *m_rt:                                                                    */
/*      nx: Model size in X direction;                                                   */
/*      nz: Model size in Z direction;                                                   */
/*      dx: Sample interval in X direction;                                              */
/*      dz: Sample interval in Z direction;                                              */
/*      x_src: Source location in X direction;                                           */
/*      z_src: Source location in Z direction;                                           */
/*      x_rec: Receiver location in X direction;                                         */
/*      z_rec: Receiver location in Z direction;                                         */
/*      thita: Ray initial direction(angle);                                             */
/*                                                                                       */
/*    ray_rt *r_rt:                                                                      */
/*      x_loc_ray: Ray location in X direction;                                          */
/*      z_loc_ray: Ray location in Z direction;                                          */
/*      xp_ray: Ray direction vector in X direction;                                     */
/*      zp_ray: Ray direction vector in Z direction;                                     */
/*      len_ray: Ray length at each point of ray;                                        */
/*      tau_ray: Travel time at each point of ray;                                       */
/*                                                                                       */
/*    float **vel: Velocity field;                                                       */
/*    int n_bdr: Buffer boundary size in each direction;                                 */
/*    float ds: Interval of 2 adjacent of ray;                                           */
/*    int n_point: Point number calculated on the ray;                                   */
/*    float t_period: Period of a given frequency;                                       */
/*    float i_sigma_rule: Square deviation equals pow (r, 2) / pow (i_sigma_rule, 2),    */
/*      where 'r' is radius of fresnell zone.                                            */
/*                                                                                       */
/*  Output:                                                                              */
/*    float **kernel: Kernel function for tomography.                                    */

/*
int calculate_fresnell_kernel_gauss (float **kernel, model_rt *m_rt, ray_rt *r_rt,
	float **vel, int n_bdr, float ds, int n_point, float t_period, float i_sigma_rule)
*/
int calculate_kernel_fresnel_norm (float **kernel, model_rt *m_rt, ray_rt *r_rt,
	float **vel, float *q1, float *p1, float *q2, float *p2, float *deriv_2nd_src,
	float *deriv_2nd_rec, float *radius_fresnell_2, ray_lamda r_lamda, float ds,
	int n_point, float t_period, float i_sigma_rule, float area_grid,
	float **tab_v, float **tab_vx, float **tab_vxx, int n_table, float *elev, int TFlag)
{
	int i, j, in, n_yn,jz_b;
	//float *q1, *p1, *q2, *p2, *deriv_2nd_src, *deriv_2nd_rec, *radius_fresnell_2,
	float distance2_from_ray, radius_fresnell_2_tmp,
	      xm, zm, xn, zn;


	switch(TFlag)
	{
		case 1:
		{
			/*---------------------------------------------------------------------------*/
			/*  Identity matrix as initial value, calculating propagate matrix.          */
			q1[0] = 1.;
			p1[0] = 0;
			q2[0] = 0;
			p2[0] = 1.;
	
			/*---------------------------------------------------------------------------*/
			/*  Calculating derivative of velocity field on a ray.                       */
			n_yn = get_value_along_ray (r_rt, m_rt, &r_lamda, vel, n_point, tab_v, tab_vx, tab_vxx, n_table);
			if (0 != n_yn)
			{
				printf ("Calculating velocity derivative along ray Error.\n");
				return 1;
			}

			/*---------------------------------------------------------------------------*/
			/*  Calculating propagating matrix.                                          */
			n_yn = get_qp (q1, p1, &r_lamda, n_point, ds);
			if (0 != n_yn)
			{
				printf ("Calculating 'qp' along ray Error.\n");
				return 1;
			}

			n_yn = get_qp (q2, p2, &r_lamda, n_point, ds);
			if (0 != n_yn)
			{
				printf ("Calculating 'qp' along ray Error.\n");
				return 1;
			}

			/*---------------------------------------------------------------------------*/
			/*  Calculating 2ed derivative of travel time field on a ray.                */
			n_yn = get_deri_2nd_src_rec (q1, p1, q2, p2, deriv_2nd_src, deriv_2nd_rec, n_point);
			if (0 != n_yn)
			{
				printf ("Calculating 'deriv_2nd_src' or 'deriv_2nd_rec' along ray Error.\n");
				return 1;
			}

			// /*-------------------------------------------------------------------------------------------------*/
			/*  Calculating radius' square of fresnell zone. 计算菲涅尔体的半径H(s)----------------------------*/
			n_yn = get_radius_fresnell_2 (deriv_2nd_src, deriv_2nd_rec, t_period, radius_fresnell_2, n_point);
			if (0 != n_yn)
			{
				printf ("Calculating square of fresnell radius Error.\n");
				return 1;
			}
			break;
		}
		case 2:
		{
			/*  Calculating radius' square of approximate fresnell zone. 计算近似菲涅尔体的半径H(s)----------------------------*/
			n_yn = get_radius_approximate_fresnell_2(r_rt, m_rt, &r_lamda, vel, n_point, tab_v, tab_vx, tab_vxx, n_table,radius_fresnell_2,t_period);
			if (0 != n_yn)
			{
				printf ("Calculating square of approximate fresnell radius Error.\n");
				return 1;
			}
			break;
		}
	}
	
	for (i=0; i<m_rt->nx; i++)
	//for (i=n_bdr; i<m_rt->nx-n_bdr; i++)
	{
		xm  = m_rt->dx * (float)i;
		jz_b= round(elev[i]/m_rt->dz);
		//for (j=0; j<m_rt->nz; j++)
		for (j=jz_b; j<m_rt->nz; j++)
		//for (j=n_bdr; j<m_rt->nz-n_bdr; j++)
		{
			zm = m_rt->dz * (float)j;

			/*----------------------------------------------------------------------*/
			/*  Find perpendicular foot in ray of point (xm, zm)                    */
			n_yn = get_perpendicular_foot (xm, zm, &xn, &zn, &in, m_rt,
				r_rt, n_point);
			if (0 != n_yn && 1 != n_yn)
			{
				printf ("Calculating Perpendicular foot Error. i = %d, j = %d\n",
					i, j);
				return 1;
			}

			/*----------------------------------------------------------------------*/
			/*  n_yn equals 1 means we cannot find perpendicular foot of current    */
			/*  point. Surely kernel equals 0, here we let kernel that equals       */
			/*  be a negtive number.                                                */
			if (1 == n_yn)
			{
				//kernel[i][j] = -1;
				kernel[i][j] = 0;
				continue;
			}

			/*----------------------------------------------------------------------*/
			/*  Calculate square of distance from current point to ray.             */
			distance2_from_ray = pow (xm-xn, 2.) + pow (zm-zn, 2.);

			/*----------------------------------------------------------------------*/
			/*            计算射线路径上垂足处(xn,zn)的菲涅尔半径,                  */
			/*       之前计算的是射线路径上采样点in或in+1...处的菲涅尔半径          */
			/*       Calculate radius of fresnell zone contain current point.       */
			radius_fresnell_2_tmp = radius_fresnell_2[in]
				+ (radius_fresnell_2[in+1] - radius_fresnell_2[in])
				* sqrt (pow (xn-r_rt->x_loc_ray[in], 2.) + pow (zn-r_rt->z_loc_ray[in], 2.))
				/ ds;
			/*----------------------------------------------------------------------*/

			/*----------------------------------------------------------------------*/
			if (distance2_from_ray > radius_fresnell_2_tmp)
			{
				//kernel[i][j] = -1;
				kernel[i][j] = 0;
				continue;
			}
			/*----------------------------------------------------------------------*/

			/*----------------------------------------------------------------------*/
			/*  Calculate kernel of current point. In one fresnell zone, kernel     */
			/*  function is Gaussian distribution.                                  */
			n_yn = get_kernel_gauss_norm (&kernel[i][j], radius_fresnell_2_tmp,
				distance2_from_ray, i_sigma_rule, m_rt->nz);
			if (0 != n_yn)
			{
				printf ("Calculate kernel using Gaussian error.\n");
				return 1;
			}

			kernel[i][j] *= area_grid;

		}
	}

	return 0;
}
