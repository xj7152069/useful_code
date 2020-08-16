
#ifndef PI
#define PI 3.14159265358979
#endif
/*---------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------*/
/*  Function 'get_kernel_gauss_unit' calculate kernel of a point.                        */
/*                                                                                       */
/*  Input:                                                                               */
/*    float radius_fresnell_2: Radius' square of fresnell zone which contain the point   */
/*      holding calculated;                                                              */
/*    float distance2_from_ray: Square distance from current point to ray;               */
/*    float i_sigma_rule: Square deviation equals pow (r, 2) / pow (i_sigma_rule, 2),    */
/*      where 'r' is radius of fresnell zone.                                            */
/*                                                                                       */
/*  Output:                                                                              */
/*    float *kernel: Kernel value calculated.                                            */
int get_kernel_gauss_unit (float *kernel, float radius_fresnell_2, float distance2_from_ray,
	float i_sigma_rule)
{
	float sigma_2;

	/*if (0 == radius_fresnell_2)
	{
		//(*kernel) = 1.;
		*kernel = 0.1;
		return 0;
	}

	sigma_2 = radius_fresnell_2 / pow(i_sigma_rule, 2.);
	*kernel = exp (-distance2_from_ray / 2. / sigma_2) / sqrt (2. * PI * sigma_2);*/

	if (0 == radius_fresnell_2)
	{
		*kernel = 1.;
		return 0;
	}

	sigma_2 = radius_fresnell_2 / pow(i_sigma_rule, 2.);
	//*kernel = exp (-distance2_from_ray / 2. / sigma_2) / sqrt (2. * PI * sigma_2);
	*kernel = exp (-distance2_from_ray / 2. / sigma_2);

	return 0;
}

/*---------------------------------------------------------------------------------------*/
/*  Function 'get_kernel_gauss_norm' calculate kernel of a point.                        */
/*                                                                                       */
/*  Input:                                                                               */
/*    float radius_fresnell_2: Radius' square of fresnell zone which contain the point   */
/*      holding calculated;                                                              */
/*    float distance2_from_ray: Square distance from current point to ray;               */
/*    float i_sigma_rule: Square deviation equals pow (r, 2) / pow (i_sigma_rule, 2),      */
/*      where 'r' is radius of fresnell zone.                                            */
/*                                                                                       */
/*  Output:                                                                              */
/*    float *kernel: Kernel value calculated.                                            */
int get_kernel_gauss_norm (float *kernel, float radius_fresnell_2, float distance2_from_ray,
	float i_sigma_rule, float dz)
{
	float sigma_2;

	/*if (0 == radius_fresnell_2)
	{
		//(*kernel) = 1.;
		*kernel = 0.1;
		return 0;
	}

	sigma_2 = radius_fresnell_2 / pow(i_sigma_rule, 2.);
	*kernel = exp (-distance2_from_ray / 2. / sigma_2) / sqrt (2. * PI * sigma_2);*/

	if (dz > radius_fresnell_2)
	{
		//(*kernel) = 0.1;
		(*kernel) = 0;
		return 0;
	}

	sigma_2 = radius_fresnell_2 / pow(i_sigma_rule, 2.);
	*kernel = exp (-distance2_from_ray / 2. / sigma_2) / sqrt (2. * PI * sigma_2);
	//*kernel = exp (-distance2_from_ray / 2. / sigma_2);

	return 0;
}

/*---------------------------------------------------------------------------------------*/
/*  Function 'get_radius_fresnell_2' calculate radius squares of fresnell zone along a   */
/*    ray.                                                                               */
/*                                                                                       */
/*  Input:                                                                               */
/*    float *deriv_2nd_src: 2ed devirative of travel time from source to perpendicular foot of    */
/*      current point;                                                                   */
/*    float *deriv_2nd_rec: 2ed devirative of travel time from receiver to perpendicular foot of  */
/*      current point;                                                                   */
/*    float t_period: Wave period for certain frequency;                                 */
/*    int n_point: Number of point calculated on a ray.                                  */
/*                                                                                       */
/*  Output:                                                                              */
/*    float radius_fresnell_2: Radius' square of fresnell zone along the ray.            */
int get_radius_fresnell_2 (float *deriv_2nd_src, float *deriv_2nd_rec, float t_period,
	float *radius_fresnell_2, int n_point)
{
	int i;
	float ftmp1, ftmp2;

	ftmp1 = pow (t_period, 2.);

	for (i=0; i<n_point; i++)
	{
		if (1e19 < deriv_2nd_src[i] || 1e19 < deriv_2nd_rec[i])
		{
			radius_fresnell_2[i] = 0;
			continue;
		}
		ftmp2 = pow (deriv_2nd_src[i] - deriv_2nd_rec[i], 2.);
		if (1e-19 > ftmp2)
		{
			printf ("i = %d, deriv_2nd_src = %f, deriv_2nd_rec = %f\n", i, deriv_2nd_src[i], deriv_2nd_rec[i]);
			printf ("Second derivative from source equals from receiver.\n");
			return 1;
		}
		radius_fresnell_2[i] = sqrt (ftmp1 / ftmp2);

		// printf ("i = %d, deriv_2nd_src = %f, deriv_2nd_rec = %f\n", i, deriv_2nd_src[i], deriv_2nd_rec[i]);
		// printf("%d %f\n",i,sqrt(radius_fresnell_2[i]));

	}
	return 0;
}

/*---------------------------------------------------------------------------------------*/
/*  Function 'get_deri_2nd_src_rec' calculate 2ed derivative of travel time field on a ray, both    */
/*    from source and receiver.                                                          */
/*                                                                                       */
/*  Input:                                                                               */
/*    float *q1, float *p1, float *q2, float *p2: Propagating matrix;                    */
/*    int n_point: Point number calculated on a ray.                                     */
/*                                                                                       */
/*  Output:                                                                              */
/*    float *deriv_2nd_src: Second derivative of travel time field on a ray from source to        */
/*      receiver.                                                                        */
/*    float *deriv_2nd_rec: Second derivative of travel time field on a ray from receiver to      */
/*      source.                                                                          */
int get_deri_2nd_src_rec (float *q1, float *p1, float *q2, float *p2, float *deriv_2nd_src, float *deriv_2nd_rec,
	int n_point)
{
	int i;
	float ftmp1, ftmp2;

	for (i=0; i<n_point; i++)
	{
		/*if (q2[i] < 1e-20 && q2[i] > 0)
			deriv_2nd_src[i] = 1e20;
		else if (q2[i] > -1e-20 && q2[i] < 0)
			deriv_2nd_src[i] = -1e20;*/
		if (fabs (q2[i]) < 1e-20)
			deriv_2nd_src[i] = 1e20;
		else
			deriv_2nd_src[i] = p2[i] / q2[i];

//		printf("i=%d p2=%f q2=%f src=%f\n",i,p2[i],q2[i],deriv_2nd_src[i]);
		// printf("i=%d q1=%f q2=%f\n",i,q1[i],q2[i]);

		ftmp1 = - q1[i] * q2[n_point-1] + q2[i] * q1[n_point-1];
		if (fabs (ftmp1) < 1e-20)
			deriv_2nd_rec[i] = 1e20;
		else
		{
			ftmp2 = - p1[i] * q2[n_point-1] + p2[i] * q1[n_point-1];
			deriv_2nd_rec[i] = ftmp2 / ftmp1;
		}
		// printf("i=%d ftmp2=%f ftmp1=%f rec=%f\n",i,ftmp2,ftmp1,deriv_2nd_rec[i]);
	}

//	printf("<=============================================================================>\n");
//	printf("<================================src_rec over=================================>\n");
//	printf("<=============================================================================>\n");

	return 0;
}

/*---------------------------------------------------------------------------------------*/
/*  Function 'get_perpendicular_foot' used to find the perpendicular foot of a point in  */
/*    a ray.                                                                             */
/*                                                                                       */
/*  Input:                                                                               */
/*    float xm, float zm: Location of the point which we want to find perpendicular;     */
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
/*    int n_point: Number of point calculated on a ray.                                  */
/*                                                                                       */
/*  Output:                                                                              */
/*    float *xn, float *zn: Location of perpendicular foot;                              */
/*    int *in: Perpendicular lacated between the (in)th and the (in+1)th point of the    */
/*      ray.                                                                             */
int get_perpendicular_foot (float xm, float zm, float *xn, float *zn, int *in, model_rt *m_rt,
	ray_rt *r_rt, int n_point)
{
	int i, i0;
	float xa, za, xb, zb, dx_ray, dz_ray, ftmp1, ftmp2;

	i0 = 0;

	/*-------------------------------------------------------------------------------*/
	/*  First, we should find the closest point on the ray to the current point.     */
	ftmp1 = pow (xm-r_rt->x_loc_ray[0], 2.) + pow (zm-r_rt->z_loc_ray[0], 2.);
	for (i=1; i<n_point; i++)
	{
		ftmp2 = pow (xm-r_rt->x_loc_ray[i], 2.)
			+ pow (zm-r_rt->z_loc_ray[i], 2.);
		if (ftmp2 < ftmp1)
		{
			i0 = i;
			ftmp1 = ftmp2;
		}
	}

	//n_yn = 0;

	/*-------------------------------------------------------------------------------*/
	/*  If the closest point is token by 0, this means the first point of the ray.   */
	if (0 == i0)
	{
		xa = r_rt->x_loc_ray[0];
		za = r_rt->z_loc_ray[0];
		xb = r_rt->x_loc_ray[1];
		zb = r_rt->z_loc_ray[1];
		dx_ray = xb - xa;
		dz_ray = zb - za;

		/*  If line AB perpendicular with X coordinate.  */
		if (fabs(dx_ray) < 1e-10)
		{
			*xn = xa;
			if ((*xn >= xa && *xn <= xb) || (*xn <= xa && *xn >=xb))
			{
				*zn = zm;
				//n_yn = 1;
				*in = 0;
				return 0;
			}
		}
		else
		{
			*xn = (dx_ray * dx_ray * xm + dz_ray * dz_ray * xa
				+ dx_ray * dz_ray * (zm - za))
				/ (dx_ray * dx_ray + dz_ray * dz_ray);
			/*  If this 'if' do not right, means current point do not have  */
			/*  perpendicular foot.                                         */
			if ((*xn >= xa && *xn <= xb) || (*xn <= xa && *xn >= xb))
			{
				*zn = dz_ray * (*xn - xa) / dx_ray + za;
				//n_yn = 1;
				*in = 0;
				return 0;
			}
		}
	}
	/*  If the closest point is token by n_point-1, this means the last point of the ray.   */
	else if (n_point-1 == i0)
	{
		xa = r_rt->x_loc_ray[n_point-2];
		za = r_rt->z_loc_ray[n_point-2];
		xb = r_rt->x_loc_ray[n_point-1];
		zb = r_rt->z_loc_ray[n_point-1];
		dx_ray = xb - xa;
		dz_ray = zb - za;

		/*  If line AB perpendicular with X coordinate.  */
		if (fabs(dx_ray) < 1e-10)
		{
			*xn = xa;
			if ((*xn >= xa && *xn <= xb) || (*xn <= xa && *xn >=xb))
			{
				*zn = zm;
				//n_yn = 1;
				*in = n_point-2;
				return 0;
			}
		}
		else
		{
			*xn = (dx_ray * dx_ray * xm + dz_ray * dz_ray * xa
				+ dx_ray * dz_ray * (zm - za))
				/ (dx_ray * dx_ray + dz_ray * dz_ray);
			/*  If this 'if' do not right, means current point do not have  */
			/*  perpendicular foot.                                         */
			if ((*xn >= xa && *xn <= xb) || (*xn <= xa && *xn >= xb))
			{
				*zn = dz_ray * (*xn - xa) / dx_ray + za;
				//n_yn = 1;
				*in = n_point-2;
				return 0;
			}
		}
	}
	/*  Closest point is a common point on a ray.   */
	else
	{
		/*  If line AB perpendicular with X coordinate.  */
		for (i=i0-1; i<i0+1; i++)
		{
			xa = r_rt->x_loc_ray[i];
			za = r_rt->z_loc_ray[i];
			xb = r_rt->x_loc_ray[i+1];
			zb = r_rt->z_loc_ray[i+1];
			dx_ray = xb - xa;
			dz_ray = zb - za;

			if (fabs(dx_ray) < 1e-10)
			{
				*xn = xa;
				if ((*xn >= xa && *xn <= xb) || (*xn <= xa && *xn >=xb))
				{
					*zn = zm;
					//n_yn = 1;
					*in = i;
					return 0;
				}
			}
			else
			{
				*xn = (dx_ray * dx_ray * xm + dz_ray * dz_ray * xa
					+ dx_ray * dz_ray * (zm - za))
					/ (dx_ray * dx_ray + dz_ray * dz_ray);
				if ((*xn >= xa && *xn <= xb) || (*xn <= xa && *xn >= xb))
				{
					*zn = dz_ray * (*xn - xa) / dx_ray + za;
					//n_yn = 1;
					*in = i;
					return 0;
				}
			}

			/*  If we cannot find perpendicular foot of current point to broken-line  */
			/*  ray, the point 'A' is just the perpendicular foot.                    */
			if (i == i0)
			{
				*xn = xa;
				*zn = za;
				//n_yn = 1;
				*in = i;
				return 0;
			}
		}
	}

	return 1;
}

int get_value_along_ray(ray_rt *r_rt, model_rt *m_rt, ray_lamda *r_lamda, float **vel,
	int n_point, float **tab_v, float **tab_vx, float **tab_vxx, int n_table)
{
	int k, n_yn;
	float lamda_xx_s, lamda_xz_s, lamda_zz_s;
	float h_s_x, h_s_z, h_n_x, h_n_z, ftmp1;

	for (k=0; k<n_point; k++)
	{

		n_yn = get_vel_differential_vel_2O (&r_lamda->vel_s[k], &lamda_xx_s,
			&lamda_xz_s, &lamda_zz_s, r_rt->x_loc_ray[k], r_rt->z_loc_ray[k], m_rt->nx, m_rt->nz, m_rt->dx, m_rt->dz,
			m_rt->x_min, m_rt->z_min, m_rt->x_max, m_rt->z_max, vel, tab_v, tab_vx, tab_vxx, n_table);
		if (0 != n_yn)
		{
			printf ("Calculating vel_differential_xz error.\n");
			return 1;
		}


		ftmp1 = sqrt (r_rt->xp_ray[k] * r_rt->xp_ray[k] + r_rt->zp_ray[k] * r_rt->zp_ray[k]);
	/*	luofei modify	*/
		// h_s_x = r_rt->xp_ray[k] / ftmp1;
		// h_s_z = r_rt->zp_ray[k] / ftmp1;
		h_s_x = r_rt->xp_ray[k] / (ftmp1+0.000001);
		h_s_z = r_rt->zp_ray[k] / (ftmp1+0.000001);
	/********************/
		h_n_x = h_s_z;
		h_n_z = -h_s_x;
		ftmp1 = h_n_x * h_s_z - h_s_x * h_n_z;
		if (0 > ftmp1)
		{
			h_n_x = -h_s_z;
			h_n_z = h_s_x;
		}

		/*r_lamda->vel_s_nn[k] = r_rt->zp_ray[k] * r_rt->zp_ray[k] * lamda_xx_s
			- 2.0 * r_rt->xp_ray[k] * r_rt->zp_ray[k] * lamda_xz_s
			+ r_rt->xp_ray[k] * r_rt->xp_ray[k] * lamda_zz_s;*/
		r_lamda->vel_s_nn[k] = h_n_x * h_n_x * lamda_xx_s + 2. * h_n_x * h_n_z * lamda_xz_s
			+ h_n_z * h_n_z * lamda_zz_s;

		// printf("i=%d x=%f z=%f vel_s=%f vel_nn=%f\n ",k,r_rt->x_loc_ray[k],r_rt->z_loc_ray[k],r_lamda->vel_s[k],r_lamda->vel_s_nn[k]);
	}

//	printf("<==================================================================================>\n");
//	printf("<============================= tabxx over   =======================================>\n");
//	printf("<==================================================================================>\n");


	return 0;
}

int get_radius_approximate_fresnell_2(ray_rt *r_rt, model_rt *m_rt, ray_lamda *r_lamda, float **vel,
	int n_point, float **tab_v, float **tab_vx, float **tab_vxx, int n_table ,float *radius_fresnell_2,float t_period)
{
	int k, n_yn;
	float lamda_xx_s, lamda_xz_s, lamda_zz_s, Total_L, Total_X;
	
	Total_L = 0.0; 
	Total_X = 0.0;

	for (k=0; k<n_point-1; k++)
	{
		Total_L = Total_L+ sqrtf((r_rt->x_loc_ray[k+1]-r_rt->x_loc_ray[k])*(r_rt->x_loc_ray[k+1]-r_rt->x_loc_ray[k])+
								 (r_rt->z_loc_ray[k+1]-r_rt->z_loc_ray[k])*(r_rt->z_loc_ray[k+1]-r_rt->z_loc_ray[k]));
	}
	
	
	for (k=1; k<n_point; k++)
	{

		n_yn = get_vel_differential_vel_2O (&r_lamda->vel_s[k], &lamda_xx_s,
			&lamda_xz_s, &lamda_zz_s, r_rt->x_loc_ray[k], r_rt->z_loc_ray[k], m_rt->nx, m_rt->nz, m_rt->dx, m_rt->dz,
			m_rt->x_min, m_rt->z_min, m_rt->x_max, m_rt->z_max, vel, tab_v, tab_vx, tab_vxx, n_table);
		if (0 != n_yn)
		{
			printf ("Calculating vel_differential_xz error.\n");
			return 1;
		}

		Total_X = Total_X + sqrtf((r_rt->x_loc_ray[k-1]-r_rt->x_loc_ray[k])*(r_rt->x_loc_ray[k-1]-r_rt->x_loc_ray[k])+
								 (r_rt->z_loc_ray[k-1]-r_rt->z_loc_ray[k])*(r_rt->z_loc_ray[k-1]-r_rt->z_loc_ray[k])); 

		radius_fresnell_2[k] = (t_period * r_lamda->vel_s[k] * Total_X * (Total_L-Total_X))/Total_L;

		
	}


	// for(k=0;k<n_point;k++){
	// 	printf("%d %f\n",k,sqrt(radius_fresnell_2[k]));
	// }

	return 0;
}


int get_qp (float *q, float *p, ray_lamda *r_lamda, int n_point, float ds)
{
	int i;
	float kq1, kp1, kq2, kp2, kq3, kp3, kq4, kp4;
	float vel_tmp, vel_s_nn_tmp;


/*-----------------------------RK Method------------------------------*/
	for (i=0; i<n_point-1; i++)
	{
		// printf("i=%d vel_s=%f vel_nn=%f vel_s1=%f vel_nn1=%f\n ",i,r_lamda->vel_s[i],r_lamda->vel_s_nn[i],r_lamda->vel_s[i+1],r_lamda->vel_s_nn[i+1]);
		kq1 = r_lamda->vel_s[i] * p[i];
		kp1 = (-1.0) / (r_lamda->vel_s[i] * r_lamda->vel_s[i]) * r_lamda->vel_s_nn[i] * q[i];
		vel_tmp = (r_lamda->vel_s[i] + r_lamda->vel_s[i+1]) / 2.0;
		vel_s_nn_tmp = (r_lamda->vel_s_nn[i] + r_lamda->vel_s_nn[i+1]) / 2.0;
		kq2 = vel_tmp * (p[i] + kp1 * ds / 2.);
		kp2 = (-1.0) / (vel_tmp * vel_tmp) * vel_s_nn_tmp *
			(q[i] + kq1 * ds / 2.);
		kq3 = vel_tmp * (p[i] + kp2 * ds / 2.);
		kp3 = (-1.0) / (vel_tmp * vel_tmp) * vel_s_nn_tmp *
			(q[i] + kq2 * ds / 2.);
		kq4 = r_lamda->vel_s[i+1] * (p[i] + kp3 * ds);
		kp4 = (-1.0) / (r_lamda->vel_s[i+1] * r_lamda->vel_s[i+1]) * r_lamda->vel_s_nn[i+1] *
			(q[i] + kq3 * ds);
		
		q[i+1] = q[i] + ds / 6.0 *
			(kq1 + kq2 * 2.0 + kq3 * 2.0 + kq4);
		p[i+1] = p[i] + ds / 6.0 *
			(kp1 + kp2 * 2.0 + kp3 * 2.0 + kp4);
		// printf("i=%d q=%f p=%f\n",i,q[i+1],p[i+1]);
	}
	return 0;
}
