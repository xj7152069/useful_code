/*-----------------------------------------------------------------------------------*/
/*                   Main Function: "ray_tracing_hamilton()"                         */
/*  For: Initial Parameters Ray Tracing                                              */
/*                                                                                   */
/*  return n_point ---- Number of Points Calculated On a Ray                         */
/*                                                                                   */
/*  InPut:                                                                           */
/*    Struct model_rt *m_rt;                                                         */
/*      int nx, nz ---- Size of Model in X and Z Coordinates;                        */
/*      float dx, dz ---- Sample Interval in X and Z Coordinates;                    */
/*      float x_src, z_src ---- Source Location;                                     */
/*      float thita ---- Angle of Initial Direction.                                 */
/*                                                                                   */
/*    float **vel ---- Velocity Field;                                               */
/*    int n_bdr ---- Expand Boundary Number in Each Direction;                       */
/*    float ds ---- Ray Length Interval for Numerical Calculating;                   */
/*                                                                                   */
/*  OutPut:                                                                          */
/*    Struct ray_rt *r_rt;                                                           */
/*      double *x_loc_ray, *z_loc_ray ---- Location of n_point Ray Points;           */
/*      double *xp_ray, *zp_ray ---- Ray Parameters of n_point Ray Points;           */
/*      double *len_ray ---- Ray Length From Source to n_point Ray Points;           */
/*      double *tau_ray ---- Travel Time From Source to n_point Ray Points.          */
/*                                                                                   */
/*    int *i_edge ---- Token Which Boundary Ray Hit.                                 */
/*                                                                                   */
/*-----------------------------------------------------------------------------------*/

/*typedef struct
{
	int nx;
	int nz;
	float dx;
	float dz;
	float x_src;
	float z_src;
	float thita;
}model_rt;

typedef struct
{
	double *x_loc_ray;
	double *z_loc_ray;
	double *xp_ray;
	double *zp_ray;
	double *len_ray;
	double *tau_ray;
}ray_rt;*/

int cell_ray_tracing (model_rt *m_rt, ray_rt *r_rt, float **vel, float ds, int i, float **tab_v, float **tab_vx, int n_table);

int ray_tracing_hamilton (model_rt *m_rt, ray_rt *r_rt, float **vel, float ds, int *i_edge,
	float t_slice, float *l_ray, float **tab_v, float **tab_vx, int n_table, int *n_point)
{
	int i, n_yn;
	double v0, x_min, x_max, z_min, z_max;
	float ftmp1;

	x_min = 0;
	x_max = (double)m_rt->dx * (double)(m_rt->nx - 1);
	z_min = 0;
	z_max = (double)m_rt->dz * (double)(m_rt->nz - 1);
	/*
	x_min = (double)m_rt->dx * (double)n_bdr;
	x_max = (double)m_rt->dx * (double)(m_rt->nx - 1) - x_min;
	z_min = (double)m_rt->dz * (double)n_bdr;
	z_max = (double)m_rt->dz * (double)(m_rt->nz - 1) - z_min;
	*/

	n_yn = get_vel (&ftmp1, m_rt->x_src, m_rt->z_src, m_rt->nx, m_rt->nz, m_rt->dx, m_rt->dz, m_rt->x_min,
		m_rt->z_min, m_rt->x_max, m_rt->z_max, vel, tab_v, tab_vx, n_table);
	if (0 != n_yn)
	{
		printf ("Calculating velocity of source location error.\n");
		return 1;
	}
	v0 = (double)ftmp1;

	r_rt->xp_ray[0] = cos ((double)m_rt->thita) / v0;
	r_rt->zp_ray[0] = sin ((double)m_rt->thita) / v0;
	r_rt->x_loc_ray[0] = (double)m_rt->x_src;
	r_rt->z_loc_ray[0] = (double)m_rt->z_src;
	r_rt->tau_ray[0] = 0;
	r_rt->len_ray[0] = 0;

	//n_point = 0;
	i = 0;
	while (1)
	{
		if (r_rt->x_loc_ray[i] < x_min)
		{
			*i_edge = -1;
			break;
		}
		else if (r_rt->x_loc_ray[i] > x_max)
		{
			*i_edge = 1;
			break;
		}
		else if (r_rt->z_loc_ray[i] < z_min)
		{
			*i_edge = -2;
			break;
		}
		else if (r_rt->z_loc_ray[i] > z_max)
		{
			*i_edge = 2;
			break;
		}

		cell_ray_tracing (m_rt, r_rt, vel, ds, i, tab_v, tab_vx, n_table);

		/*
		if (t_slice <= r_rt->tau_ray[i+1] && t_slice > r_rt->tau_ray[i])
		{
			*l_ray = (float)i + (t_slice - r_rt->tau_ray[i]) / (r_rt->tau_ray[i+1] - r_rt->tau_ray[i]);
		}
		*/
		
		i ++;
	}

	*n_point = i;	//The Last Point Out of Model Already.

	return 0;

}

int cell_ray_tracing (model_rt *m_rt, ray_rt *r_rt, float **vel, float ds, int i, float **tab_v, float **tab_vx, int n_table)
{
	int n_yn;
	double h, v0, k1[4], k2[4], k3[4], k4[4], lamda_x, lamda_z, dtmp1, dtmp2;
	float ftmp1, ftmp2, ftmp3;

	h = (double)ds;

	/*ix = (int)(r_rt->x_loc_ray[i] / m_rt->dx);
	dx0 = (double)(r_rt->x_loc_ray[i] - m_rt->dx * (float)ix);
	iz = (int)(r_rt->z_loc_ray[i] / m_rt->dz);
	dz0 = (double)(r_rt->z_loc_ray[i] - m_rt->dz * (float)iz);

	v0 = (double)vel[ix][iz]
		+ (double)(vel[ix+1][iz] - vel[ix][iz]) * dx0 / (double)m_rt->dx
		+ (double)(vel[ix][iz+1] - vel[ix][iz]) * dz0 / (double)m_rt->dz;

	dtmp1 = (double)(vel[ix+1][iz] - vel[ix-1][iz]) / (double)m_rt->dx / 2.;
	dtmp2 = (double)(vel[ix+2][iz] - vel[ix][iz]) / (double)m_rt->dx / 2.;
	lamda_x = dtmp1 + (dtmp2 - dtmp1) * dx0 / (double)m_rt->dx;

	dtmp1 = (double)(vel[ix][iz+1] - vel[ix][iz-1]) / (double)m_rt->dz / 2.;
	dtmp2 = (double)(vel[ix][iz+2] - vel[ix][iz]) / (double)m_rt->dz / 2.;
	lamda_z = dtmp1 + (dtmp2 - dtmp1) * dz0 / (double)m_rt->dz;*/

	n_yn = get_vel_differential_vel (&ftmp1, &ftmp2, &ftmp3, r_rt->x_loc_ray[i], r_rt->z_loc_ray[i],
		m_rt->nx, m_rt->nz, m_rt->dx, m_rt->dz, m_rt->x_min, m_rt->z_min, m_rt->x_max, m_rt->z_max,
		vel, tab_v, tab_vx, n_table);
	if (0 != n_yn)
	{
		printf ("Calculating vel_differential_vel error.\n");
		return 1;
	}
	lamda_x = (double)ftmp1;
	lamda_z = (double)ftmp2;
	v0 = (double)ftmp3;

	k1[0] = v0 * r_rt->xp_ray[i];
	k1[1] = v0 * r_rt->zp_ray[i];
	k1[2] = - lamda_x / v0 / v0;
	k1[3] = - lamda_z / v0 / v0;

	//k2[0] = v0 * (r_rt->xp_ray[i] + 0.5 * h * k1[0]);
	k2[0] = v0 * (r_rt->xp_ray[i] + 0.5 * h * k1[2]);
		//Here we should not use "v0=v(sn)", but "v0=v(sn+0.5h)"
	//k2[1] = v0 * (r_rt->zp_ray[i] + 0.5 * h * k1[1]);
	k2[1] = v0 * (r_rt->zp_ray[i] + 0.5 * h * k1[3]);
	k2[2] = - lamda_x / v0 / v0;
	k2[3] = - lamda_z / v0 / v0;

	k3[0] = v0 * (r_rt->xp_ray[i] + 0.5 * h * k2[2]);
		//Here we should not use "v0=v(sn)", but "v0=v(sn+0.5h)"
	k3[1] = v0 * (r_rt->zp_ray[i] + 0.5 * h * k2[3]);
	k3[2] = - lamda_x / v0 / v0;
	k3[3] = - lamda_z / v0 / v0;

	k4[0] = v0 * (r_rt->xp_ray[i] + h * k3[2]);
		//Here we should not use "v0=v(sn)", but "v0=v(sn+h)"
	k4[1] = v0 * (r_rt->zp_ray[i] + h * k3[3]);
	k4[2] = - lamda_x / v0 / v0;
	k4[3] = - lamda_z / v0 / v0;

	/*printf ("k1_0 = %f, k1_1 = %f, k1_2 = %f, k1_3 = %f\n", k1[0], k1[1], k1[2], k1[3]);
	printf ("k2_0 = %f, k2_1 = %f, k2_2 = %f, k2_3 = %f\n", k2[0], k2[1], k2[2], k2[3]);
	printf ("k3_0 = %f, k3_1 = %f, k3_2 = %f, k3_3 = %f\n", k3[0], k3[1], k3[2], k3[3]);
	printf ("k4_0 = %f, k4_1 = %f, k4_2 = %f, k4_3 = %f, h = %f\n\n", k4[0], k4[1], k4[2], k4[3], h);*/
	r_rt->x_loc_ray[i+1] = r_rt->x_loc_ray[i] + (k1[0] + 2.*k2[0] + 2.*k3[0] + k4[0]) * h / 6.;
	r_rt->z_loc_ray[i+1] = r_rt->z_loc_ray[i] + (k1[1] + 2.*k2[1] + 2.*k3[1] + k4[1]) * h / 6.;
	r_rt->xp_ray[i+1] = r_rt->xp_ray[i] + (k1[2] + 2.*k2[2] + 2.*k3[2] + k4[2]) * h / 6.;
	r_rt->zp_ray[i+1] = r_rt->zp_ray[i] + (k1[3] + 2.*k2[3] + 2.*k3[3] + k4[3]) * h / 6.;

	dtmp1 = r_rt->x_loc_ray[i+1] - r_rt->x_loc_ray[i];
	dtmp2 = r_rt->z_loc_ray[i+1] - r_rt->z_loc_ray[i];
	r_rt->len_ray[i+1] = r_rt->len_ray[i] + sqrt (dtmp1*dtmp1 + dtmp2*dtmp2);
	dtmp1 = r_rt->len_ray[i+1] - r_rt->len_ray[i];
	r_rt->tau_ray[i+1] = r_rt->tau_ray[i] + dtmp1 / v0;
	//Here we should not use "v0=v(sn)", but average value between sn and an+h

	return 0;
}
