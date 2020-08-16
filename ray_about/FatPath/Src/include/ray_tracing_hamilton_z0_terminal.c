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

int ray_tracing_hamilton_z0_terminal (model_rt *m_rt, ray_rt *r_rt, float **vel, float ds, int *i_edge,
	float **tab_v, float **tab_vx, int n_table, float z0, int *n_point)
{
	int i, n_yn;
	float v0;

	/*
	x_min = 0;
	x_max = (double)m_rt->dx * (double)(m_rt->nx - 1);
	z_min = 0;
	z_max = (double)m_rt->dz * (double)(m_rt->nz - 1);
	*/

	n_yn = get_vel (&v0, m_rt->x_src, m_rt->z_src, m_rt->nx, m_rt->nz, m_rt->dx, m_rt->dz, m_rt->x_min,
		m_rt->z_min, m_rt->x_max, m_rt->z_max, vel, tab_v, tab_vx, n_table);
	if (0 != n_yn)
	{
		printf ("Calculating velocity of source location error.\n");
		return 1;
	}

	r_rt->xp_ray[0] = cos (m_rt->thita) / v0;
	r_rt->zp_ray[0] = sin (m_rt->thita) / v0;
	r_rt->x_loc_ray[0] = m_rt->x_src;
	r_rt->z_loc_ray[0] = m_rt->z_src;
	r_rt->tau_ray[0] = 0;
	r_rt->len_ray[0] = 0;

	//n_point = 0;
	i = 0;
	while (1)
	{
		if (r_rt->x_loc_ray[i] <= m_rt->x_min)
		{
			*i_edge = -1;
			break;
		}
		else if (r_rt->x_loc_ray[i] >= m_rt->x_max)
		{
			*i_edge = 1;
			break;
		}
		else if (r_rt->z_loc_ray[i] <= m_rt->z_min)
		{
			*i_edge = -2;
			break;
		}
		//else if (r_rt->z_loc_ray[i] > z_max)
		else if (r_rt->z_loc_ray[i] >= z0)
		{
			*i_edge = 2;
			break;
		}

		cell_ray_tracing (m_rt, r_rt, vel, ds, i, tab_v, tab_vx, n_table);
		
		i ++;
	}

	if (2 == *i_edge)
		*n_point = i + 1;
	else
		*n_point = i;	//The Last Point Out of Model Already.

	return 0;

}

