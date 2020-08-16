
/*----------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/

/*Subroutine for P2P ray tracing.*/
/*This is a self recursion function.*/

int p2p_ray_tracing_hamilton (model_rt *m_rt, ray_rt *r_rt, float **vel, float ds,
	float **tab_v, float **tab_vx, int n_table, float x_alpha, float z_alpha, float *thita_left, float *thita_right,
	int n_thita, float *distance_left, float *distance_right,
	float *epsi_vicinity, int *i_recursion, int n_recursion)
{
	/*
	Reciever located on a horizontal surface, rays terminated when they hit this surface.
	In each recursion, there are minimum and maximum take-off angles, between this 2 angles, rays
		are shot regularly. The ray whose terminal point located at left side of reciever and
		closest to reciever is then found, so as the ray that terminal point lacated at right
		side of reciever. If any of this 2 rays is close enough with reciever, the P2P ray is
		just this ray, if not, this 2 rays are passed to next recursion.
	*/
	int i_thita, n_yn, i_edge, n_point;
	float thita_min, thita_max, d_thita, ftmp1, ftmp2;

	/*Left ray*/
	if (fabs (*distance_left) < (*epsi_vicinity))
	{
		m_rt->thita = *thita_left;
		return 0;
	}
	/*Right ray*/
	else if (fabs (*distance_right) < (*epsi_vicinity))
	{
		m_rt->thita = *thita_right;
		return 0;
	}

	thita_min = *thita_left;
	thita_max = *thita_right;
	d_thita = (thita_max - thita_min) / (float)(n_thita - 1);

	/*Ray tracing n_thita times between letf and right take-off angles.*/
	for (i_thita=0; i_thita<n_thita; i_thita++)
	{
		m_rt->thita = thita_min + d_thita * (float)i_thita;

		n_yn = ray_tracing_hamilton_z0_terminal (m_rt, r_rt, vel, ds, &i_edge, tab_v, tab_vx,
			n_table, z_alpha, &n_point);
		if (0 != n_yn)
		{
			printf ("FUNCTION 'ray_tracing_hamilton_z0_terminal' error.\n");
			return 1;
		}
		
		if (2 != i_edge)
		{

			ftmp1 = r_rt->x_loc_ray[n_point-1] - x_alpha;
			ftmp2 = r_rt->z_loc_ray[n_point-1] - z_alpha;
			ftmp2 = sqrt (ftmp1 * ftmp1 + ftmp2 * ftmp2);
			if (0 > ftmp1 && ftmp2 < fabs(*distance_left))
			{
				*distance_left = - ftmp2;
				*thita_left = m_rt->thita;
			}
			else if (0 < ftmp1 && ftmp2 < *distance_right)
			{
				*distance_right = ftmp2;
				*thita_right = m_rt->thita;
			}
			
		}
		else
		{

			ftmp1 = r_rt->x_loc_ray[n_point-2] + (r_rt->x_loc_ray[n_point-1] - r_rt->x_loc_ray[n_point-2])
				* (z_alpha - r_rt->z_loc_ray[n_point-2]) / (r_rt->z_loc_ray[n_point-1] - r_rt->z_loc_ray[n_point-2]);
			r_rt->x_loc_ray[n_point-1] = ftmp1;
			r_rt->z_loc_ray[n_point-1] = z_alpha;

			ftmp1 = r_rt->x_loc_ray[n_point-1] - x_alpha;
			if (0 > ftmp1 && ftmp1 > *distance_left)
			{
				*distance_left = ftmp1;
				*thita_left = m_rt->thita;
			}
			else if (0 < ftmp1 && ftmp1 < *distance_right)
			{
				*distance_right = ftmp1;
				*thita_right = m_rt->thita;
			}
			
		}
	}

	(*i_recursion) ++;

	/*If recursion time is too large, the P2P condition is relaxed.*/
	if ((*i_recursion) >= n_recursion)
	{
		(*epsi_vicinity) *= 2.;
		*i_recursion = 0;
	}

	/*Subroutine call itself, so this is called "self recursion".*/
	n_yn = p2p_ray_tracing_hamilton (m_rt, r_rt, vel, ds, tab_v, tab_vx, n_table, x_alpha, z_alpha,
		thita_left, thita_right, n_thita, distance_left, distance_right,
		epsi_vicinity, i_recursion, n_recursion);
	if (0 != n_yn)
	{
		printf ("FUNCTION 'ray_tracing_point_to_vicinity' error.\n");
		return 1;
	}

	return 0;
}

