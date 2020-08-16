
/*
int solver_sirt (float **kernel, model_rt *m_rt, int **n_cover,
	ray_info_ini *rii, float **d_slw, int i_src, int i_rec,
	float d_travel_time_with_d_slw)
*/
int solver_sirt (float **kernel, model_rt *m_rt, int **n_cover,
	float **data_diff, float **d_slw, int i_src, int i_rec,
	float data_diff_here)
{
	int i, j;
	float kernel_quadratic_sum, d_travel_time, d_slw_tmp;

	//area_grid = m_rt->dx * m_rt->dz;

	kernel_quadratic_sum = 0;
	for (i=0; i<m_rt->nx; i++)
	{
		for (j=0; j<m_rt->nz; j++)
		{
			if (1e-10 > kernel[i][j])
				continue;

			kernel_quadratic_sum = kernel_quadratic_sum + pow (kernel[i][j], 2.);
		}
	}

	/*
	d_travel_time = rii->csi[i_src].travel_time_cal[i_rec]
		- d_travel_time_with_d_slw;
	*/
	d_travel_time = data_diff[i_src][i_rec]
		- data_diff_here;

	for (i=0; i<m_rt->nx; i++)
	{
		for (j=0; j<m_rt->nz; j++)
		{
			if (1e-10 > kernel[i][j])
				continue;

			d_slw_tmp = d_travel_time * kernel[i][j] / kernel_quadratic_sum
				/ (float)n_cover[i][j];
			d_slw[i][j] = d_slw[i][j] + d_slw_tmp;
			/*if (d_slw[i][j] > 0)
				d_slw[i][j] = 0;
			else if (d_slw[i][j] < -0.00018)
				d_slw[i][j] = -0.00018;*/
		}
	}

	return 0;
}

int solver_sirt_low_computation (float **kernel, model_rt *m_rt, int **n_cover,
	float **data_diff, float **d_slw, int i_src, int i_rec, float **kernel_quadratic_sum,
	float data_diff_here)
{
	int i, j;
	float d_travel_time, d_slw_tmp;

	//area_grid = m_rt->dx * m_rt->dz;

	/*
	d_travel_time = rii->csi[i_src].travel_time_cal[i_rec]
		- d_travel_time_with_d_slw;
	*/
	d_travel_time = data_diff[i_src][i_rec]
		- data_diff_here;

	for (i=0; i<m_rt->nx; i++)
	{
		for (j=0; j<m_rt->nz; j++)
		{
			if (1e-10 > kernel[i][j])
				continue;

			d_slw_tmp = d_travel_time * kernel[i][j] / kernel_quadratic_sum[i_src][i_rec]
				/ (float)n_cover[i][j];
			d_slw[i][j] = d_slw[i][j] + d_slw_tmp;
			/*if (d_slw[i][j] > 0)
				d_slw[i][j] = 0;
			else if (d_slw[i][j] < -0.00018)
				d_slw[i][j] = -0.00018;*/
		}
	}

	return 0;
}
