int calculate_data_diff_here (float *data_diff_here, float **kernel, float **d_slw, model_rt m_rt)
{
	int ix, iz;

	*data_diff_here = 0;

	for (ix=0; ix<m_rt.nx; ix++)
		for (iz=0; iz<m_rt.nz; iz++)
			*data_diff_here += kernel[ix][iz] * d_slw[ix][iz];

	return 0;
}
