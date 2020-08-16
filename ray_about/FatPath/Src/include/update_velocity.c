int update_velocity (model_rt m_rt, float **vel, float **d_slw, float epsi_update)
{
	int ix, iz;
	float ftmp1;

	for (ix=0; ix<m_rt.nx; ix++)
		for (iz=0; iz<m_rt.nz; iz++)
		{
			//vel[ix][iz] = 1. / (d_slw[ix][iz] + 1. / vel[ix][iz]);

			ftmp1 = 1. / (d_slw[ix][iz] + 1. / vel[ix][iz]);

			/*
			*/
			if ((ftmp1 - vel[ix][iz])/vel[ix][iz] > epsi_update)
				vel[ix][iz] = vel[ix][iz] + vel[ix][iz] * epsi_update;
			else if ((vel[ix][iz] - ftmp1)/vel[ix][iz] > epsi_update)
				vel[ix][iz] = vel[ix][iz] - vel[ix][iz] * epsi_update;
			else
				vel[ix][iz] = ftmp1;
		}

	return 0;
}
