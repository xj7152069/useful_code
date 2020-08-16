
int calculate_kernel_ray (float **kernel, model_rt m_rt, ray_rt r_rt,
	int n_point, float ds)
{
	int ix, iz, i, ia, ja, ib, jb, itmp1, itmp2;
	float cx, cz, ex, ez, ftmp1, ftmp2, ftmp3, ftmp4;

	for (ix=0; ix<m_rt.nx; ix++)
		for (iz=0; iz<m_rt.nz; iz++)
			kernel[ix][iz] = 0;
	
	ia = (int)(r_rt.x_loc_ray[0] / m_rt.dx);
	ja = (int)(r_rt.z_loc_ray[0] / m_rt.dz);
	for (i=1; i<n_point; i++)
	{
		ib = (int)(r_rt.x_loc_ray[i] / m_rt.dx);
		jb = (int)(r_rt.z_loc_ray[i] / m_rt.dz);


		if (ia == ib && ja == jb)
			kernel[ia][ja] += ds;
		else if (ia == ib && ja != jb)
		{
			itmp1 = (ja>jb)?ja:jb;
			cz = m_rt.dz * (float)itmp1;
			cx = r_rt.x_loc_ray[i] - (r_rt.x_loc_ray[i-1] - r_rt.x_loc_ray[i]) *
				(r_rt.z_loc_ray[i] - cz) / (r_rt.z_loc_ray[i-1] - r_rt.z_loc_ray[i]);
			ftmp1 = cx - r_rt.x_loc_ray[i-1];
			ftmp2 = cz - r_rt.z_loc_ray[i-1];
			ftmp1 = sqrt (ftmp1 * ftmp1 + ftmp2 * ftmp2);
			//kernel[ia][ja] += sqrt (ftmp1 * ftmp1 + ftmp2 * ftmp2);
			kernel[ia][ja] += ftmp1;

			kernel[ib][jb] += ds - ftmp1;
		}
		else if (ia != ib && ja == jb)
		{
			itmp1 = (ia>ib)?ia:ib;
			cx = m_rt.dx * (float)itmp1;
			cz = r_rt.z_loc_ray[i] - (r_rt.z_loc_ray[i-1] - r_rt.z_loc_ray[i]) *
				(r_rt.x_loc_ray[i] - cx) / (r_rt.x_loc_ray[i-1] - r_rt.x_loc_ray[i]);
			ftmp1 = cx - r_rt.x_loc_ray[i-1];
			ftmp2 = cz - r_rt.z_loc_ray[i-1];
			ftmp1 = sqrt (ftmp1 * ftmp1 + ftmp2 * ftmp2);
			kernel[ia][ja] += ftmp1;

			kernel[ib][jb] += ds - ftmp1;
		}
		else
		{
			itmp1 = (ia>ib)?ia:ib;
			itmp2 = (ja>jb)?ja:jb;
			cx = m_rt.dx * (float)itmp1;
			ez = m_rt.dz * (float)itmp2;
			cz = r_rt.z_loc_ray[i] - (r_rt.z_loc_ray[i-1] - r_rt.z_loc_ray[i]) *
				(r_rt.x_loc_ray[i] - cx) / (r_rt.x_loc_ray[i-1] - r_rt.x_loc_ray[i]);
			ex = r_rt.x_loc_ray[i] - (r_rt.x_loc_ray[i-1] - r_rt.x_loc_ray[i]) *
				(r_rt.z_loc_ray[i] - ez) / (r_rt.z_loc_ray[i-1] - r_rt.z_loc_ray[i]);
			itmp1 = (int)(0.5 * (cx + ex) / m_rt.dx);
			itmp2 = (int)(0.5 * (cz + ez) / m_rt.dz);
			ftmp1 = cx - r_rt.x_loc_ray[i-1];
			ftmp2 = cz - r_rt.z_loc_ray[i-1];
			ftmp3 = sqrt (ftmp1 * ftmp1 + ftmp2 * ftmp2);   //|AC|
			ftmp1 = ex - r_rt.x_loc_ray[i-1];
			ftmp2 = ez - r_rt.z_loc_ray[i-1];
			ftmp4 = sqrt (ftmp1 * ftmp1 + ftmp2 * ftmp2);   //|AE|
			if (ftmp3 < ftmp4)
			{
				kernel[ia][ja] += ftmp3;
				kernel[ib][jb] += ds - ftmp4;
				kernel[itmp1][itmp2] += ftmp4 - ftmp3;
			}
			else
			{
				kernel[ia][ja] += ftmp4;
				kernel[ib][jb] += ds - ftmp3;
				kernel[itmp1][itmp2] += ftmp3 - ftmp4;
			}
		}

		ia = ib;
		ja = jb;

	}

	return 0;
}


