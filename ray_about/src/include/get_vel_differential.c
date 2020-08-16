int get_vel_differential_vel_12O (float *dVdx, float *dVdz, float *V, float *dVdxx, float *dVdxz,
	float *dVdzz, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, float **Tab_Vxx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*V = 0;
	*dVdx = 0;
	*dVdz = 0;
	*dVdxx = 0;
	*dVdxz = 0;
	*dVdzz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
		{
			*V += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_V[Lz][k];
			*dVdx += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_V[Lz][k];
			*dVdz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vx[Lz][k];
			*dVdxx += Vel[Nx0+i][Nz0+k] * Tab_Vxx[Lx][i] * Tab_V[Lz][k];
			*dVdxz += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_Vx[Lz][k];
			*dVdzz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vxx[Lz][k];
		}

	*dVdx = *dVdx / dx;
	*dVdz = *dVdz / dz;
	*dVdxx = *dVdxx / (dx * dx);
	*dVdxz = *dVdxz / (dx * dz);
	*dVdzz = *dVdzz / (dz * dz);

	return 0;
}

int get_vel_differential_xz_12O (float *dVdx, float *dVdz, float *dVdxx, float *dVdxz,
	float *dVdzz, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, float **Tab_Vxx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*dVdx = 0;
	*dVdz = 0;
	*dVdxx = 0;
	*dVdxz = 0;
	*dVdzz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
		{
			*dVdx += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_V[Lz][k];
			*dVdz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vx[Lz][k];
			*dVdxx += Vel[Nx0+i][Nz0+k] * Tab_Vxx[Lx][i] * Tab_V[Lz][k];
			*dVdxz += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_Vx[Lz][k];
			*dVdzz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vxx[Lz][k];
		}

	*dVdx = *dVdx / dx;
	*dVdz = *dVdz / dz;
	*dVdxx = *dVdxx / (dx * dx);
	*dVdxz = *dVdxz / (dx * dz);
	*dVdzz = *dVdzz / (dz * dz);

	return 0;
}

int get_vel_differential_vel_2O (float *V, float *dVdxx, float *dVdxz,
	float *dVdzz, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, float **Tab_Vxx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*V = 0;
	*dVdxx = 0;
	*dVdxz = 0;
	*dVdzz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
		{
			*V += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_V[Lz][k];
			*dVdxx += Vel[Nx0+i][Nz0+k] * Tab_Vxx[Lx][i] * Tab_V[Lz][k];
			*dVdxz += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_Vx[Lz][k];
			*dVdzz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vxx[Lz][k];
		}

	*dVdxx = *dVdxx / (dx * dx);
	*dVdxz = *dVdxz / (dx * dz);
	*dVdzz = *dVdzz / (dz * dz);

	return 0;
}

int get_vel_differential_xz_2O (float *dVdxx, float *dVdxz,
	float *dVdzz, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, float **Tab_Vxx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*dVdxx = 0;
	*dVdxz = 0;
	*dVdzz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
		{
			*dVdxx += Vel[Nx0+i][Nz0+k] * Tab_Vxx[Lx][i] * Tab_V[Lz][k];
			*dVdxz += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_Vx[Lz][k];
			*dVdzz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vxx[Lz][k];
		}

	*dVdxx = *dVdxx / (dx * dx);
	*dVdxz = *dVdxz / (dx * dz);
	*dVdzz = *dVdzz / (dz * dz);

	return 0;
}

int get_vel_differential_vel (float *dVdx, float *dVdz, float *V, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*V = 0;
	*dVdx = 0;
	*dVdz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
		{
			*V += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_V[Lz][k];
			*dVdx += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_V[Lz][k];
			*dVdz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vx[Lz][k];
		}

	*dVdx = *dVdx / dx;
	*dVdz = *dVdz / dz;

	return 0;
}

int get_vel (float *V, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*V = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
			*V += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_V[Lz][k];

	return 0;
}

int get_vel_differential_xz (float *dVdx, float *dVdz, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*dVdx = 0;
	*dVdz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
		{
			*dVdx += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_V[Lz][k];
			*dVdz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vx[Lz][k];
		}

	*dVdx = *dVdx / dx;
	*dVdz = *dVdz / dz;

	return 0;
}

int get_vel_differential_x (float *dVdx, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*dVdx = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
			*dVdx += Vel[Nx0+i][Nz0+k] * Tab_Vx[Lx][i] * Tab_V[Lz][k];

	*dVdx = *dVdx / dx;

	return 0;
}

int get_vel_differential_z (float *dVdz, float x, float z, int Nx, int Nz,
	float dx, float dz, float Xmin, float Zmin, float Xmax, float Zmax, float **Vel,
	float **Tab_V, float **Tab_Vx, int Ntable)
{
	int INx, INz, i, k;
	int Nx0, Nz0, Lx, Lz;

	INx = (int) ((x -  Xmin) / dx);
	INz = (int) ((z -  Zmin) / dz);
	Lx  = (int) (((x - Xmin) / dx - INx) * Ntable);
	Lz  = (int) (((z - Zmin) / dz - INz) * Ntable);

	if(Lx < 0)
		Lx = 0;
	if(Lz < 0)
		Lz = 0;
	if(INx > Nx - 3)
		INx = Nx - 3;
	if(INz > Nz - 3)
		INz = Nz - 3;
	if(INx < 1)
		INx = 1;
	if(INz < 1)
		INz = 1;

	Nx0 = INx - 1;
	Nz0 = INz - 1;

	*dVdz = 0;

	for(i=0; i<4; i++)
		for(k=0; k<4; k++)
			*dVdz += Vel[Nx0+i][Nz0+k] * Tab_V[Lx][i] * Tab_Vx[Lz][k];

	*dVdz = *dVdz / dz;

	return 0;
}

int extrpolation_coefficient_1O (float **Tab_V, float **Tab_Vx, int Nmax)
{
	int i;
	float x;

	for(i=0; i<Nmax; i++)
	{
		x = i * 1.0 / (Nmax - 1);
		Tab_V[i][0] = -0.5 * x * x * x + x * x - 0.5*x;
		Tab_V[i][1] =  1.5 * x * x * x - 2.5 * x * x + 1;
		Tab_V[i][2] = -1.5 * x * x * x +   2 * x * x + 0.5 * x;
		Tab_V[i][3] =  0.5 * x * x * x - 0.5 * x * x;

		Tab_Vx[i][0] = -1.5 * x * x + 2 * x - 0.5;
		Tab_Vx[i][1] =  4.5 * x * x - 5 * x;
		Tab_Vx[i][2] = -4.5 * x * x + 4 * x + 0.5;
		Tab_Vx[i][3] =  1.5 * x * x -   x;
	}

	return 0;
}

int extrpolation_coefficient_12O (float **Tab_V, float **Tab_Vx, float **Tab_Vxx, int Nmax)
{
	int i ;
	float x ;

	for(i=0;i<Nmax;i++)
	{
		x = i * 1.0 / (Nmax - 1);

		Tab_V[i][0] = - 0.5 * x * x * x +     x * x - 0.5 * x;
		Tab_V[i][1] = 1.5 * x * x * x - 2.5 * x * x + 1;
		Tab_V[i][2] = - 1.5 * x * x * x +   2 * x * x + 0.5 * x;
		Tab_V[i][3] = 0.5 * x * x * x - 0.5 * x * x;

		Tab_Vx[i][0] = - 1.5 * x * x + 2 * x - 0.5;
		Tab_Vx[i][1] = 4.5 * x * x - 5 * x;
		Tab_Vx[i][2] = - 4.5 * x * x + 4 * x + 0.5;
		Tab_Vx[i][3] = 1.5 * x * x - x;

		Tab_Vxx[i][0] = - 3 * x + 2;
		Tab_Vxx[i][1] = 9 * x - 5;
		Tab_Vxx[i][2] = - 9 * x + 4;
		Tab_Vxx[i][3] = 3 * x - 1;
	}

	return 0;
}
