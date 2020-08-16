

typedef struct
{
	int nx;
	int nz;
	float dx;
	float dz;
	float x_src;
	float z_src;
	float x_rec;
	float z_rec;
	float thita;
	float x_min;
	float x_max;
	float z_min;
	float z_max;
}model_rt;

typedef struct
{
	float *x_loc_ray;
	float *z_loc_ray;
	float *xp_ray;
	float *zp_ray;
	float *len_ray;
	float *tau_ray;
}ray_rt;


typedef struct
{
	float *vel_s;
	float *vel_s_nn;
}ray_lamda;
