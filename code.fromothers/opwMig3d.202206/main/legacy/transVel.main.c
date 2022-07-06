//Sys.
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>

//Definition of Macros.
#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef FILE_NAME_MAX_LENGTH
#define FILE_NAME_MAX_LENGTH 256
#endif

#ifndef DEBUG_READ_WRITE_MODE
#define DEBUG_READ_WRITE_MODE
#endif

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y)) 
#endif

#include"/home/fb/pw3dlib/include/fbsegy.h"
#include"/home/fb/pw3dlib/include/fballoc.h"
#include"/home/fb/pw3dlib/include/fbrw.h"
#include"/home/fb/pw3dlib/include/fballoc.c"
#include"/home/fb/pw3dlib/include/fbrw.c"
#include"/home/fb/pw3dlib/include/fbsurw.c"

int main( int argc , char *argv[] )
{
	char	velxyzfn[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_210z_676y_676x_pc";
	char	velzxyfn[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_ny676_nx676_nz210.dat";
	char	velzyxfn[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_nx676_ny676_nz210.dat";

	int	ny=676;
	int	nx=676;
	int	nz=210;
	float	dy=20.0;
	float	dx=20.0;
	float	dz=20.0;

/*
	float	***vel3din_xyz	= alloc3float(nx,ny,nz);
	read_3d_float_rb(vel3din_xyz,   nz, ny, nx, velxyzfn);

	float	***vel3dout_zxy	= alloc3float(nz,nx,ny);
	for ( int iy=0; iy<ny; ++iy)
	for ( int ix=0; ix<nx; ++ix)
	for ( int iz=0; iz<nz; ++iz)
		vel3dout_zxy[iy][ix][iz] = vel3din_xyz[iz][iy][ix];
	write_3d_float_wb(vel3dout_zxy, ny, nx, nz, velzxyfn);

	float	***vel3dout_zyx	= alloc3float(nz,ny,nx);
	for ( int ix=0; ix<nx; ++ix)
	for ( int iy=0; iy<ny; ++iy)
	for ( int iz=0; iz<nz; ++iz)
		vel3dout_zyx[ix][iy][iz] = vel3din_xyz[iz][iy][ix];
	write_3d_float_wb(vel3dout_zyx, nx, ny, nz, velzyxfn);
*/

/*
	char	veltxyfn_dat[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_ny676_nx676_ntau615_dt8ms.vrms.dat";
	char	veltxyfn_su[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_ny676_nx676_ntau615_dt8ms.vrms.su";
	int	ns=615;
	int	dt=8000;
	float	***vel3d_txy	= alloc3float(ns,nx,ny);
	read_3d_float_rb(vel3d_txy, ny, nx, ns, veltxyfn_dat);

	fbsegy	segyhdr= {0} ;
        FILE *fp = NULL;
        if( (fp=fopen(veltxyfn_su,"wb")) == NULL )
        {
                printf("can not open file:%s to write.\n" , veltxyfn_su ) ;
                exit(1) ;
        }

	for ( int iy=0; iy<ny; ++iy)
	for ( int ix=0; ix<nx; ++ix)
	{
		int	cdp	= iy + 1;
		int	line	= ix + 1;
		segyhdr.cdp	= cdp;
		segyhdr.line	= line;
		segyhdr.ns	= ns;
		segyhdr.dt	= dt;
                fwrite( &segyhdr ,   sizeof(fbsegy) ,  1 , fp ) ;
                fwrite( &vel3d_txy[iy][ix][0] , sizeof(float) ,  ns , fp ) ;
	}
        fclose(fp);
*/

	char	veltyxfn_dat[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_nx676_ny676_ntau615_dt8ms.vrms.dat";
	char	veltyxfn_su[FILE_NAME_MAX_LENGTH]="/mnt/hgfs/VMshared_Kingsoft/data/salt3d/saltv_nx676_ny676_ntau615_dt8ms.vrms.su";
	int	ns=615;
	int	dt=8000;
	float	***vel3d_tyx	= alloc3float(ns,nx,ny);
	read_3d_float_rb(vel3d_tyx, nx, ny, ns, veltyxfn_dat);

	fbsegy	segyhdr= {0} ;
        FILE *fp = NULL;
        if( (fp=fopen(veltyxfn_su,"wb")) == NULL )
        {
                printf("can not open file:%s to write.\n" , veltyxfn_su ) ;
                exit(1) ;
        }

	for ( int ix=0; ix<nx; ++ix)
	for ( int iy=0; iy<ny; ++iy)
	{
		int	cdp	= iy + 1;
		int	line	= ix + 1;
		segyhdr.cdp	= cdp;
		segyhdr.line	= line;
		segyhdr.ns	= ns;
		segyhdr.dt	= dt;
                fwrite( &segyhdr ,   sizeof(fbsegy) ,  1 , fp ) ;
                fwrite( &vel3d_tyx[ix][iy][0] , sizeof(float) ,  ns , fp ) ;
	}
        fclose(fp);

	return	0;
}
