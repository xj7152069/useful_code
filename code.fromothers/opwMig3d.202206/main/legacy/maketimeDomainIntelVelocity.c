#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"fballoc.h"
#include"fballoc.c"
#include"fbrw.h"
#include"fbrw.c"

#define LengthMax 1024

int main( int argc , char *argv[] )
{
	//wpi-200.
	char fn_vin[LengthMax]="/data2_2/fb/LJ2017_cmp/velocity/rmsVel_ncdp2000_nline800_ntau3010.dat";
	char fn_vout[LengthMax]="/data2_2/fb/LJ2017_cmp/timeDomainIntVel_ncdp_nline_ntau.dat";
	char fn_vout_txy[LengthMax]="/data2_2/fb/LJ2017_cmp/timeDomainIntVel_txy.dat";
	/*
	//lenovo-222.
	char fn_vin[LengthMax]="/home/fb/LJ2017_cmp/LJ_rms_new.dat";
	char fn_vout[LengthMax]="/ssd/ssd_raid5/fb/LJ2017_cmp/timeDomainIntVel_ncdp_nline_ntau.dat";
	char fn_vout_txy[LengthMax]="/home/fb/LJ2017_cmp/timeDomainIntVel_txy.dat";
	 */
	//char fn_vout[LengthMax]="/ssd/ssd_raid5/fb/LJ2017_cmp/timeDomainIntVel_ncdp_nline_ntau.dat";
	int	nx	= 233 ;
	int	ny	= 337 ;
	int	ns	= 3010 ;
    int cdpmin  = 2001;
    int cdpmax  = 4000;
	int	linemin	= 1800 ;
	int	linemax	= 2600 ;
	int	ncdp	= 4641 ;
	int	nline	= 3361 ;
	int	ntau	= ns ;
	ncdp	= cdpmax - cdpmin + 1 ;
	nline	= linemax - linemin + 1 ;
    printf("ncdp=%d nline=%d VelSize=%f(GB).\n", ncdp, nline, 4.0*ncdp*nline*ntau/1024./1024./1024.);

	float	***v_in = alloc3float(ns, nx, ny);
	float	***v_out= alloc3float(ncdp, nline, ntau);
	float	***v_out_txy    = alloc3float(ntau, ncdp, nline);

    printf("Reading the sparse RMs velocity.\n");
    // read the sparse RMs velocity.
	read_3d_float_rb(v_in, ny, nx, ns, fn_vin);

    // calculate the 1-D RMS velocity.
    //#pragma omp parallel for
    for ( int itau = 0 ; itau < ntau ; ++ itau )
    {
        float   vsum    = 0;
        for ( int iy = 0 ; iy < ny ; ++ iy )
            for ( int ix = 0 ; ix < nx ; ++ ix )
                vsum    += v_in[iy][ix][itau];

        float   vavg    = vsum/(1.*ny*nx);
        int iy  = ny/2;
        int ix  = nx/2;
        if( 0 == itau%200 )
            printf("itau=%d vavg=%f vsum=%f iy=%d ix=%d vrms=%f.\n", itau, vavg, vsum, iy, ix, v_in[iy][ix][itau]);
        for ( int iline = 0 ; iline < nline ; ++ iline )
            for ( int icdp = 0 ; icdp < ncdp ; ++ icdp )
                    v_out[itau][iline][icdp]    = vavg;

        /*
        for ( int iline = 0 ; iline < nline ; ++ iline )
            for ( int icdp = 0 ; icdp < ncdp ; ++ icdp )
                    v_out_txy[iline][icdp][itau]    = vavg;
        */
    }

    printf("Writing the TD-Interval velocity.\n");
	write_3d_float_wb(v_out, ntau, nline, ncdp, fn_vout);
	//write_3d_float_wb(v_out_txy, nline, ncdp, ntau, fn_vout_txy);

	free3float(v_in);
	free3float(v_out);

	return 0 ;
}
