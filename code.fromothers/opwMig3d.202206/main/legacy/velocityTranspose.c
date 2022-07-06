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
    char fn_vout[LengthMax]="/data2_2/fb/LJ2017_cmp/timeDomainIntVel_nx_ny_ns.dat";
    int	nx	= 2000 ;
    int	ny	= 800 ;
    int	ns	= 3010 ;
    printf("nx=%d ny=%d VelSize=%f(GB).\n", nx, ny, 4.0*nx*ny*ns/1024./1024./1024.);

    float	***v_in = alloc3float(ns, nx, ny);
    float	***v_out= alloc3float(nx, ny, ns);

    printf("Reading the RMs velocity.\n");
    // read the sparse RMs velocity.
    read_3d_float_rb(v_in, ny, nx, ns, fn_vin);

    // calculate the 1-D RMS velocity.
    //#pragma omp parallel for
    for ( int it = 0 ; it < ns ; ++ it )
    {
        if( 0 == it%200 )
            printf("it=%d \n",it);
        for ( int iy = 0 ; iy < ny ; ++ iy )
            for ( int ix = 0 ; ix < nx ; ++ ix )
                v_out[it][iy][ix]    = v_in[iy][ix][it];

    }

    printf("Writing the TD-Interval velocity.\n");
    write_3d_float_wb(v_out, ns, ny, nx, fn_vout);

    free3float(v_in);
    free3float(v_out);

    return 0 ;
}
