#include "fbopw3d.h"
/*
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"fballoc.h"
#include"fballoc.c"
#include"fbrw.h"
#include"fbrw.c"
*/

#define LengthMax 1024
float   ricker_t0(float t, float fm, float t0)
{
    //  Ricker FB.
    float x=pow(PI*fm*(t-t0),2);
    return exp(-x)*(1-2*x);
}

void write_1d_float_wb_long ( float *data , long n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file: %s to write\n" , fn ) ;
        exit(0) ;
    }
    fwrite( data , FLOAT_SIZE_BYTES , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}

void read_1d_float_rb_long ( float *data , long n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file: %s to read\n" , fn ) ;
        exit(0) ;
    }
    fread( data , FLOAT_SIZE_BYTES , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n",fn);
#endif
}

int main( int argc , char *argv[] )
{
    //lenovo-222.
    //char prog_dir[LengthMax]="/ssd/ssd_3.2T/fb/opwMig3d.202202/syntheticData";
    //ocean-313.
    char prog_dir[LengthMax]="/home/fb/data/syntheticData";
    char fn_vin[LengthMax]="";
    sprintf(fn_vin, "%s/velocity/fb_velocity_48mx_8mz_1000_496.dat", prog_dir);
    char fn_vout[LengthMax]="";
    sprintf(fn_vout, "%s/velocity/velocity1D_ix585_nz320_dz8m.dat", prog_dir);

    // DaQing Velocity Model.
    int nx  = 1000 ;
    int nz  = 496 ;
    float   dx  = 48.0;
    float   dz  = 8.0;
    float   **vel2d = alloc2float(nz, nx);
    read_2d_float_rb(vel2d, nx, nz, fn_vin);

    // Extract a 1-D model.
    int ivx = 585;
    int nz2 = 320 ;
    float   *vel1d_z = calloc(nz2, sizeof(float));
    for(int iz=0; iz<nz2; ++iz)
        vel1d_z[iz]   = vel2d[ivx][iz];
    //write_1d_float_wb(vel1d_z, nz2, fn_vout);

    // Get zero-offset section.
    int ns      = 4000;
    float   dt  = 1.0/1E3;  // unit of "dt" is second.
    float   fm  = 25.0;
    float   *wavelet= calloc(ns, sizeof(float));
    float   *trace  = calloc(ns, sizeof(float));
    float   *v1d_intt= calloc(ns, sizeof(float));
    float   t0  = 0;
    float   verr= 10.0;
    int     itbeg   = 0;
    int     itend   = ns;
    for(int iz=0; iz<nz2-1; ++iz)
    {
        t0  += 2*dz/vel1d_z[iz];
        if( fabs(vel1d_z[iz]-vel1d_z[iz+1])>verr )
        {
            printf("Z=%f(m) V=%f(m/s) t0=%f(ms)\n", dz*iz, vel1d_z[iz], t0*1000);
            for(int it=0; it<ns; ++it)
                wavelet[it]   = ricker_t0(dt*it, fm, t0);

            for(int it=0; it<ns; ++it)
                trace[it]   += wavelet[it];

            itend   = (int)(t0/dt + 0.5);
            for(int it=itbeg; it<itend; ++it)
                v1d_intt[it]   = vel1d_z[iz];
            itbeg = itend;
        }
    }
    for(int it=itbeg; it<ns; ++it)
        v1d_intt[it]   = vel1d_z[nz2-1];
    //write_1d_float_wb(trace, ns, "./zotrace.dat");
    char fn_vintt1d[LengthMax]="";
    sprintf(fn_vintt1d, "%s/velocity/vintt_ns4000_by_fb.dat", prog_dir);
    //write_1d_float_wb(v1d_intt, ns, fn_vintt1d);

    // output the 3-D CDP-gather and 3-D velocity model.
    int Nvtau       = ns;
    int Nvtau_first = 1;
    int Nvtau_final = Nvtau;
    int Nvtau_step  = 1;
    int Nvmy        = 800;
    int Nvline_first= 1;
    int Nvline_final= 800;
    int Nvline_step = 1;
    int Nvmx        = 800;
    int Nvcdp_first = 1;
    int Nvcdp_final = 800;
    int Nvcdp_step  = 1;
    float Dvtau     = 1.0; //ms.
    float Dvy       = 20;  //m.
    float Dvx       = 20;  //m.
    float   ***vintt3d  = alloc3float(Nvmx, Nvmy, Nvtau);
    for(int itau=0; itau<Nvtau; ++itau)
        for(int ivy=0; ivy<Nvmy; ++ivy)
            for(int ivx=0; ivx<Nvmx; ++ivx)
                vintt3d[itau][ivy][ivx] = v1d_intt[itau];

    printf("Writing the TD-Interval velocity.\n");
    char fn_vintt3d[LengthMax]="";
    sprintf(fn_vintt3d, "%s/velocity/vintt3d_ntau_4000_nvy_800_nvx_800.dat", prog_dir);
    write_3d_float_wb(vintt3d, Nvtau, Nvmy, Nvmx, fn_vintt3d);
    return 0 ;

    int Ntau        = ns;
    int Ntau_first  = 1;
    int Ntau_final  = Ntau;
    int Ntau_step   = 1;
    int Nline       = 500;
    int Nline_first = 151;
    int Nline_final = 650;
    int Nline_step  = 1;
    int Ncdp        = 500;
    int Ncdp_first  = 151;
    int Ncdp_final  = 650;
    int Ncdp_step   = 1;
    float dcdp      = 20;  //m.
    float dline     = 20;  //m.

    char fn_opw3d[LengthMax]="";
    float   ***opwData3d    = alloc3float(Ntau, Ncdp, Nline);

    for(int it=0; it<Ntau; ++it)
    {
        wavelet[it] = ricker_t0(dt*it, fm, dt*1000);
        trace[it]   = 0;
    }
    //write_1d_float_wb(wavelet, Ntau, "./wavelet_orig.dat");
    /*
    for(int it=0; it<Ntau; ++it)
        for(int jt=0; jt<=it; ++jt)
            trace[it] += wavelet[jt];
    write_1d_float_wb(trace, Ntau, "./wavelet_integration.dat");
     */
    float   dw  = 2.0*PI/(dt*Ntau);
    for(int it=0; it<Ntau; ++it)
        trace[it] = wavelet[it]*exp(-dw*dt*it);
    //write_1d_float_wb(trace, Ntau, "./wavelet_dampped.dat");
    // for impulse response!
    memset( (void*)&opwData3d[0][0][0], 0, sizeof(float)*1L*Ntau*Ncdp*Nline );
    for(int it=0; it<Ntau; ++it)
        opwData3d[Nline/2][Ncdp/2][it]   = wavelet[it];
        //opwData3d[Nline/2][Ncdp/2][it]   = trace[it];
    sprintf(fn_opw3d, "%s/data/opw3d_nline_500_ncdp_500_ntau_4000_impulseResponse.dat", prog_dir);
    printf("Writing the impulse-response OPW data.\n");
    write_3d_float_wb(opwData3d, Nline, Ncdp, Ntau, fn_opw3d);

    // constant velocity.
    for(int itau=0; itau<Nvtau; ++itau)
        for(int ivy=0; ivy<Nvmy; ++ivy)
            for(int ivx=0; ivx<Nvmx; ++ivx)
                vintt3d[itau][ivy][ivx] = 2500.;

    printf("Writing the constant-velocity.\n");
    sprintf(fn_vintt3d, "%s/velocity/vconst3d_ntau_4000_nvy_800_nvx_800.dat", prog_dir);
    //write_3d_float_wb(vintt3d, Nvtau, Nvmy, Nvmx, fn_vintt3d);
    //return 0 ;

    // produce the Ph=0 OPW data.
    memset( (void*)&opwData3d[0][0][0], 0, sizeof(float)*1L*Ntau*Ncdp*Nline );
    for(int iline=Nline_first; iline<=Nline_final; ++iline)
        for(int icdp=Ncdp_first; icdp<=Ncdp_final; ++icdp)
            for(int it=0; it<Ntau; ++it)
            {
                int iy = iline - Nline_first;
                int ix = icdp - Ncdp_first;
                opwData3d[iy][ix][it] = trace[it];
            }
    sprintf(fn_opw3d, "%s/opw3d_nline_500_ncdp_500_ntau_4000.dat", prog_dir);
    printf("Writing the Ph=0 OPW data.\n");
    //write_3d_float_wb(opwData3d, Nline, Ncdp, Ntau, fn_opw3d);


    // synthetic the non-zero offset CMP gather.
    float   *v1d_rmst= calloc(Ntau, sizeof(float));
    char fn_v1d_rmst[LengthMax]="";
    sprintf(fn_v1d_rmst, "%s/velocity/vrmst_ns4000_by_velconv.dat", prog_dir);
    read_1d_float_rb(v1d_rmst, Ntau, fn_v1d_rmst);

    float   ***vrmst3d  = alloc3float(Nvmx, Nvmy, Nvtau);
    for(int itau=0; itau<Nvtau; ++itau)
        for(int ivy=0; ivy<Nvmy; ++ivy)
            for(int ivx=0; ivx<Nvmx; ++ivx)
                vrmst3d[itau][ivy][ivx] = v1d_rmst[itau];

    printf("Writing the TD-RMS velocity.\n");
    char fn_vrmst3d[LengthMax]="";
    sprintf(fn_vrmst3d, "%s/velocity/vrmst3d_ntau_4000_nvy_800_nvx_800.dat", prog_dir);
    //write_3d_float_wb(vrmst3d, Nvtau, Nvmy, Nvmx, fn_vrmst3d);

    float   offset_x_min    = -2500.;
    float   offset_x_max    = +2500.;
    float   offset_y_min    = -2500.;
    float   offset_y_max    = +2500.;
    float   doffset_x   = 50.0;
    float   doffset_y   = 50.0;
    offset_y_min    = 0;
    offset_y_max    = 0;
    int ntr_x   = (int)( (offset_x_max-offset_x_min)/doffset_x + 1.5 );
    int ntr_y   = (int)( (offset_y_max-offset_y_min)/doffset_y + 1.5 );
    int ntr     = ntr_x*ntr_y;
    float   *offset_1d  = calloc(ntr, sizeof(float));
    printf("ntr_x=%d ntr_y=%d ntr=%d\n", ntr_x, ntr_y, ntr);

    float   **cmpGather    = alloc2float(Ntau, ntr);
    memset( (void*)&cmpGather[0][0], 0, sizeof(float)*1L*Ntau*ntr );

    int itr_g   = 0;
    for(int itr_y = 0; itr_y < ntr_y; ++ itr_y )
    {
        //printf("itr_y=%d ntr_y=%d \n", itr_y, ntr_y);
        for(int itr_x = 0; itr_x < ntr_x; ++ itr_x )
        {
            float   offset_x    = offset_x_min + doffset_x*itr_x;
            float   offset_y    = offset_y_min + doffset_y*itr_y;
            float   offset      = sqrt(offset_x*offset_x + offset_y*offset_y);
            offset_1d[itr_g]    = offset;
            float   x2      = offset*offset;
            for(int it=0; it<Ntau-1; ++it)
            {
                float   vrms    = v1d_rmst[it];
                float   v2      = vrms*vrms;
                float   t0      = dt*it;
                float   tnow    = sqrt(t0*t0+x2/v2);

                if( fabs(v1d_intt[it]-v1d_intt[it+1])>verr )
                {
                    if ( 0 == itr_y && 0 == itr_x )
                    {
                        printf("t0=%f tnow=%f vrms=%f\n", t0, tnow, vrms);
                    }
                    for(int it=0; it<ns; ++it)
                        wavelet[it]   = ricker_t0(dt*it, fm, tnow);

                    for(int it=0; it<ns; ++it)
                        cmpGather[itr_g][it]    += wavelet[it];
                }
            }
            ++ itr_g;
        }
    }
    printf("Writing the CMP data.\n");
    char fn_cmp[LengthMax]="";
    sprintf(fn_cmp, "%s/cmpgather_ns_4000_ntr_10201.dat", prog_dir);
    write_2d_float_wb(cmpGather, ntr, Ntau, fn_cmp);

    // tau-phr transform
    int nphr, ntau, itaus, itaue;
    float   dphr, phrmin, phrmax;
    float   dt_ms   = 1.0;
    nphr    = 61;
    dphr    = 10.0/1E6;
    phrmin  = 0;
    phrmax  = phrmin + dphr*(nphr-1);
    ntau    = 4000;
    itaus   = 0;
    itaue   = ntau;
    float   **taup2d    = alloc2float(ntau, nphr);
    memset( (void*)&taup2d[0][0], 0, sizeof(float)*1L*ntau*nphr );
    // the tau-phr transform proposed by Feng et al., 20111.
    tauph_transform_3d_to_2d_robust(ns, ntau, nphr,
            ntr, itaus, itaue,
            cmpGather, taup2d,
            v1d_rmst, offset_1d, dt_ms, phrmin, dphr);
    printf("Writing the classical tau-phr data.\n");
    char fn_taup2d[LengthMax]="";
    sprintf(fn_taup2d, "%s/taup2d_ns_4000_nphr_101.fb.dat", prog_dir);
    write_2d_float_wb(taup2d, nphr, ntau, fn_taup2d);

    // the theoretical tau-p domain ellipse curve.
    memset( (void*)&taup2d[0][0], 0, sizeof(float)*1L*ntau*nphr );
    for(int it=0; it<ntau; ++it)
    {
        float   vrms    = v1d_rmst[it];
        float   v2      = vrms*vrms;
        float   t0      = dt*it;
        if( fabs(v1d_intt[it]-v1d_intt[it+1])>verr )
        {
            printf("t0=%f vrms=%f\n", t0, vrms);
            for(int iphr=0; iphr<nphr; ++iphr)
            {
                float   phr     = phrmin + dphr*(iphr-1);
                float   tmp     = 1 - phr*phr*v2/4.0;
                if ( tmp > 0 )
                {
                    float   tnow    = t0*sqrt(tmp);
                    for(int it=0; it<ns; ++it)
                        wavelet[it]   = ricker_t0(dt*it, fm, tnow);

                    for(int it=0; it<ns; ++it)
                        taup2d[iphr][it]    += wavelet[it];
                }
            }
        }
    }
    printf("Writing the theoretical tau-phr data.\n");
    sprintf(fn_taup2d, "%s/taup2d_ns_4000_nphr_101.exact.dat", prog_dir);
    write_2d_float_wb(taup2d, nphr, ntau, fn_taup2d);
    //return 0 ;

    float   Byte_GB = 1.0*1024*1024*1024;
    printf("Forming the 4-D OPW data.\n");
    printf("ntau=%d Ncdp=%d Nline=%d nphr=%d\n", ntau, Ncdp, Nline, nphr);
    printf("Size of 4-D OPW data is %f(GB).\n", sizeof(float)*1L*ntau*Ncdp*Nline*nphr/Byte_GB);
    float   ****taup4d    = alloc4float(ntau, Ncdp, Nline, nphr);
    memset( (void*)&taup4d[0][0][0][0], 0, sizeof(float)*1L*ntau*Ncdp*Nline*nphr );
    for(int iphr=0; iphr<nphr; ++iphr)
    {
        for(int iline=Nline_first; iline<=Nline_final; ++iline)
        {
            for(int icdp=Ncdp_first; icdp<=Ncdp_final; ++icdp)
            {
                int iy = iline - Nline_first;
                int ix = icdp - Ncdp_first;
                for(int it=0; it<ntau; ++it)
                    taup4d[iphr][iy][ix][it] = taup2d[iphr][it];
            }
        }
    }
    printf("Writing the 4-D OPW data.\n");
    char fn_opw4d[LengthMax]="";
    sprintf(fn_opw4d, "%s/opw4d_nphr_61_nline_600_ncdp_600_ntau_4000.exact.dat", prog_dir);
    write_1d_float_wb_long ( &taup4d[0][0][0][0], 1L*ntau*Ncdp*Nline*nphr, fn_opw4d);

    return 0 ;
}
