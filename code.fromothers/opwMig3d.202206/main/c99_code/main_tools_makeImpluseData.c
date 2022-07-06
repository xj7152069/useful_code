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
    //char prog_dir[LengthMax]="/ssd/ssd_raid0/fb/syntheticData";
    //ocean-313.
    char prog_dir[LengthMax]="/home/fb/data/syntheticData";
    // Get zero-offset section.
    int     ns  = 4000;
    float   dt  = 1.0/1E3;  // unit of "dt" is second.
    float   fm  = 25.0;
    float   t0  = 1.0;
    float   *wavelet= calloc(ns, sizeof(float));
    float   *trace  = calloc(ns, sizeof(float));
    for(int it=0; it<ns; ++it)
    {
        wavelet[it] = ricker_t0(dt*it, fm, t0);
        trace[it]   = 0;
    }
    //write_1d_float_wb(wavelet, Ntau, "./wavelet_orig.dat");
    for(int it=0; it<ns; ++it)
        for(int jt=0; jt<=it; ++jt)
            trace[it] += wavelet[jt];
    //write_1d_float_wb(trace, ns, "./wavelet_integration.dat");
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


    /*
    float   dw  = 2.0*PI/(dt*Ntau);
    for(int it=0; it<Ntau; ++it)
        trace[it] = wavelet[it]*exp(-dw*dt*it);
    //write_1d_float_wb(trace, Ntau, "./wavelet_dampped.dat");
    //
     */
    // for impulse response!
    char fn_opw3d[LengthMax]="";
    float   ***opwData3d    = alloc3float(Ntau, Ncdp, Nline);

    /*
    // for single impulse response!
    memset( (void*)&opwData3d[0][0][0], 0, sizeof(float)*1L*Ntau*Ncdp*Nline );
    for(int it=0; it<Ntau; ++it)
        opwData3d[Nline/2][Ncdp/2][it]   = trace[it];
        //opwData3d[Nline/2][Ncdp/2][it]   = wavelet[it];
    sprintf(fn_opw3d, "%s/data/opw3d_nline_500_ncdp_500_ntau_4000_single_impulseResponse.dat", prog_dir);
    printf("Writing the single impulse-response OPW data.\n");
    //write_1d_float_wb(trace, ns, "./zotrace.dat");
    //write_3d_float_wb(opwData3d, Nline, Ncdp, Ntau, fn_opw3d);
    //return 0 ;
     */

    // for multiple impulse response!
    memset( (void*)&opwData3d[0][0][0], 0, sizeof(float)*1L*Ntau*Ncdp*Nline );
    memset( (void*)&trace[0], 0, sizeof(float)*1L*ns);
    int nref = 4;
    float dtref = 0.5;
    for(int iref=1; iref<nref; ++iref)
    {
        float   t0now  = dtref*iref;
        for(int it=0; it<ns; ++it)
            wavelet[it]   = ricker_t0(dt*it, fm, t0now);

        memset( (void*)&trace[0], 0, sizeof(float)*1L*ns);

        for(int it=0; it<ns; ++it)
            for(int jt=0; jt<=it; ++jt)
                trace[it] += wavelet[jt];

        for(int it=0; it<Ntau; ++it)
            opwData3d[Nline/2][Ncdp/2][it]   += wavelet[it];
            //opwData3d[Nline/2][Ncdp/2][it]   += trace[it];
    }
    write_1d_float_wb(opwData3d[Nline/2][Ncdp/2], ns, "./zotrace.dat");
    sprintf(fn_opw3d, "%s/data/opw3d_nline_500_ncdp_500_ntau_4000_multi_impulseResponse.Ricker.dat", prog_dir);
    //sprintf(fn_opw3d, "%s/data/opw3d_nline_500_ncdp_500_ntau_4000_multi_impulseResponse.trueAmp.dat", prog_dir);
    printf("Writing the multiple impulse-response OPW data.\n");
    //write_3d_float_wb(opwData3d, Nline, Ncdp, Ntau, fn_opw3d);
    //return 0 ;

    // constant velocity.
    float   vconst  = 2500.0;
    for(int itau=0; itau<Nvtau; ++itau)
        for(int ivy=0; ivy<Nvmy; ++ivy)
            for(int ivx=0; ivx<Nvmx; ++ivx)
                vintt3d[itau][ivy][ivx] = vconst;

    char fn_vintt3d[LengthMax]="";
    printf("Writing the constant-velocity.\n");
    sprintf(fn_vintt3d, "%s/velocity/vconst3d_ntau_4000_nvy_800_nvx_800.dat", prog_dir);
    //write_3d_float_wb(vintt3d, Nvtau, Nvmy, Nvmx, fn_vintt3d);
    //return 0 ;

    // the theoretical tau-p domain ellipse curve.
    int ntau    = Ntau;
    int nphr    = 61;
    float   phrmin  = 0.;
    float   dphr    = 10/1E6;
    float   **taup2d    = alloc2float(ntau, nphr);
    memset( (void*)&taup2d[0][0], 0, sizeof(float)*1L*ntau*nphr );
    for(int iref=1; iref<nref; ++iref)
    {
        float   vrms    = vconst;
        float   v2      = vrms*vrms;
        float   t0      = dtref*iref;
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
    char fn_taup2d[LengthMax]="";
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
