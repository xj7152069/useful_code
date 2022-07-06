#include "fbopw3d.h"

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
    char prog_dir[FILE_NAME_MAX_LENGTH]="/ssd/ssd_4T/fb/q1j/new_cigs_PhaseShift";
    char fn_stack[FILE_NAME_MAX_LENGTH]="";
    sprintf(fn_stack, "%s/stackedRawImage.dat", prog_dir);
    int nline       = 500;
    int nline_first = 151;
    int nline_final = 650;
    int nline_step  = 1;
    int ncdp        = 500;
    int ncdp_first  = 151;
    int ncdp_final  = 650;
    int ncdp_step   = 1;
    float dcdp      = 20;  //m.
    float dline     = 20;  //m.
// for Q1J 3-D data
    nline       = 210;
    nline_first = 1;
    nline_final = 210;
    nline_step  = 1;
    ncdp        = 751;
    ncdp_first  = 1;
    ncdp_final  = 751;
    ncdp_step   = 1;
    dcdp      = 12.5;  //m.
    dline     = 12.5;  //m.

    int nphr, ntau, itaus, itaue;
    int iphr1   = 11;
    int iphr2   = 50;
    float   dphr, phrmin, phrmax;
    float   dt_ms   = 1.0;
    //nphr    = 101;
    nphr    = 75;
    dphr    = 10.0/1E6;
    phrmin  = 0;
    phrmax  = phrmin + dphr*(nphr-1);
    ntau    = 2001;
    int     cdp_skip    = 10;
    int     line_skip   = 10;
    cdp_skip    = 1;
    line_skip   = 1;
    int     ncdp_out    = (int)(ncdp/cdp_skip);
    int     nline_out   = (int)(nline/line_skip);
    printf("ntau=%d ncdp_out=%d nline_out=%d nphr=%d\n", ntau, ncdp_out, nline_out, nphr);

    float   ****cig4d    = alloc4float(ntau, nphr, ncdp_out, nline_out);
    memset( (void*)&cig4d[0][0][0][0], 0, sizeof(float)*1L*ntau*nphr*ncdp_out*nline_out);
    float   Byte_GB = 1.0*1024*1024*1024;
    printf("Extracting the CIGs.\n");
    printf("Size of CIG file is %f(GB).\n", sizeof(float)*1L*ntau*nphr*ncdp_out*nline_out/Byte_GB);
    float   ***image  = alloc3float(ntau, ncdp, nline);
    float   ***stack  = alloc3float(ntau, ncdp, nline);
    memset( (void*)&image[0][0][0], 0, sizeof(float)*1L*ntau*ncdp*nline);

    iphr1   = 1;
    iphr2   = 70;
    for(int iphr = iphr1; iphr <= iphr2; ++iphr)
    {
        //if ( iphr > 5 && iphr < 12 ) continue;
        char fn_commPhrImage[FILE_NAME_MAX_LENGTH]="";
        //sprintf(fn_commPhrImage, "%s/opwMig3D_PhaseShift_Iphr_%d.dat", prog_dir, iphr);
        sprintf(fn_commPhrImage, "%s/opwMig3D_PhaseShift_Iphr_%d.dat_hfp", prog_dir, iphr);
        printf("Now, reading file: %s\n", fn_commPhrImage);
        read_3d_float_rb(image, nline, ncdp, ntau, fn_commPhrImage);

        for(int iy=0; iy<nline_out; ++iy)
            for(int ix=0; ix<ncdp_out; ++ix)
                for(int it=0; it<ntau; ++it)
                {
                    int iline = iy*line_skip;
                    int icdp = ix*cdp_skip;
                    if( iline < nline && icdp < ncdp)
                        cig4d[iy][ix][iphr-iphr1][it] = image[iline][icdp][it];
                }
        for(int iy=0; iy<nline; ++iy)
            for(int ix=0; ix<ncdp; ++ix)
                for(int it=0; it<ntau; ++it)
                    stack[iy][ix][it]   += image[iy][ix][it];
    }
    printf("Writing the 4-D sparse CIGs. ntau=%d ncdp_out=%d nline_out=%d\n", ntau, ncdp_out, nline_out);
    char fn_cig[FILE_NAME_MAX_LENGTH]="";
    sprintf(fn_cig, "%s/cigs_dline%d_dcdp%d_nline%d_ncdp%d.dat", prog_dir, line_skip, cdp_skip, nline_out, ncdp_out);
    write_1d_float_wb_long ( &cig4d[0][0][0][0], 1L*ntau*nphr*ncdp_out*nline_out, fn_cig);

    printf("Writing the stacked image.\n");
    write_1d_float_wb_long ( &stack[0][0][0], 1L*ntau*ncdp*nline, fn_stack);

    return 0 ;
}
