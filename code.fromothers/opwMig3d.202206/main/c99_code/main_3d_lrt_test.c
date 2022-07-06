//include
#include "fbopw3d.h"

//  Global variable: sx = sx/Coor_S_R_Scale; gx = gx/Coor_S_R_Scale
float    Coor_S_R_Scale;

// function prototypes begin.
int    getTraceNumber ( int ns, int *ntr, char *fn );
// function prototypes end.

int    getTraceNumber ( int ns, int *ntr, char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file: %s to read.\n" , fn ) ;
        return  -1;
    }
    long    itr = 0;
    fbsegy_nonstd   segy_hdr;
    float   *trace  = calloc(ns, FLOAT_SIZE_BYTES);
    int flag_eof = 0;
    while ( 1 )
    {
        fread( &segy_hdr, sizeof(fbsegy_nonstd) , 1 , fp ) ;
        flag_eof = feof(fp);
        if (flag_eof)
        {
            // read header failed.
            break;
        }
        fread( trace, FLOAT_SIZE_BYTES, ns , fp ) ;
        if (flag_eof)
        {
            // read trace error.
            break;
        }
        ++itr;
    }

    fclose(fp);
    free(trace);
    *ntr    = (int)itr;
    return  0;
}

int main( int argc , char *argv[] )
{
    int mystat  = 0;
    /*
    if( argc != 3 )
    {
        printf(" Wrong Parameters!\n");
        return 1 ;
    }
     */
    char    parFile[FILE_NAME_MAX_LENGTH]="";
    char    cmpfn[FILE_NAME_MAX_LENGTH]="";
    strcpy(parFile, argv[1]);
    strcpy(cmpfn,   argv[2]);

    // init. of global variables.
    Coor_S_R_Scale  = 1.0;      // set the default value of global variable.

    // the maximum size of single CMP file.
    long    fileszieMaxByte = 10L*1024*1024*1024*1024;    // 10TB.
    int cmpfileOffsetBeg = 3600;    // = 0, SU-format; =3600, SEGY-format.
    int conv    = 1;                // = 0 ; assume data is in native format.
    int endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    int ns, dt_us;                  // unit of dt_us is microsecond (1E-6*second).

    // parameters for tau-phr transform.
    int ntrMin_PWD   = 21;
    int nphr, ntau, itaus, itaue;
    float   dphr, phrmin, phrmax;

    // parameters for parallel-calculation.
    int nline_loc   = 10;
    int firstCMPLineNo_PWD, lastCMPLineNo_PWD;

    int DisPlayCMPInfo  = 1;    // =1. display the CMP file information on the screen.
    int verbose = 1;
    int myid    = 0;
    int nfile   = -1;

    char    **CMPFileNameList   = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    **CMPIndexNameList  = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    opwdata_dir[FILE_NAME_MAX_LENGTH]="./";
    int   OPWD_flag = 1;            // =1, only read parameters for making CMP index files.
    readPar_MakeIndex_OPWD(myid, parFile,
            &nfile, CMPFileNameList, CMPIndexNameList,
            &fileszieMaxByte, &cmpfileOffsetBeg, &conv, &endian,
            &ns, &dt_us, OPWD_flag,
            opwdata_dir, &Coor_S_R_Scale,
            &ntrMin_PWD, &ntau, &itaus, &itaue,
            &nphr, &phrmin, &phrmax, &dphr,
            &nline_loc, &firstCMPLineNo_PWD, &lastCMPLineNo_PWD,
            verbose);
    float   dt_s= dt_us/1E6;

    // tau-ph transform of a single CDP gather.
    int ntrCMP = 0;
    getTraceNumber(ns, &ntrCMP, cmpfn);
    printf(" cmp filename is: %s\t ntrCMP=%d ns=%d\n", cmpfn, ntrCMP, ns);

    // read a CMP gather.
    float   **gather    = alloc2float(ns, ntrCMP);
    float   **taup2d    = alloc2float(ns, nphr);
    fbsegy_std  *suhdr  = calloc(ntrCMP, sizeof(fbsegy_std));
    read_2d_float_rb_su(gather, suhdr, ntrCMP, ns, cmpfn);

    // calculate the offset.
    float   *offset = calloc(ntrCMP, FLOAT_SIZE_BYTES);
    float   offsetMin   = 100000.;
    float   offsetMax   = 0.;
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
    {
        float   sx  = (float)suhdr[itr].sx/Coor_S_R_Scale;
        float   sy  = (float)suhdr[itr].sy/Coor_S_R_Scale;
        float   gx  = (float)suhdr[itr].gx/Coor_S_R_Scale;
        float   gy  = (float)suhdr[itr].gy/Coor_S_R_Scale;
        offset[itr] = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) );
        if (offset[itr] <offsetMin)
            offsetMin   = offset[itr];
        if (offset[itr] >offsetMax)
            offsetMax   = offset[itr];
    }
    fprintf(stdout, "offsetMin=%f offsetMax=%f\n", offsetMin, offsetMax);

    // mute the refraction signals.
    float   vel_surf    = 1800.0;
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
    {
        float   tsr = offset[itr]/vel_surf;
        int ntmute   = (int)(tsr/dt_s);
        if ( ntmute >= ns )
            continue;
        //for ( int it = ntmute ; it < ns; ++it )
        for ( int it = 0 ; it < ntmute; ++it )
            gather[itr][it] = 0;
    }
    write_2d_float_wb(gather, ntrCMP, ns, "gather_muteRefraction.dat");

    // type 1: calculate the local-LRT of the input CMP gather.
    // i.e., separating the CMP gather by local spatial windows,
    // and do LRT to each local gather, and sum up the tau-p spectrum.
    float   offsetWidth  = 400.0;
    offsetWidth  = 800.0;
    int	reconstrutDataFlag = 0;
    verbose = 1;
    float   **gather_out   = alloc2float(ns, ntrCMP);
    float   **taup2d_lslrt    = alloc2float(ns, nphr);
    memset( (void*)&taup2d_lslrt[0][0], 0, FLOAT_SIZE_BYTES*1L*ntau*nphr);
    fprintf(stdout, "calling apply_LSLRT_to_CMPgather()\n");
    apply_LocalLRT_to_CMPgather(ntrCMP, offset, ns, dt_s, gather, offsetWidth,
            nphr, phrmin, dphr, ntau, taup2d_lslrt, verbose);
    write_2d_float_wb(taup2d_lslrt, nphr, ntau, "taup2d_localLRT.dat");

    // type 2: perform the adaptive finite-aperture tau-phr transform for a CMP gather.
    // i.e., apply the global LRT to the CMP gather by adaptively
    // selecting the spatial integration aperture, see Feng et al., 2009.
    float   dt_ms  = dt_us/1E3;
    float   *vrms1d = calloc(ns, FLOAT_SIZE_BYTES);
    float   vstk_min    = 1500.;
    float   vstk_max    = 4500.;
    for ( int it = 0 ; it < ns; ++it )
        vrms1d[it]  = vstk_min + (vstk_max-vstk_min)*it/ns;
    float   **taup2d_frss    = alloc2float(ns, nphr);
    memset( (void*)&taup2d_frss[0][0], 0, FLOAT_SIZE_BYTES*1L*ntau*nphr);
    fprintf(stdout, "calling tauph_transform_3d_to_2d_robust()\n");
    tauph_transform_3d_to_2d_robust(ns, ntau, nphr,
            ntrCMP, itaus, itaue,
            gather, taup2d_frss,
            vrms1d, offset, dt_ms, phrmin, dphr);
    write_2d_float_wb(taup2d_frss, nphr, ntau, "taup2d_frss.dat");

    postProcessingTauPspectrum(taup2d_frss, nphr, ntau, taup2d_lslrt);
    write_2d_float_wb(taup2d_lslrt, nphr, ntau, "taup2d_localLRT_pp_noOMP.dat");

    // type 3: calculate the global-LRT of the input CMP gather.
    // i.e., apply the classical LRT to the CMP gather.
    float   x0  = 0;
    float   **gather1   = alloc2float(ntau, ntrCMP);
    float   *halfoffset = calloc(ntrCMP, FLOAT_SIZE_BYTES);
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
    {
        halfoffset[itr]  = 0.5*offset[itr];
        for (int it = 0 ; it < ntau; ++it)
            gather1[itr][it]    = gather[itr][it];
    }
    fprintf(stdout, "calling slantStack_TD_2D_omp()\n");
    memset( (void*)&taup2d[0][0], 0, sizeof(float)*nphr*ntau );
    slantStack_TD_2D_omp(gather1, ntrCMP, ntau, dt_s, halfoffset, x0,
            taup2d, nphr, phrmin, dphr);
    write_2d_float_wb(taup2d, nphr, ntau, "taup2d_globalLRT.dat");

    // type 4: For XiangJian, calculate the time-domain Least-squares LRT code.
    // i.e., separating the CMP gather by local spatial windows,
    // and do least-squares LRT with sparsity promotion to each local gather, and sum up the tau-p spectrum.
    /*
    apply_LSLRT_to_CMPgather(ntrCMP, offset, suhdr, ns, dt_s, gather, offsetWidth,
           nphr, phrmin, dphr, ntau, taup2d_lslrt, reconstrutDataFlag, gather_out, verbose);
    */

    // free memory.
    free2char(CMPFileNameList);
    free2char(CMPIndexNameList);

    return mystat;
}
