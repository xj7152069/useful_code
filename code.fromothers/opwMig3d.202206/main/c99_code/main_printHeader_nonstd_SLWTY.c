
//include
#include "fbopw3d.h"

int printHeader_nonstd(
// Input parameters.
    char    *cmpGatherFilename,
    int fileOffsetBeg,      // =0, SU-format; =3600, SEGY-format.
    int conv,               // convert data to native format.
                            // = 0 ; assume data is in native format.
    int endian,             // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max,
    int verbose             // =1, print debug information.
    )
{
    int stat    = 0;
    FILE *cmpGatherFilefp = NULL;   // file pointer for operating the cmp-gather file.
    // open the cmp gather file.
    if ( NULL == (cmpGatherFilefp = fopen(cmpGatherFilename,"rb")) )
    {
        // open file failed.
        stat    = +1;
    }

    unsigned short  ns, dtus;
    int lineno, cdpno, fldr, ep, cdp;
    int sx, sy, gx, gy;
    int mxi, myi;
    int kw181;
    int kw185;
    int kw189;
    int kw193;
    int kw197;
    int kw201;

    // 1. read the segy-header of the first trace and get some parameters.
    fbsegy_nonstd   segy_hdr    = {0} ;
    fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
    fread( &segy_hdr, sizeof(fbsegy_nonstd) , 1 , cmpGatherFilefp );
    if ( 0 == endian && 1 == conv ) // swapbytes of header keywords.
    {
        printf("conv=%d endian=%d \n", conv, endian);
        swap_u_short_2(&segy_hdr.ns);
        swap_u_short_2(&segy_hdr.dt);
        swap_int_4(&segy_hdr.fldr);
        swap_int_4(&segy_hdr.ep);
        swap_int_4(&segy_hdr.cdp);
        swap_int_4(&segy_hdr.lineno);
        swap_int_4(&segy_hdr.cdpno);
        swap_int_4(&segy_hdr.byte189);
        swap_int_4(&segy_hdr.byte193);
        swap_int_4(&segy_hdr.byte197);
        swap_int_4(&segy_hdr.byte201);
    }

    ns      = segy_hdr.ns;
    dtus    = segy_hdr.dt;
    fldr    = segy_hdr.fldr;
    ep      = segy_hdr.ep;
    cdp     = segy_hdr.cdp;
    lineno  = segy_hdr.lineno;
    cdpno   = segy_hdr.cdpno;
    kw189   = segy_hdr.byte189;
    kw193   = segy_hdr.byte193;
    kw197   = segy_hdr.byte197;
    kw201   = segy_hdr.byte201;
    if ( 1 == verbose )
    {
        printf("ns=%d dtus=%d lineno=%d cdpno=%d fldr=%d ep=%d cdp=%d kw189=%d kw193=%d\n",
            ns, dtus, lineno, cdpno, fldr, ep, cdp, kw189, kw193);
    }
    //exit(1);

    int line1   = lineno;
    int cdp1    = cdpno;
    int MAX_INT = (unsigned)(-1)>>1;
    int cdpmin=MAX_INT, cdpmax=-1;
    int linemin=MAX_INT, linemax=-1;
    int ntrmin  = MAX_INT;
    int ntrmax  = -1;
    float   offX_Max = 0;
    float   offY_Max = 0;
    float   off_Max = 0;
    float   off_Min =1E10;
    if ( 1 == verbose )
printf("ns=%d dt=%d(us) line1=%d cdp1=%d MAX_INT=%d\n", ns, dtus, line1, cdp1, MAX_INT);

    // 2. read the segy-header of each trace.
    int flag_eof = 0;
    fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
    long    itr = 0;
    int ntrace_cmp  = 0;
    while ( 1 )
    {
        fread( &segy_hdr, sizeof(fbsegy_nonstd) , 1 , cmpGatherFilefp ) ;
        flag_eof = feof(cmpGatherFilefp);
        if (flag_eof)
        {
            // read header failed.
            break;
        }
        //if( 0 == itr%10000 )  printf("itr=%ld\n", itr);

        if ( 0 == endian && 1 == conv ) // swapbytes of header keywords.
        {
            swap_u_short_2(&segy_hdr.ns);
            swap_u_short_2(&segy_hdr.dt);
            swap_int_4(&segy_hdr.fldr);
            swap_int_4(&segy_hdr.ep);
            swap_int_4(&segy_hdr.cdp);
            swap_int_4(&segy_hdr.lineno);
            swap_int_4(&segy_hdr.cdpno);
            swap_int_4(&segy_hdr.sx);
            swap_int_4(&segy_hdr.sy);
            swap_int_4(&segy_hdr.gx);
            swap_int_4(&segy_hdr.gy);
            swap_int_4(&segy_hdr.byte189);
            swap_int_4(&segy_hdr.byte193);
            swap_int_4(&segy_hdr.byte197);
            swap_int_4(&segy_hdr.byte201);
        }

        ns  = segy_hdr.ns;
        dtus    = segy_hdr.dt;
        fldr    = segy_hdr.fldr;
        ep      = segy_hdr.ep;
        cdp     = segy_hdr.cdp;
        lineno  = segy_hdr.lineno;
        cdpno   = segy_hdr.cdpno;
        sx  = segy_hdr.sx;
        sy  = segy_hdr.sy;
        gx  = segy_hdr.gx;
        gy  = segy_hdr.gy;
        mxi = segy_hdr.byte189;
        myi = segy_hdr.byte193;
//        if(1292==lineno)
        printf("ns=%d dtus=%d mx=%f mxi=%d my=%f myi=%d 181-184(int)=%d 185-188(int)=%d 189-192(int)=%d 193-197(int)=%d fldr=%d ep=%d cdp=%d\n",
            ns, dtus, 0.5*(sx+gx), mxi, 0.5*(sy+gy), myi, lineno, cdpno, segy_hdr.byte189, segy_hdr.byte193, fldr, ep, cdp);
        float   offset_x    = fabs(gx - sx);
        float   offset_y    = fabs(gy - sy);
        float   offset_xy   = sqrt(offset_x*offset_x + offset_y*offset_y);
        //printf("itr=%ld fldr=%d ep=%d cdp=%d lineno=%d cdpno=%d sx=%d sy=%d gx=%d gy=%d offset=%f\n", itr, fldr, ep, cdp, lineno, cdpno, sx, sy, gx, gy, offset_xy);

        // Skip a trace.
        fseek( cmpGatherFilefp, (long)ns*4L, SEEK_CUR );

        ++ itr;
        ++ ntrace_cmp;

        /*
         */
        if ( itr > ntr_cmpfile_max )
        {
// trace number exceeds the maximum trace number defined in the main program.
fprintf(stderr, "Error! Trace number exceeds the maximum trace number defined in the main program.");
fprintf(stderr, "ntr_cmpfile_max = %ld\n", ntr_cmpfile_max);
fprintf(stderr, "Trace number scanned: %ld\n", itr);
            stat    = -1;
            break;
        }
    }
printf("ntrmax=%d\n", ntrmax);
    fclose(cmpGatherFilefp);

    return stat;
}

int main( int argc , char *argv[] )
{
    int mystat  = 0;

    char    segyFile[FILE_NAME_MAX_LENGTH]="";
    strcpy(segyFile, argv[1]);
    int cmpfileOffsetBeg = 3600;    // = 0, SU-format; =3600, SEGY-format.
    int conv    = 1;                // = 0 ; assume data is in native format.
    int endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max = 12641L;

    ntr_cmpfile_max = 20;
    printHeader_nonstd(segyFile, cmpfileOffsetBeg, conv, endian, ntr_cmpfile_max, 1);

    return  mystat;
}

