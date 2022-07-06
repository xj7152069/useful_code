
//include
#include "fbCommon.h"
#include "fbswapbyte.h"
#include "fbsegy_nonstd_NWGI.h"

int printHeader_nonstd(
// Input parameters.
    char    *cmpGatherFilename,
    int fileOffsetBeg,      // =0, SU-format; =3600, SEGY-format.
    int conv,               // convert data to native format.
                            // = 0 ; assume data is in native format.
    int endian,             // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max,
    int outputAcq3D_flag,   // =1, output the acquisition file for tomo.
    char    *acqFilename,   // the acquisition filename.
    int outputshot_flag,    // =1, output the shot map.
    char    *shotmapFilename,// the shot map filename.
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
		fprintf(stderr, "Can not open file: %s to read.\n" , cmpGatherFilename);
		return stat ;
    }

	FILE	*acqFilefp = NULL ;
    if ( 1 == outputAcq3D_flag )
    {
        if( NULL == (acqFilefp=fopen(acqFilename,"w")) )
        {
            fprintf(stderr, "Can not open file: %s to write.\n" , acqFilename);
            stat    = -1;
            return stat ;
        }
    }

	FILE	*shatmapFilefp = NULL ;
    if ( 1 == outputshot_flag )
    {
        if( NULL == (shatmapFilefp=fopen(shotmapFilename,"w")) )
        {
            fprintf(stderr, "Can not open file: %s to write.\n" , shotmapFilename);
            stat    = -1;
            return stat ;
        }
    }

    unsigned short  ns, dtus;
    int lineno, cdpno, fldr, ep, cdp;
    int sx, sy, sz, gx, gy, gz, tfb;
    int kw205;  // first-break;
    int kw237;  // trace number within shotgather.
    int shotno_prev, shotno_curr;
    int shot_num    = 0;
    int trace_num    = 0;
    float   offset_min  = 50.0; //m.
    float   tfb_min     = 50.0; //ms.

    // 1. read the segy-header of the first trace and get some parameters.
    fbsegy_nonstd_NWGI   segy_hdr    = {0} ;
    fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
    fread( &segy_hdr, sizeof(fbsegy_nonstd_NWGI) , 1 , cmpGatherFilefp );
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
        swap_int_4(&segy_hdr.byte205);
        swap_int_4(&segy_hdr.byte237);
    }

    ns      = segy_hdr.ns;
    dtus    = segy_hdr.dt;
    fldr    = segy_hdr.fldr;
    ep      = segy_hdr.ep;
    cdp     = segy_hdr.cdp;
    lineno  = segy_hdr.lineno;
    cdpno   = segy_hdr.cdpno;
    kw205   = segy_hdr.byte205;
    kw237   = segy_hdr.byte237;
    shotno_prev = fldr-1;
    if ( 1 == verbose )
    {
        printf("ns=%d dtus=%d fldr=%d tfb=%d traceid=%d\n", ns, dtus, fldr, kw205, kw237);
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
        fread( &segy_hdr, sizeof(fbsegy_nonstd_NWGI) , 1 , cmpGatherFilefp ) ;
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
            swap_int_4(&segy_hdr.selev);
            swap_int_4(&segy_hdr.gx);
            swap_int_4(&segy_hdr.gy);
            swap_int_4(&segy_hdr.gelev);
            swap_int_4(&segy_hdr.byte205);
            swap_int_4(&segy_hdr.byte237);
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
        sz  = segy_hdr.selev;
        gx  = segy_hdr.gx;
        gy  = segy_hdr.gy;
        gz  = segy_hdr.gelev;
        kw237   = segy_hdr.byte237;
        kw205   = segy_hdr.byte205;
        tfb     = kw205;
        float   offset_x    = fabs(gx - sx);
        float   offset_y    = fabs(gy - sy);
        float   offset_z    = fabs(gz - sz);
        float   offset_xy   = sqrt(offset_x*offset_x + offset_y*offset_y);
        float   offset      = sqrt(offset_x*offset_x + offset_y*offset_y + offset_z*offset_z);
        /*
        printf("itr=%ld fldr=%d sx=%d sy=%d sz=%d gx=%d gy=%d gz=%d offset=%f traceid=%d tfb=%d(ms)\n",
                itr+1, fldr, sx, sy, sz, gx, gy, gz, offset_xy, kw237, tfb);
        */

        // output the acquisition file for tomo.
        if ( 1 == outputAcq3D_flag )
        {
            shotno_curr = fldr;
            float	coordinate_z, coordinate_x, coordinate_y;
            float	tmp1, tmp2;
            int     srflag		= 0;
            if ( shotno_curr != shotno_prev)
            {   // the next shotgather.
                ++shot_num;
                srflag		= 0;
                coordinate_z= (float)sz;
                coordinate_x= (float)sx;
                coordinate_y= (float)sy;
                tmp1        = 0;
                tmp2        = 0;
                fprintf(acqFilefp, "%f\t%f\t%f\t%f\t%f\t%d\n",
                    coordinate_z, coordinate_x, coordinate_y, tmp1, tmp2, srflag);
                shotno_prev = shotno_curr;

                printf("fldr=%d sx=%d sy=%d sz=%d gx=%d gy=%d gz=%d offset=%f traceid=%d tfb=%d(ms)\n",
                        fldr, sx, sy, sz, gx, gy, gz, offset, kw237, tfb);

                if ( 1 == outputshot_flag )
                {
                    fprintf(shatmapFilefp, "%d\t%d\n", sx, sy);
                }
            }

            srflag		= 1;
            coordinate_z= (float)gz;
            coordinate_x= (float)gx;
            coordinate_y= (float)gy;
            tmp1        = tfb/1000.0;   // unit of first-arrival is "second".
            tmp2        = tfb/1000.0;   // unit of first-arrival is "second".
            //if ( tfb > 0 )
            if ( tfb > tfb_min && offset > offset_min )
            {
                fprintf(acqFilefp, "%f\t%f\t%f\t%f\t%f\t%d\n",
                        coordinate_z, coordinate_x, coordinate_y, tmp1, tmp2, srflag);
                ++trace_num;
            }
        }

        // Skip a trace.
        fseek( cmpGatherFilefp, (long)ns*4L, SEEK_CUR );

        ++ itr;
        ++ ntrace_cmp;

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

    if ( 1 == outputAcq3D_flag )
    {
        fprintf(stderr, "Total Trace number scanned: %ld\n", itr);
        fprintf(stderr, "Shot number scanned:  %d\n", shot_num);
        fprintf(stderr, "Valid trace number scanned:  %d\n", trace_num);
        fclose(acqFilefp);
    }

    return stat;
}

int main( int argc , char *argv[] )
{
    int mystat  = 0;

    char    segyFilename[FILE_NAME_MAX_LENGTH]="";
    char    acqFilename[FILE_NAME_MAX_LENGTH]="";
    char    shotmapFilename[FILE_NAME_MAX_LENGTH]="";
    strcpy(segyFilename,    argv[1]);
    strcpy(acqFilename,     argv[2]);
    strcpy(shotmapFilename, argv[3]);
    int cmpfileOffsetBeg = 3600;    // = 0, SU-format; =3600, SEGY-format.
    int conv    = 1;                // = 0 ; assume data is in native format.
    int endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max = 4210087;
    int outputAcq3D_flag    = 1;    // =1, output the acquisition file for tomo.
    int outputshot_flag     = 1;    // =1, output the shot map.

    //ntr_cmpfile_max = 2000;
    /*
    printHeader_nonstd(segyFilename, cmpfileOffsetBeg, conv, endian, ntr_cmpfile_max,
        outputAcq3D_flag, acqFilename, outputshot_flag, shotmapFilename, 1);
    */

    int rotateCoorFlag     = 1;    // =1, rotate the coordinates.
    float   angle   = 20.742;
	FILE	*shatmapFilefp_r = NULL ;
	FILE	*shatmapFilefp_w = NULL ;
    char    shotmapFilename_rot[FILE_NAME_MAX_LENGTH]="./shotmap.rot.txt";
    if ( 1 == rotateCoorFlag )
    {
        if( NULL == (shatmapFilefp_r=fopen(shotmapFilename,"r")) )
        {
            fprintf(stderr, "Can not open file: %s to read.\n" , shotmapFilename);
            mystat    = -1;
        }

        if( NULL == (shatmapFilefp_w=fopen(shotmapFilename_rot,"w")) )
        {
            fprintf(stderr, "Can not open file: %s to write.\n" , shotmapFilename_rot);
            mystat    = -1;
        }
        float   cos_sita    = cos(angle*PI/180.0);
        float   sin_sita    = sin(angle*PI/180.0);

        while ( 1 )
        {
            float   sx, sy;
            fscanf(shatmapFilefp_r, "%f%f\n", &sx, &sy);
            int flag_eof = feof(shatmapFilefp_r);
            if (!flag_eof)
            {
                float   sx1 = cos_sita*sx - sin_sita*sy;
                float   sy1 = sin_sita*sx + cos_sita*sy;
                fprintf(shatmapFilefp_w, "%f\t%f\n", sx1, sy1);
            }
            else
            {
                // read header failed.
                break;
            }
        }
    }
    fclose(shatmapFilefp_r);
    fclose(shatmapFilefp_w);


    return  mystat;
}
