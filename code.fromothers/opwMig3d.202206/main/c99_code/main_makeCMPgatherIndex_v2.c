/*=====================================================================*
 *                                                                      *
 *                 Make index for multiple CMP Files.                   *
 *                                                                      *
 *======================================================================*
 *   Summary:                                                           *
 *      The program will scan CMP files and creates index files.        *
 *      Format of binary index files:                                   *
 *      typedef struct{                                                 *
 *         int   line;                                                  *
 *         int   cdp;                                                   *
 *         int   sx;                                                    *
 *         int   sy;                                                    *
 *         int   gx;                                                    *
 *         int   gy;                                                    *
 *         int   offset;                                                *
 *      }cmpFileIndex;                                                  *
 *======================================================================*
 *       Author: Feng Bo                                                *
 *       Date  : 2022/01/25.                                            *
 *======================================================================*/

//include
#include "fbopw3d.h"

int main( int argc , char *argv[] )
{
    int mystat  = 0;
    int nline_first, nline_final, nline_step;
    int ncdp_first,  ncdp_final,  ncdp_step;
    int cdpMin, cdpMax, lineMin, lineMax, nfoldMax, ntr;
    float   Offset_X_Max, Offset_Y_Max, Offset_Max, Offset_Min;

    char    parFile[FILE_NAME_MAX_LENGTH]="";
    strcpy(parFile, argv[1]);

    // the maximum size of single CMP file.
    long    fileszieMaxByte = 10L*1024*1024*1024*1024;    // 10TB.
    int cmpfileOffsetBeg = 3600;    // = 0, SU-format; =3600, SEGY-format.
    int conv    = 1;                // = 0 ; assume data is in native format.
    int endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    int ns, dt_us;                  // unit of dt_us is microsecond (1E-6*second).

    int DisPlayCMPInfo  = 1;        // =1. display the CMP file information on the screen.
    int verbose = 1;
    int nfile   = -1;
    char    **CMPFileNameList   = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    **CMPIndexNameList  = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);

    char    opwdata_dir[FILE_NAME_MAX_LENGTH]="./";
    float tmpf;
    int   tmpi;
    int   OPWD_flag = 0;            // =1, only read parameters for making CMP index files.
    readPar_MakeIndex_OPWD(0, parFile,
            &nfile, CMPFileNameList, CMPIndexNameList,
            &fileszieMaxByte, &cmpfileOffsetBeg, &conv, &endian,
            &ns, &dt_us, OPWD_flag,
            opwdata_dir, &tmpf,
            &tmpi, &tmpi, &tmpi, &tmpi,
            &tmpi, &tmpf, &tmpf, &tmpf,
            &tmpi, &tmpi, &tmpi,
            verbose);

    long    ntr_cmpfile_max = 1;
    int traceBytes   = 4*(ns+60);
    ntr_cmpfile_max = (long)( (fileszieMaxByte-cmpfileOffsetBeg) / traceBytes ) + 10;
    fprintf(stdout, " maximum trace number is %ld\n", ntr_cmpfile_max);
    //return mystat;

    cmpFileIndex    *MycmpFileIndex = calloc(ntr_cmpfile_max, sizeof(cmpFileIndex));

    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        fprintf(stdout, "\n\t(%d/%d) Scanning file: %s\n", ifile+1, nfile, CMPFileNameList[ifile]);

        // open and init. cmp-gather file.
        int ostat   = open_cmp_gather( CMPFileNameList[ifile] );
        if ( 0 != ostat )
        {
            mystat  = 11;   // open cmp gather file failed.
            goto    mainProgStop;
        }

        long    ntr_cmpfile_get = -1;
        int ntr_gather_max  = -1;
        int istat   = init_cmp_gather_nonstd(
                cmpfileOffsetBeg, conv, endian, ntr_cmpfile_max,
                &ncdp_first, &ncdp_final,
                &nline_first, &nline_final,
                &ntr_gather_max, &ntr_cmpfile_get,
                &Offset_X_Max, &Offset_Y_Max,
                &Offset_Max, &Offset_Min,
                MycmpFileIndex, verbose);
        if ( 0 != istat )
        {
            mystat  = 12;   // read cmp gather file failed.
            goto    mainProgStop;
        }

        // close the cmp-gather file.
        close_cmp_gather();

        fprintf(stdout, " Initialize file %s done.\n", CMPFileNameList[ifile]);
        nline_step  = 1;
        ncdp_step   = 1;
        ntr     = ntr_gather_max;

        if( 1 == DisPlayCMPInfo )
        {
            fprintf(stdout, " ------------------------------------\n");
            fprintf(stdout, " Current CMP File is %s. \n", CMPFileNameList[ifile]);
            fprintf(stdout, " total trace number is %ld\n", ntr_cmpfile_get);
            fprintf(stdout, " max trace number with CMP-gather is %d\n", ntr_gather_max);
            fprintf(stdout, " nline_first   = %d\n", nline_first);
            fprintf(stdout, " nline_final   = %d\n", nline_final);
            fprintf(stdout, " ncdp_first    = %d\n", ncdp_first);
            fprintf(stdout, " ncdp_final    = %d\n", ncdp_final);
            //fprintf(stdout, " offset_min(m) = %f\n", offset_min);
            //fprintf(stdout, " offset_max(m) = %f\n", offset_max);
            fprintf(stdout, " ------------------------------------\n");
        }

        if ( 0 == ifile )
        {
            cdpMin  = ncdp_first;
            cdpMax  = ncdp_final;
            lineMin = nline_first;
            lineMax = nline_final;
            nfoldMax= ntr_gather_max;
        }
        else    // update the CDP- and Line- range.
        {
            if (ncdp_first < cdpMin)    cdpMin  = ncdp_first;
            if (ncdp_final > cdpMax)    cdpMax  = ncdp_final;
            if (nline_first < lineMin)  lineMin = nline_first;
            if (nline_final > lineMax)  lineMax = nline_final;
            if (ntr_gather_max > nfoldMax)  nfoldMax= ntr_gather_max;
        }

        write_cmpFileIndex(CMPIndexNameList[ifile], ntr_cmpfile_get, MycmpFileIndex);
    }

    // test of reading the CMP-index files.
    fprintf(stdout, " \n------------------------------------\n");
    cmpFileIndex    **MycmpFileIndex_2d = alloc2cmpFileIndex(ntr_cmpfile_max, nfile);
    long *ntr_files  = calloc(nfile, sizeof(long));
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr_get;
        read_cmpFileIndex(CMPIndexNameList[ifile], &ntr_get, MycmpFileIndex_2d[ifile]);
        ntr_files[ifile]    = ntr_get;
        fprintf(stdout, " total trace number is %ld\n", ntr_get);
    }

    /*
    // test of look up (cdp/line).
    //int cdpNo=3200, lineNo=2421;  // for BN3D data.
    //int cdpNo=2800, lineNo=2210;  // for Luojia3D data.
    int cdpNo=2100, lineNo=1250;    // for Q1j3D data.
    cmpFileIndex    *MycmpFileIndex_single   = calloc(nfoldMax, sizeof(cmpFileIndex));
    int itr_single  = 0;
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr = ntr_files[ifile];
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int iline    = MycmpFileIndex_2d[ifile][itr].line;
            int icdp    = MycmpFileIndex_2d[ifile][itr].cdp;
            if( lineNo == iline && cdpNo == icdp)
            {
                MycmpFileIndex_single[itr_single].line  = iline;
                MycmpFileIndex_single[itr_single].cdp   = icdp;
                MycmpFileIndex_single[itr_single].sx   = MycmpFileIndex_2d[ifile][itr].sx;
                MycmpFileIndex_single[itr_single].sy   = MycmpFileIndex_2d[ifile][itr].sy;
                MycmpFileIndex_single[itr_single].gx   = MycmpFileIndex_2d[ifile][itr].gx;
                MycmpFileIndex_single[itr_single].gy   = MycmpFileIndex_2d[ifile][itr].gy;
                ++ itr_single;
            }
        }
    }

    // display information.
    fprintf(stdout, " \n------------------------------------\n");
    fprintf(stdout, " lineNo = %d\n", lineNo);
    fprintf(stdout, " cdpNo = %d\n", cdpNo);
    fprintf(stdout, " ntrCMP = %d\n", itr_single);
    for ( int itr = 0 ; itr < itr_single; ++itr )
    {
        int cdp  = MycmpFileIndex_single[itr].cdp;
        int line = MycmpFileIndex_single[itr].line;
        int sx  = MycmpFileIndex_single[itr].sx;
        int sy  = MycmpFileIndex_single[itr].sy;
        int gx  = MycmpFileIndex_single[itr].gx;
        int gy  = MycmpFileIndex_single[itr].gy;
        int mx  = 0.5*(sx + gx);
        int my  = 0.5*(sy + gy);
        fprintf(stderr, " line=%d cdp=%d sx=%d sy=%d gx=%d gy=%d mx=%d my=%d \n", line, cdp, sx, sy, gx, gy, mx, my);
    }

    free(MycmpFileIndex_single);
    */

    free2cmpFileIndex(MycmpFileIndex_2d);

    if( 1 == DisPlayCMPInfo )
    {
        fprintf(stdout, " ------------------------------------\n");
        fprintf(stdout, " Cheers!\n");
        fprintf(stdout, " CMP files scanned %d\n", nfile);
        fprintf(stdout, " Statistical information: \n");
        fprintf(stdout, " lineMin = %d\n", lineMin);
        fprintf(stdout, " lineMax = %d\n", lineMax);
        fprintf(stdout, " cdpMin  = %d\n", cdpMin);
        fprintf(stdout, " cdpMax  = %d\n", cdpMax);
        fprintf(stdout, " nfoldMax= %d\n", nfoldMax);
        fprintf(stdout, " Offset_X_Max = %f(m)\n", Offset_X_Max);
        fprintf(stdout, " Offset_Y_Max = %f(m)\n", Offset_Y_Max);
        fprintf(stdout, " Offset_Max = %f(m)\n", Offset_Max);
        fprintf(stdout, " Offset_Min = %f(m)\n", Offset_Min);
        fprintf(stdout, " ------------------------------------\n");
    }

mainProgStop:
    fprintf(stdout, "mystat = %d\n", mystat);

    free2char(CMPFileNameList);
    free2char(CMPIndexNameList);
    free(MycmpFileIndex);

    return mystat;
}