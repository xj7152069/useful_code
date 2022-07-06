//include
#include "fbopw3d.h"

int getMinMaxLineNo(
        int nfile,
        cmpFileIndex    **MycmpFileIndex_2d,
        long *ntr_files,
        int *LineMin_get,
        int *LineMax_get,
        int *NfoldMin_get,
        int *NfoldMax_get,
        int verbose
        )
{
    int mystat  = 0;

    int lineMin = 100000000;
    int lineMax = -1;
    int ntrMin = 100000000;
    int ntrMax = -1;
    int ntr_cmp = 0;
    int line1   = MycmpFileIndex_2d[0][0].line;
    int cdp1    = MycmpFileIndex_2d[0][0].cdp;
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr = ntr_files[ifile];
        ntr_cmp = 0;
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int ilineNo   = MycmpFileIndex_2d[ifile][itr].line;
            if (ilineNo < lineMin)    lineMin  = ilineNo;
            if (ilineNo > lineMax)    lineMax  = ilineNo;

            int icdpNo   = MycmpFileIndex_2d[ifile][itr].cdp;
            if ( icdpNo != cdp1 )  // a new CMP gather.
            {
                if (ntr_cmp < ntrMin)    ntrMin  = ntr_cmp;
                if (ntr_cmp > ntrMax)    ntrMax  = ntr_cmp;
                cdp1    = icdpNo;
                ntr_cmp  = 0;
                if ( ilineNo != line1 )    // a new CMP Line.
                {
                    line1   = ilineNo;
                    cdp1    = icdpNo;
                }
            }
            ++ ntr_cmp;
        }
    }

    if ( 1 == verbose )
    {
        fprintf(stdout, " lineMin = %d\n", lineMin);
        fprintf(stdout, " lineMax = %d\n", lineMax);
        fprintf(stdout, " nfoldMin= %d\n", ntrMin);
        fprintf(stdout, " nfoldMax= %d\n", ntrMax);
    }

    *NfoldMin_get   = ntrMin;
    *NfoldMax_get   = ntrMax;
    *LineMin_get    = lineMin;
    *LineMax_get    = lineMax;

    return mystat;
}

void getMinMaxcdpNo(
        int nfile,
        long *ntr_files,
        int nline,
        CMPLineInfo *MyCMPLineInfo,
        cmpFileIndex    **MycmpFileIndex_2d,
        int lineMin,
        int lineMax,
        int *cdpMin_proj_get,
        int *cdpMax_proj_get,
        int verbose
        )
{
    int cdpMin_proj = MycmpFileIndex_2d[0][0].cdp;
    int cdpMax_proj = MycmpFileIndex_2d[0][0].cdp;

#pragma omp parallel for shared(cdpMin_proj, cdpMax_proj)
    for ( int lineNo = lineMin; lineNo <= lineMax; ++ lineNo )
    {
        int cdpMin_loc  = 100000000;
        int cdpMax_loc  = -1;
        int ntrace_line = 0;
        for ( int ifile = 0 ; ifile < nfile; ++ifile )
        {
            long    ntr = ntr_files[ifile];
            for ( long itr = 0 ; itr < ntr; ++itr )
            {
                int lineNo_now  = MycmpFileIndex_2d[ifile][itr].line;
                int icdpNo      = MycmpFileIndex_2d[ifile][itr].cdp;
                if( lineNo == lineNo_now )
                {
                    if (icdpNo < cdpMin_loc)    cdpMin_loc  = icdpNo;
                    if (icdpNo > cdpMax_loc)    cdpMax_loc  = icdpNo;
                    ++ ntrace_line;
                }
            }
        }

        int iline   = lineNo - lineMin;
        MyCMPLineInfo[iline].lineNo         = lineNo;
        MyCMPLineInfo[iline].ncdp_line      = cdpMax_loc - cdpMin_loc + 1;
        MyCMPLineInfo[iline].ntrace_line    = ntrace_line;
        MyCMPLineInfo[iline].cdpMin_line    = cdpMin_loc;
        MyCMPLineInfo[iline].cdpMax_line    = cdpMax_loc;

        //fprintf(stdout, " lineNo=%d cdpMin=%d cdpMax=%d\n", iline, cdpMin_loc, cdpMax_loc);
        // update the global minimum and maximum cdp number in the project.
#pragma omp critical (inFunc_getMinMaxcdpNo)
        {
            if (cdpMin_loc < cdpMin_proj)    cdpMin_proj  = cdpMin_loc;
            if (cdpMax_loc > cdpMax_proj)    cdpMax_proj  = cdpMax_loc;
        }
    }

    // display project information.
    if ( 1 == verbose )
    {
        for ( int lineNo = lineMin; lineNo <= lineMax; ++ lineNo )
        {
            int iline   = lineNo - lineMin;
            int cdpMin_loc  = MyCMPLineInfo[iline].cdpMin_line;
            int cdpMax_loc  = MyCMPLineInfo[iline].cdpMax_line;
            int lineNo      = MyCMPLineInfo[iline].lineNo;
            int ncdp_line   = MyCMPLineInfo[iline].ncdp_line;
            int ntrace_line = MyCMPLineInfo[iline].ntrace_line;
            fprintf(stdout, " lineNo=%d ncdp=%d ntrace=%d cdpMin=%d cdpMax=%d\n",
                    lineNo, ncdp_line, ntrace_line, cdpMin_loc, cdpMax_loc);
        }
    }

    *cdpMin_proj_get    = cdpMin_proj;
    *cdpMax_proj_get    = cdpMax_proj;

}

int main( int argc , char *argv[] )
{
    /*
    int mystat  = 0;
    char    parFile[FILE_NAME_MAX_LENGTH]="";
    strcpy(parFile, argv[1]);

    int DisPlayCMPInfo  = 1;    // =1. display the CMP file information on the screen.
    int verbose = 1;
    int nfile   = -1;
    char    **CMPFileNameList   = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    **CMPIndexNameList  = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    readpar(parFile, &nfile, CMPFileNameList, CMPIndexNameList, verbose);

    // for Luojia3D data.
    int cmpfileOffsetBeg= 0; // =0, SU-format; =3600, SEGY-format.
    int conv    = 0;        // = 0 ; assume data is in native format.
    int endian  = 0;        // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max = 100000000;
    int ns  = 2901;
    int dt  = 2000;
     */
    int mystat  = 0;
    //int nline_first, nline_final, nline_step;
    //int ncdp_first,  ncdp_final,  ncdp_step;
    //float   Offset_X_Max, Offset_Y_Max, Offset_Max, Offset_Min;
    //int cdpMin, cdpMax;
    int lineMin, lineMax, nfoldMin, nfoldMax, ntr;

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

    char    opwdata_dir[FILE_NAME_MAX_LENGTH]=".";
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

    cmpFileIndex    **MycmpFileIndex_2d = alloc2cmpFileIndex(ntr_cmpfile_max, nfile);
    long *ntr_files  = calloc(nfile, sizeof(long));

    // read multiple CMP-index files.
    fprintf(stdout, " \n------------------------------------\n");
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr_cmpfile;
        read_cmpFileIndex(CMPIndexNameList[ifile], &ntr_cmpfile, MycmpFileIndex_2d[ifile]);
        ntr_files[ifile]    = ntr_cmpfile;
        fprintf(stdout, " total trace number is %ld\n", ntr_cmpfile);
    }

    // get the minimum and maximum line number in the project.
    int nline;
    getMinMaxLineNo(
            nfile, MycmpFileIndex_2d, ntr_files,
            &lineMin, &lineMax, &nfoldMin, &nfoldMax,
            verbose);
    nline   = lineMax - lineMin + 1;

    // get the minimum and maximum cdp number of each CMP line.
    CMPLineInfo *MyCMPLineInfo  = calloc(nline, sizeof(CMPLineInfo));
    int ncdp_proj, cdpMin_proj, cdpMax_proj;
    getMinMaxcdpNo(
            nfile, ntr_files, nline, MyCMPLineInfo,
            MycmpFileIndex_2d, lineMin, lineMax,
            &cdpMin_proj, &cdpMax_proj, verbose);
    ncdp_proj    = cdpMax_proj - cdpMin_proj + 1;

    ProjCMPInfo MyProjCMPInfo;
    MyProjCMPInfo.lineMin   = lineMin;
    MyProjCMPInfo.lineMax   = lineMax;
    MyProjCMPInfo.cdpMin_proj   = cdpMin_proj;
    MyProjCMPInfo.cdpMax_proj   = cdpMax_proj;

    if( 1 == DisPlayCMPInfo )
    {
        fprintf(stdout, "\n\t ----- display project information. -----\n");
        fprintf(stdout, " ------------------------------------\n");
        fprintf(stdout, " Cheers!\n");
        fprintf(stdout, " CMP-index files scanned: %d\n", nfile);
        fprintf(stdout, " Statistical information: \n");
        fprintf(stdout, " nline     = %d\n", nline);
        fprintf(stdout, " lineMin   = %d\n", lineMin);
        fprintf(stdout, " lineMax   = %d\n", lineMax);
        fprintf(stdout, " ncdp_proj = %d\n", ncdp_proj);
        fprintf(stdout, " cdpMin    = %d\n", cdpMin_proj);
        fprintf(stdout, " cdpMax    = %d\n", cdpMax_proj);
        fprintf(stdout, " nfoldMax  = %d\n", nfoldMax);
        fprintf(stdout, " ------------------------------------\n");
    }
    fflush(stdout);
    fflush(stderr);

    // find a CMP-gather using a given cdp+line.
    // for Q1j3D data: cdp:[1085,1835], line:[1191,1400].
    int cdpNo_min=1085, lineNo_min=1191;
    int cdpNo_max=1835, lineNo_max=1400;
    float   mx1 = 0;
    float   my1 = 0;
    mx1 = 753581;
    my1 = 4949006;
    cdpNo_min=1251; lineNo_min=1281;
    cdpNo_max=1255; lineNo_max=1281;
    /*
    // for Luojia3D data: cdp:[2401,3700], line:[1991,2460].
    int cdpNo_min=3001, lineNo_min=2211;
    int cdpNo_max=3004, lineNo_max=2212;
    float   mx1 = 0;
    float   my1 = 0;
     */
    float   **foldMap   = alloc2float(ncdp_proj, nline);
    memset((void *)&foldMap[0][0], 0, sizeof(float)*ncdp_proj*nline);

    //#pragma omp parallel for shared(foldMap, MycmpFileIndex_2d)
    for ( int iline_no = lineNo_min ; iline_no <= lineNo_max; ++iline_no )
    {
        for ( int icdp_no  = cdpNo_min  ; icdp_no  <= cdpNo_max ; ++icdp_no  )
        {
            int itr_single  = 0;
            int lineNo  = iline_no;
            int cdpNo   = icdp_no;
            long   mx   = 0.;
            long   my   = 0.;
            CMPInfo MyCMPInfo;
            MyCMPInfo.itr_cdmpfile  = (long *)calloc(nfoldMax, sizeof(long));
            fbsegy_std  *suhdr  = calloc(nfoldMax, sizeof(fbsegy_std));
            memset((void *)MyCMPInfo.itr_cdmpfile, 0, sizeof(long)*nfoldMax);
            memset((void *)suhdr, 0, sizeof(fbsegy_std)*nfoldMax);

            // find the CDP gather with the given CDP&Line No. (lineNo, cdpNo)
            for ( int ifile = 0 ; ifile < nfile; ++ifile )
            {
                long    ntr = ntr_files[ifile];
                for ( long itr = 0 ; itr < ntr; ++itr )
                {
                    int iline   = MycmpFileIndex_2d[ifile][itr].line;
                    int icdp    = MycmpFileIndex_2d[ifile][itr].cdp;
                    if( lineNo == iline && cdpNo == icdp )
                    {
                        suhdr[itr_single].cdp    = cdpNo;
                        suhdr[itr_single].line   = lineNo;
                        suhdr[itr_single].cdpt   = itr_single+1;
                        suhdr[itr_single].tracl  = itr_single+1;
                        suhdr[itr_single].tracr  = itr_single+1;
                        suhdr[itr_single].sx = MycmpFileIndex_2d[ifile][itr].sx;
                        suhdr[itr_single].gx = MycmpFileIndex_2d[ifile][itr].gx;
                        suhdr[itr_single].sy = MycmpFileIndex_2d[ifile][itr].sy;
                        suhdr[itr_single].gy = MycmpFileIndex_2d[ifile][itr].gy;
                        suhdr[itr_single].ns = (unsigned short)ns;
                        suhdr[itr_single].dt = (unsigned short)dt_us;
                        suhdr[itr_single].f1 = 0;
                        suhdr[itr_single].d1 = 1.*dt_us/1000000;
                        suhdr[itr_single].f2 = 1;
                        suhdr[itr_single].d2 = 1.;

                        MyCMPInfo.ifile = ifile;
                        MyCMPInfo.itr_cdmpfile[itr_single] = itr;

                        float   sx  = (float)suhdr[itr_single].sx;
                        float   gx  = (float)suhdr[itr_single].gx;
                        float   sy  = (float)suhdr[itr_single].sy;
                        float   gy  = (float)suhdr[itr_single].gy;
                        float   offset  = sqrtf( (sx-gx)*(sx-gx) + (sy-gy)*(sy-gy) );
                        suhdr[itr_single].offset    = offset;
                        mx  += 0.5*(sx+gx);
                        my  += 0.5*(sy+gy);

                        ++ itr_single;
                    }
                }
            }
            int ntrCMP  = itr_single;
            MyCMPInfo.ntr_cdp   = ntrCMP;
            // calculate the fold number.
            foldMap[lineNo-lineMin][cdpNo-cdpMin_proj] = (float)ntrCMP;
            //printf("lineNo-lineMin=%d cdpNo-cdpMin_proj=%d \n", lineNo-lineMin, cdpNo-cdpMin_proj);

            // calculate the averaging mid-point coordinate.
            mx  /= ntrCMP;
            my  /= ntrCMP;
            /*
            if (icdp_no==cdpNo_min && iline_no==lineNo_min)
            {
                mx1 = mx;
                my1 = my;
            }
            //printf("ntrCMP=%d filename is %s\n", ntrCMP, CMPFileNameList[MyCMPInfo.ifile]);
            //printf("lineNo=%d cdpNo=%d ntrCMP=%d mx=%ld mxsft=%f my=%ld mysft=%f\n",
            //   lineNo, cdpNo, ntrCMP, mx, mx-mx1, my, my-my1);
             */
            printf("%ld\t%ld\n", mx, my);
            if ( 0 == ntrCMP )
            {
                fprintf(stderr, "Error! lineNo=%d cdpNo=%d ntr=%d\n", lineNo, cdpNo, ntrCMP);
                //mystat  = -1;
                //return mystat;
            }
            printf("ntr=%d\n", ntrCMP);
            fflush(stdout);
            fflush(stderr);
            //continue;

            char    cdpfilename_su[FILE_NAME_MAX_LENGTH]="";
            char    cdpfilename_dat[FILE_NAME_MAX_LENGTH]="";
            sprintf(cdpfilename_su, "%s/cdpgather_line%d_cdp%d.su", opwdata_dir, iline_no, icdp_no);
            sprintf(cdpfilename_dat, "%s/cdpgather_line%d_cdp%d.dat", opwdata_dir, iline_no, icdp_no);
            //fprintf(stdout, " cdpfilename is: %s\n", cdpfilename_su);
            // open the CMP file.
            open_cmp_gather(CMPFileNameList[MyCMPInfo.ifile]);

            float   **cmpgather = alloc2float(ns, ntrCMP);
            fbsegy_std  suhdr_tmp;
            for ( int itr = 0 ; itr < ntrCMP; ++itr )
            {
                long    jtr = MyCMPInfo.itr_cdmpfile[itr];
                long    fileOffset  = cmpfileOffsetBeg + 4L*(ns+60)*jtr;
                printf("itr=%d ntr=%d jtr=%ld\n", itr, ntrCMP, jtr);
                // Remember that "cmpGatherFilefp" is a global file pointer!
                fseek( cmpGatherFilefp, fileOffset, 0 );

                fread( &suhdr_tmp, sizeof(fbsegy_std), 1, cmpGatherFilefp );
                if ( ns != fread( &cmpgather[itr][0], 4, ns, cmpGatherFilefp ) )
                {
                    fprintf(stderr, "Read file: %s error! itr=%d\n", CMPFileNameList[MyCMPInfo.ifile], itr);
                    break;
                }

                if ( 0 == endian && 1 == conv )	// swapbytes of trace.
                {
                    for ( int it = 0 ; it < ns; ++it )
                        swap_float_4(&cmpgather[itr][it]);
                }
            }

            // output the current cdp gather.
            write_2d_float_wb(cmpgather, ntrCMP, ns, cdpfilename_dat);
            write_2d_float_wb_suhdr(cmpgather, suhdr, ntrCMP, ns, cdpfilename_su);

            // close the CMP file.
            close_cmp_gather();
            free2float(cmpgather);

            // free memory space.
            free(MyCMPInfo.itr_cdmpfile);
            free(suhdr);

            // display information.
            //fprintf(stdout, " \n------------------------------------\n");
            //fprintf(stdout, " lineNo = %d\t", iline_no);
            //fprintf(stdout, " cdpNo  = %d\t", icdp_no);
            //fprintf(stdout, " ntrCMP = %d\n", ntrCMP);

        }
    }
    //write_2d_float_wb(foldMap, nline, ncdp_proj, "./foldMap.dat");

    if( 1 == DisPlayCMPInfo )
    {
        fprintf(stdout, " ------------------------------------\n");
        fprintf(stdout, " Cheers!\n");
        fprintf(stdout, " CMP-index files scanned %d\n", nfile);
        fprintf(stdout, " Statistical information: \n");
        fprintf(stdout, " lineMin = %d\n", lineMin);
        fprintf(stdout, " lineMax = %d\n", lineMax);
        fprintf(stdout, " cdpMin  = %d\n", cdpMin_proj);
        fprintf(stdout, " cdpMax  = %d\n", cdpMax_proj);
        fprintf(stdout, " nfoldMax= %d\n", nfoldMax);
        fprintf(stdout, " ------------------------------------\n");
    }

    // free memory.
    free2char(CMPFileNameList);
    free2char(CMPIndexNameList);
    free2cmpFileIndex(MycmpFileIndex_2d);
    free(ntr_files);

    return mystat;
}
