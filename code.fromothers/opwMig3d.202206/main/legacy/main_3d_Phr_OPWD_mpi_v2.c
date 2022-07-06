//include
#include "fbopw3d.h"
// MPI
#include "mpi.h"

//  Global variable: sx = sx/Coor_S_R_Scale; gx = gx/Coor_S_R_Scale
float    Coor_S_R_Scale;

// function prototypes begin.
int outputOPWdata(
        int myid,
        int firstCMPLineNo_PWD, // the first CMP Line-No. for tau-phr transform.
        int lastCMPLineNo_PWD,  // the last CMP Line-No. for tau-phr transform.
        int nline_loc,          // the total line number for serially processing in each MPI process.
        // parameters of the project.
        ProjCMPInfo *MyProjCMPInfo,
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    *dir,
        // parameters for tau-phr transform.
        int nphr,
        int ntau,
        char    *OPWdataFilename,
        int verbose
        );
void write_1d_float_wb_long ( float *data , long n , char *fn );
void read_1d_float_rb_long ( float *data , long n , char *fn );
int getMinMaxLineNo(
        int nfile,
        cmpFileIndex    **MycmpFileIndex_2d,
        long *ntr_files,
        int *LineMin_get,
        int *LineMax_get,
        int *NfoldMin_get,
        int *NfoldMax_get,
        int verbose
        );
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
        );
int processingOneCMPLine(
        int lineNo,         // the given Line-No.
        ProjCMPInfo *MyProjCMPInfo,
        int nfile,          // number of CMP files.
        long *ntr_files,    // total trace number within each CMP file.
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    **CMPFileNameList,
        // specific parameters of the CMP files (might be adjusted).
        int fileOffsetBeg,  // =0, SU-format; =3600, SEGY-format.
        int conv,           // =0 ; assume data is in native format.
        int endian,         // set =0 for little-endian machines(PC's,DEC,etc.).
        int ns,
        int dt,
        // parameters for tau-phr transform.
        int ntrMin_PWD,
        int ntau,
        int itaus,
        int itaue,
        int nphr,
        float   phrmax,
        float   phrmin,
        float   dphr,
        float   ***taup3d,  //taup3d[ncdp_proj][nphr][ntau].
        int verbose);
int processingMuitipleCMPLines(
        int firstLineNo,    // the first Line-No.
        int lastLineNo,     // the last Line-No.
        int nline_loc,      // = lastLineNo - firstLineNo + 1
        ProjCMPInfo *MyProjCMPInfo,
        int nfile,          // number of CMP files.
        long *ntr_files,    // total trace number within each CMP file.
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    **CMPFileNameList,
        // specific parameters of the CMP files (might be adjusted).
        int fileOffsetBeg,  // =0, SU-format; =3600, SEGY-format.
        int conv,           // =0 ; assume data is in native format.
        int endian,         // set =0 for little-endian machines(PC's,DEC,etc.).
        int ns,
        int dt,
        // parameters for tau-phr transform.
        int ntrMin_PWD,
        int ntau,
        int itaus,
        int itaue,
        int nphr,
        float   phrmax,
        float   phrmin,
        float   dphr,
        char    *dir,
        int verbose);
int parallelProcessing(
        int nproc,
        int myid,
        int firstCMPLineNo_PWD, // the first CMP Line-No. for tau-phr transform.
        int lastCMPLineNo_PWD,  // the last CMP Line-No. for tau-phr transform.
        int nline_loc,          // the total line number for serially processing in each MPI process.
        // parameters of the project.
        ProjCMPInfo *MyProjCMPInfo,
        int nfile,          // number of CMP files.
        long *ntr_files,    // total trace number within each CMP file.
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    **CMPFileNameList,
        // specific parameters of the CMP files (might be adjusted).
        int fileOffsetBeg,  // =0, SU-format; =3600, SEGY-format.
        int conv,           // =0 ; assume data is in native format.
        int endian,         // set =0 for little-endian machines(PC's,DEC,etc.).
        int ns,
        int dt,
        // parameters for tau-phr transform.
        int ntrMin_PWD,
        int ntau,
        int itaus,
        int itaue,
        int nphr,
        float   phrmax,
        float   phrmin,
        float   dphr,
        char    *dir,
        int verbose
        );
        // function prototypes end.


        // parallel calculation of multiple CMP Lines with MPI+OpenMP.
        int parallelProcessing(
                int nproc,
                int myid,
                int firstCMPLineNo_PWD, // the first CMP Line-No. for tau-phr transform.
                int lastCMPLineNo_PWD,  // the last CMP Line-No. for tau-phr transform.
                int nline_loc,          // the total line number for serially processing in each MPI process.
                // parameters of the project.
                ProjCMPInfo *MyProjCMPInfo,
                int nfile,          // number of CMP files.
                long *ntr_files,    // total trace number within each CMP file.
                cmpFileIndex    **MycmpFileIndex_2d,
                CMPLineInfo *MyCMPLineInfo,
                char    **CMPFileNameList,
                // specific parameters of the CMP files (might be adjusted).
                int fileOffsetBeg,  // =0, SU-format; =3600, SEGY-format.
                int conv,           // =0 ; assume data is in native format.
                int endian,         // set =0 for little-endian machines(PC's,DEC,etc.).
                int ns,
                int dt,
                // parameters for tau-phr transform.
                int ntrMin_PWD,
                int ntau,
                int itaus,
                int itaue,
                int nphr,
                float   phrmax,
                float   phrmin,
                float   dphr,
                char    *dir,
                int verbose
                )
{
    MPI_Status Status ;
    double  stime, etime;
    stime   = MPI_Wtime();

    int ntask;      // number of task for parallel-computing.
    ntask = (int)( (lastCMPLineNo_PWD-firstCMPLineNo_PWD+1)/nline_loc ) + 1;

    int R_table[R_table_Max], S_table[S_table_Max] ;
    fflush(stdout);
    fflush(stderr);
    fprintf(stdout, "myid=%d nproc=%d ntask=%d\n", myid, nproc, ntask);

    // The parallel mode.
    if ( 0 == myid )    //  Master Node.
    {
        if ( 1 == verbose )
            printf("\t****Enter Func: %s().\tTask number: %6d Processor number: %4d \n\n", __FUNCTION__, ntask, nproc);

        for ( int   itask = 1 ; itask <= ntask ; ++ itask )
        {
            MPI_Recv(R_table, R_table_Max, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            int isend = Status.MPI_SOURCE ;

            S_table[0]  = 1 ;
            S_table[1]  = itask ;
            MPI_Send(S_table, S_table_Max, MPI_INT, isend, itask, MPI_COMM_WORLD);
        }

        for ( int inode = 1 ; inode < nproc ; ++ inode )
        {
            MPI_Recv(R_table, R_table_Max, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            int isend = Status.MPI_SOURCE ;

            S_table[0]  = 0 ;
            S_table[1]  = 0 ;
            MPI_Send(S_table, S_table_Max, MPI_INT, isend, inode, MPI_COMM_WORLD);
        }
    }
    else                //  Working Node.
    {
        S_table[0]  = 0 ;
        S_table[1]  = 0 ;
        MPI_Send(S_table, S_table_Max, MPI_INT, 0, 0, MPI_COMM_WORLD);

        while(1)
        {
            MPI_Recv(R_table, R_table_Max, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

            if ( 0 != R_table[0] )
            {
                S_table[0]  = R_table[0];
                S_table[1]  = R_table[1];
                int itask   = R_table[1];

                double  stime   = MPI_Wtime();
                int firstLineNo = firstCMPLineNo_PWD + (itask-1)*nline_loc;
                int lastLineNo  = firstLineNo + nline_loc - 1;
                if ( ntask == itask )
                    lastLineNo  = lastCMPLineNo_PWD;
                fprintf(stdout, "myid=%d Processing Line-No.[%d,%d]\n", myid, firstLineNo, lastLineNo);

                processingMuitipleCMPLines(
                        firstLineNo, lastLineNo, nline_loc,
                        MyProjCMPInfo, nfile, ntr_files,
                        MycmpFileIndex_2d, MyCMPLineInfo, CMPFileNameList,
                        // specific parameters of the CMP files (might be adjusted).
                        fileOffsetBeg, conv, endian, ns, dt,
                        // parameters for tau-phr transform.
                        ntrMin_PWD, ntau, itaus, itaue,
                        nphr, phrmax, phrmin, dphr, dir,
                        verbose);

                MPI_Send(S_table, S_table_Max, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }   // end of R_table
            else
                break;

        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    etime   = MPI_Wtime();
    if ( 0 == myid && 1 == verbose )
        printf("\t****Run Time of Func: %s() is %f\n",  __FUNCTION__, etime-stime);

    return 0;
}

// read a whole CMP-Line of a given Line-No from the multi-CMP files.
int processingOneCMPLine(
        int lineNo,         // the given Line-No.
        ProjCMPInfo *MyProjCMPInfo,
        int nfile,          // number of CMP files.
        long *ntr_files,    // total trace number within each CMP file.
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    **CMPFileNameList,
        // specific parameters of the CMP files (might be adjusted).
        int fileOffsetBeg,  // =0, SU-format; =3600, SEGY-format.
        int conv,           // =0 ; assume data is in native format.
        int endian,         // set =0 for little-endian machines(PC's,DEC,etc.).
        int ns,
        int dt,
        // parameters for tau-phr transform.
        int ntrMin_PWD,
        int ntau,
        int itaus,
        int itaue,
        int nphr,
        float   phrmax,
        float   phrmin,
        float   dphr,
        float   ***taup3d,  //taup3d[ncdp_proj][nphr][ntau].
        int verbose)
{
    int mystat = 0;
    int lineMin = MyProjCMPInfo->lineMin;
    int lineMax = MyProjCMPInfo->lineMax;
    //int nline   = MyProjCMPInfo->nline;
    int cdpMin_proj = MyProjCMPInfo->cdpMin_proj;
    int cdpMax_proj = MyProjCMPInfo->cdpMax_proj;
    int ncdp_proj   = cdpMax_proj - cdpMin_proj + 1;
    //int nfoldMax= MyProjCMPInfo->nfoldMax_proj;

    // (1) find the specific CMP-file where the given Line-No is stored.
    // P.S. Each CMP-Line must be stored in a single CMP-file.
    int ifile_get = -1;
    int ifile_old = -1;
    int ntrace_line_get = 0;
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr = ntr_files[ifile];
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int lineNo_now   = MycmpFileIndex_2d[ifile][itr].line;
            if( lineNo == lineNo_now )
            {
                ++ ntrace_line_get;
                ifile_get = ifile;

                if ( 1 == ntrace_line_get )
                {
                    ifile_old = ifile;
                }
                else
                {
                    if ( ifile_get != ifile_old )
                    {
                        fprintf(stderr, " \n\n Error! Each CMP-Line must be stored in a single CMP-file.\n");
                        fprintf(stderr, " Line-No = %d\n",   lineNo);
                        fprintf(stderr, " ifile_old = %d\n", ifile_old);
                        fprintf(stderr, " ifile_now = %d\n", ifile_get);
                        mystat = -1;
                        break;
                    }
                }
            }
        }
    }

    // (2) find the CDP range: [cdpMin_line, cdpMax_line]
    int iline_get   = -1;
    int ncdp_line   = -1;
    int ntrace_line = -1;
    int cdpMin_line = -1;
    int cdpMax_line = -1;
    for ( int ilineNo = lineMin; ilineNo <= lineMax; ++ ilineNo )
    {
        int iline   = ilineNo - lineMin;
        if( lineNo == MyCMPLineInfo[iline].lineNo )
        {
            ncdp_line   = MyCMPLineInfo[iline].ncdp_line;
            ntrace_line = MyCMPLineInfo[iline].ntrace_line;
            cdpMin_line = MyCMPLineInfo[iline].cdpMin_line;
            cdpMax_line = MyCMPLineInfo[iline].cdpMax_line;
            iline_get   = iline;
            fprintf(stdout, " lineNo=%d ncdp_line=%d ntrace=%d cdpMin=%d cdpMax=%d\n",
                    lineNo, ncdp_line, ntrace_line, cdpMin_line, cdpMax_line);
            break;
        }
    }

    if ( ntrace_line != ntrace_line_get )
    {
        fprintf(stderr, " \n\n Error! Trace number of the given CMP-Line does not match the index file.\n");
        fprintf(stderr, " Current Line-No is: %d\n", lineNo);
        fprintf(stderr, " Trace number scanned is %d\n", ntrace_line_get);
        fprintf(stderr, " Trace number in the index is %d\n", ntrace_line);
    }

    if ( 1 == verbose )
    {
        fprintf(stderr, " File index is %d\n", ifile_get);
        fprintf(stderr, " File name is %s\n", CMPFileNameList[ifile_get]);
        fprintf(stderr, " Trace number in the index is %d\n", ntrace_line);
        fprintf(stderr, " Trace number scanned is %d\n", ntrace_line_get);
    }

    // (3) search the fold number of each CMP gather and allocate space.
    CMPData *MyCMPData_line = calloc(ncdp_line, sizeof(CMPData));
    for ( int cdpNo = cdpMin_line ; cdpNo <= cdpMax_line; ++cdpNo )
    {
        int icdp    = cdpNo - cdpMin_line;
        int nfold_cdp   = 0;
        long    ntr = ntr_files[ifile_get];
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int lineNo_now   = MycmpFileIndex_2d[ifile_get][itr].line;
            int cdpNo_now    = MycmpFileIndex_2d[ifile_get][itr].cdp;
            if( lineNo == lineNo_now && cdpNo == cdpNo_now )
            {
                nfold_cdp   ++;
            }
        }

        MyCMPData_line[icdp].ntr_cdp    = nfold_cdp;
        if ( nfold_cdp > 0 )
        {
            MyCMPData_line[icdp].sx = calloc(nfold_cdp, sizeof(int));
            MyCMPData_line[icdp].sy = calloc(nfold_cdp, sizeof(int));
            MyCMPData_line[icdp].gx = calloc(nfold_cdp, sizeof(int));
            MyCMPData_line[icdp].gy = calloc(nfold_cdp, sizeof(int));
            MyCMPData_line[icdp].offset = calloc(nfold_cdp, sizeof(int));
            MyCMPData_line[icdp].itr_cdmpfile = calloc(nfold_cdp, sizeof(long));
            MyCMPData_line[icdp].gather = alloc2float(ns, nfold_cdp);
        }
    }

    // (4) read a whole CMP-line of the given Line-No.
    // open the CMP file.
    open_cmp_gather(CMPFileNameList[ifile_get]);

    /*
    // parameters for velocity semblance spectrum.
    float   dt_s= dt/1E6;
    float   fv  = 1500.0;
    float   dv  = 25.0;
    int     nv  = 121;
    float   ***semb3d   = alloc3float(ns, nv, ncdp_line);
    memset( (void*)&semb3d[0][0][0], 0, sizeof(float)*1L*ns*nv*ncdp_line );

    float   ***taupsection  = alloc3float(ntau, ncdp_proj, nphr);
    memset( (void*)&taupsection[0][0][0], 0, sizeof(float)*1L*ntau*ncdp_proj*nphr );
    float   **commonOffset   = alloc2float(ns, ncdp_proj);
    memset( (void*)&commonOffset[0][0], 0, sizeof(float)*1L*ns*ncdp_proj);
    */

    float   dt_ms  = dt/1E3;
    float   *vrms1d = calloc(ns, sizeof(float));
    float   vstk_min    = 1500.;
    float   vstk_max    = 4500.;
    for ( int it = 0 ; it < ns; ++it )
        vrms1d[it]  = vstk_min + (vstk_max-vstk_min)*it/ns;

    // loop for CDP point of the given Line-No.
#pragma omp parallel for shared(MyCMPData_line, MycmpFileIndex_2d)
    for ( int cdpNo = cdpMin_line ; cdpNo <= cdpMax_line; ++cdpNo )
    {
        fbsegy_std  Myfbsegy_std_tmp;
        int icdp    = cdpNo - cdpMin_line;
        int ntrCMP  = 0;
        long        ntr = ntr_files[ifile_get];
        int nfold_cdp   = MyCMPData_line[icdp].ntr_cdp;
        if ( nfold_cdp > 0 )
        {
            for ( long itr = 0 ; itr < ntr; ++itr )
            {
                int lineNo_now   = MycmpFileIndex_2d[ifile_get][itr].line;
                int cdpNo_now    = MycmpFileIndex_2d[ifile_get][itr].cdp;
                if( lineNo == lineNo_now && cdpNo == cdpNo_now )
                {
                    MyCMPData_line[icdp].sx[ntrCMP]   = MycmpFileIndex_2d[ifile_get][itr].sx;
                    MyCMPData_line[icdp].sy[ntrCMP]   = MycmpFileIndex_2d[ifile_get][itr].sy;
                    MyCMPData_line[icdp].gx[ntrCMP]   = MycmpFileIndex_2d[ifile_get][itr].gx;
                    MyCMPData_line[icdp].gy[ntrCMP]   = MycmpFileIndex_2d[ifile_get][itr].gy;
                    MyCMPData_line[icdp].itr_cdmpfile[ntrCMP]   = itr;
                    ++ ntrCMP;
                }
            }
        }
        else
            continue;
        fprintf(stderr, "lineNo=%d cdpNo=%d ntrCMP=%d nfold_cdp=%d\n", lineNo, cdpNo, ntrCMP, nfold_cdp);

        // read a single CMP-gather.
#pragma omp critical (read_single_CMP_gather)
        for ( int itr = 0 ; itr < ntrCMP; ++itr )
        {
            //printf("itr=%d ntr=%d\n", itr, ntrCMP);
            long    jtr = MyCMPData_line[icdp].itr_cdmpfile[itr];
            long    fileOffset  = fileOffsetBeg + 4L*(ns+60)*jtr;
            fseek( cmpGatherFilefp, fileOffset, 0 );

            fread( &Myfbsegy_std_tmp, sizeof(fbsegy_std), 1, cmpGatherFilefp );
            if ( ns != fread( &(MyCMPData_line[icdp].gather[itr][0]), 4, ns, cmpGatherFilefp ) )
            {
                fprintf(stderr, "Read file: %s error!\n", CMPFileNameList[ifile_get]);
                break;
            }

            if ( 0 == endian && 1 == conv ) // swapbytes of trace if necessary.
            {
                for ( int it = 0 ; it < ns; ++it )
                    swap_float_4(&(MyCMPData_line[icdp].gather[itr][it]));
            }
        }

        float   *offset = calloc(ntrCMP, sizeof(float));
        for ( int itr = 0 ; itr < ntrCMP; ++itr )
        {
            float   sx  = (float)MyCMPData_line[icdp].sx[itr]/Coor_S_R_Scale;
            float   sy  = (float)MyCMPData_line[icdp].sy[itr]/Coor_S_R_Scale;
            float   gx  = (float)MyCMPData_line[icdp].gx[itr]/Coor_S_R_Scale;
            float   gy  = (float)MyCMPData_line[icdp].gy[itr]/Coor_S_R_Scale;
            offset[itr] = sqrt( (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy) );
        }

        /*
        // output the common-offset section.
        int icdp_proj   = cdpNo - cdpMin_proj;
        float   offsetCenter    = 500.0;
        float   offsetErr       = 100.0;
        for ( int itr = 0 ; itr < ntrCMP; ++itr )
        {
        float   offsetNow  = offset[itr];
        if ( fabs(offsetNow - offsetCenter) < offsetErr )
        {
        for ( int it = 0 ; it < ns; ++it )
        {
        commonOffset[icdp_proj][it] = MyCMPData_line[icdp].gather[itr][it];
        }
        }
        }
        */

        // perform the tau-phr transform for a CMP gather.
        if (ntrCMP >= ntrMin_PWD)
        {
            //for field data, scaning for nan-value and zeroing it.
            for ( int itr = 0 ; itr < ntrCMP; ++itr )
            {
                for( int it  = 0 ; it  < ns ; ++ it  )
                {
                    if( isnan(MyCMPData_line[icdp].gather[itr][it]) )
                        MyCMPData_line[icdp].gather[itr][it] = 0.0 ;
                }
            }

            /*
            // calculate the velocity semblance spectrum.
            weightedSemblance(fv, dv, nv,
            MyCMPData_line[icdp].gather, ntrCMP, offset,
            ns, dt_s, semb3d[icdp]);
            */

            int icdp_proj   = cdpNo - cdpMin_proj;
            tauph_transform_3d_to_2d_robust(ns, ntau, nphr,
                    ntrCMP, itaus, itaue,
                    MyCMPData_line[icdp].gather, taup3d[icdp_proj],
                    vrms1d, offset, dt_ms, phrmin, dphr);

        }

        // free buffers.
        if ( nfold_cdp > 0 )
        {
            free(MyCMPData_line[icdp].sx);
            free(MyCMPData_line[icdp].sy);
            free(MyCMPData_line[icdp].gx);
            free(MyCMPData_line[icdp].gy);
            free(MyCMPData_line[icdp].offset);
            free(MyCMPData_line[icdp].itr_cdmpfile);
            free2float(MyCMPData_line[icdp].gather);
        }
        free(offset);

    }// end loop for CDP point of the given Line-No.

    // close the CMP file.
    close_cmp_gather();

    /*
       for ( int cdpNo = cdpMin_line ; cdpNo <= cdpMax_line; ++cdpNo )
       {
       int icdp_proj   = cdpNo - cdpMin_proj;
       for ( int iphr = 0 ; iphr < nphr; ++iphr )
       {
       for( int it  = 0 ; it  < ntau ; ++ it  )
       {
       taupsection[iphr][icdp_proj][it]    = taup3d[icdp_proj][iphr][it];
       }
       }
       }
       char    taupSectionFilename[FILE_NAME_MAX_LENGTH]="./taup3dsection_nphr_ncdp_ntau.dat";
       printf("ncdp_proj=%d nphr=%d ntau=%d\n", ncdp_proj, nphr, ntau);
       write_3d_float_wb(taupsection, nphr, ncdp_proj, ntau, taupSectionFilename);
       free3float(taupsection);

       char    COSectionFilename[FILE_NAME_MAX_LENGTH]="./commonOffsetsection.dat";
       write_2d_float_wb(commonOffset, ncdp_proj, ns, COSectionFilename);
       free2float(commonOffset);

       char    sembCMPLineFilename[FILE_NAME_MAX_LENGTH]="./sembLine.dat";
       write_3d_float_wb(semb3d, ncdp_line, nv, ns, sembCMPLineFilename);
       free3float(semb3d);

       char    taupCMPLineFilename[FILE_NAME_MAX_LENGTH]="./taup3dLine.dat";
       write_3d_float_wb(taup3d, ncdp_line, nphr, ntau, taupCMPLineFilename);
       free3float(taup3d);
       */

    /*
    // for debug: output the CMP line.
    long trace_idx_beg   = MyCMPData_line[0].itr_cdmpfile[0];
    int ntrtmp  = MyCMPData_line[ncdp_line-1].ntr_cdp;
    long trace_idx_end   = MyCMPData_line[ncdp_line-1].itr_cdmpfile[ntrtmp-1];
    printf("trace index begin: %ld\n", trace_idx_beg);
    printf("trace index end: %ld\n", trace_idx_end);

    char    cmpLineFilename[FILE_NAME_MAX_LENGTH]="./cmpLine.dat";
    float   **cmpgather_line   = alloc2float(ns, ntrace_line);
    int itrace_line = 0;
    for ( int cdpNo = cdpMin_line ; cdpNo <= cdpMax_line; ++cdpNo )
    {
    int icdp    = cdpNo - cdpMin_line;
    int nfold_cdp   = MyCMPData_line[icdp].ntr_cdp;
    for ( int itr = 0 ; itr < nfold_cdp; ++itr )
    {
    for ( int it = 0 ; it < ns; ++it )
    cmpgather_line[itrace_line][it] = MyCMPData_line[icdp].gather[itr][it];
    ++ itrace_line;
    }
    }
    printf("ntrace_line=%d itrace_line=%d\n", ntrace_line, itrace_line);
    write_2d_float_wb(cmpgather_line, ntrace_line, ns, cmpLineFilename);
    free2float(cmpgather_line);
    */

    free(MyCMPData_line);
    free(vrms1d);

    return mystat;
}

// processing multiple CMP-Lines serially.
int processingMuitipleCMPLines(
        int firstLineNo,    // the first Line-No.
        int lastLineNo,     // the last Line-No.
        int nline_loc,      // = lastLineNo - firstLineNo + 1
        ProjCMPInfo *MyProjCMPInfo,
        int nfile,          // number of CMP files.
        long *ntr_files,    // total trace number within each CMP file.
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    **CMPFileNameList,
        // specific parameters of the CMP files (might be adjusted).
        int fileOffsetBeg,  // =0, SU-format; =3600, SEGY-format.
        int conv,           // =0 ; assume data is in native format.
        int endian,         // set =0 for little-endian machines(PC's,DEC,etc.).
        int ns,
        int dt,
        // parameters for tau-phr transform.
        int ntrMin_PWD,
        int ntau,
        int itaus,
        int itaue,
        int nphr,
        float   phrmax,
        float   phrmin,
        float   dphr,
        char    *dir,
        int verbose)
{
    int mystat = 0;
    int cdpMin_proj = MyProjCMPInfo->cdpMin_proj;
    int cdpMax_proj = MyProjCMPInfo->cdpMax_proj;
    int ncdp_proj   = cdpMax_proj - cdpMin_proj + 1;

    char    taupSectionFilename[FILE_NAME_MAX_LENGTH]="";
    sprintf(taupSectionFilename,"%s/tauphrdata_nphr_nline_ncdp_ntau_lineNo%d_lineNo%d.dat", dir, firstLineNo, lastLineNo);
    fprintf(stdout, "output filename is %s\n", taupSectionFilename);

    float   ****taup4d   = alloc4float(ntau, ncdp_proj, nline_loc, nphr);
    memset( (void*)&taup4d[0][0][0][0], 0, sizeof(float)*1L*ntau*ncdp_proj*nline_loc*nphr);
    float   ***taup3d   = alloc3float(ntau, nphr, ncdp_proj);

    for ( int lineNo = firstLineNo; lineNo <= lastLineNo; ++ lineNo )
    {
        int iline   = lineNo - firstLineNo;
        memset( (void*)&taup3d[0][0][0], 0, sizeof(float)*1L*ntau*nphr*ncdp_proj );

        fflush(stdout);
        fprintf(stdout, " Processing Line-No: %d\n", lineNo);
        processingOneCMPLine(
                lineNo, MyProjCMPInfo, nfile, ntr_files,
                MycmpFileIndex_2d, MyCMPLineInfo, CMPFileNameList,
                // specific parameters of the CMP files (might be adjusted).
                fileOffsetBeg, conv, endian, ns, dt,
                // parameters for tau-phr transform.
                ntrMin_PWD, ntau, itaus, itaue,
                nphr, phrmax, phrmin, dphr, taup3d,
                verbose);

        for ( int cdpNo = cdpMin_proj ; cdpNo <= cdpMax_proj; ++cdpNo )
        {
            int icdp_proj   = cdpNo - cdpMin_proj;
            for ( int iphr = 0 ; iphr < nphr; ++iphr )
            {
                for( int it  = 0 ; it  < ntau ; ++ it  )
                {
                    taup4d[iphr][iline][icdp_proj][it]    = taup3d[icdp_proj][iphr][it];
                }
            }
        }
    }

    fprintf(stdout, "output filename is %s\n", taupSectionFilename);
    fprintf(stdout, "firstLineNo=%d lastLineNo=%d ncdp_proj=%d nphr=%d ntau=%d\n",
            firstLineNo, lastLineNo, ncdp_proj, nphr, ntau);
    write_1d_float_wb_long ( &taup4d[0][0][0][0], 1L*ntau*ncdp_proj*nline_loc*nphr, taupSectionFilename);

    // free memory.
    free3float(taup3d);
    free4float(taup4d);

    return mystat;
}

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

void write_1d_float_wb_long ( float *data , long n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file: %s to write\n" , fn ) ;
        exit(0) ;
    }
    fwrite( data , sizeof(float) , n , fp ) ;
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
    fread( data , sizeof(float) , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n",fn);
#endif
}

int main( int argc , char *argv[] )
{
    int myid, nproc;
    MPI_Init (&argc , &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int mystat  = 0;
    if( argc != 2 )
    {
        printf(" Wrong Parameters!\n");
        return 1 ;
    }
    char    parFile[FILE_NAME_MAX_LENGTH]="";
    strcpy(parFile, argv[1]);

    // init. of global variables.
    Coor_S_R_Scale  = 1.0;      // set the default value of global variable.

    char    dir[FILE_NAME_MAX_LENGTH]="./";
    // the maximum size of single CMP file.
    long    fileszieMaxByte = 1L*1024*1024*1024*1024;    // 1TB.
    int cmpfileOffsetBeg = 3600;    // = 0, SU-format; =3600, SEGY-format.
    int conv    = 1;                // = 0 ; assume data is in native format.
    int endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    int ns, dt;

    char    **CMPFileNameList   = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    **CMPIndexNameList  = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);

    int DisPlayCMPInfo  = 1;    // =1. display the CMP file information on the screen.
    int verbose = 1;
    int nfile   = -1;
    if ( 0 == myid )
        verbose = 1;
    else
        verbose = 0;

    readPar_OPWD(myid, parFile, &nfile, CMPFileNameList, CMPIndexNameList,
            dir, &fileszieMaxByte, &cmpfileOffsetBeg, &conv, &endian,
            &ns, &dt, &Coor_S_R_Scale, verbose);

    long    ntr_cmpfile_max = 100000000;
    int traceBytes   = 4*(ns+60);
    ntr_cmpfile_max = (long)( (fileszieMaxByte-cmpfileOffsetBeg) / traceBytes ) + 10;
    fprintf(stdout, " maximum trace number is %ld\n", ntr_cmpfile_max);

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
    int nline, lineMin, lineMax, nfoldMin, nfoldMax;
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

    if( 1 == DisPlayCMPInfo && 0 == myid )
    {
        fprintf(stdout, "\n\t ----- display project information. -----\n");
        fprintf(stdout, " ------------------------------------\n");
        fprintf(stdout, " Cheers!\n");
        fprintf(stdout, " CMP-index files scanned %d\n", nfile);
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
    //goto    jobDone;

    /*
    // for BN3D data.
    char    dir[FILE_NAME_MAX_LENGTH]="/ssd/ssd_raid0/fb/BN3D_test";
    int cmpfileOffsetBeg= 3600; // =0, SU-format; =3600, SEGY-format.
    int conv    = 1;            // = 0 ; assume data is in native format.
    int endian  = 0;            // set =0 for little-endian machines(PC's,DEC,etc.).
    int ns  = 6001;
    int dt  = 2000;             // unit of "dt" is us.
    Coor_S_R_Scale  = 1.0;      // set value of global variable.

    // for Luojia3D data: LJ2017_cmp_Line2201.su
    char    dir[FILE_NAME_MAX_LENGTH]="/ssd/ssd_raid0/fb/tmpdir/shengli_lj";
    int cmpfileOffsetBeg= 0; // =0, SU-format; =3600, SEGY-format.
    int conv    = 0;        // = 0 ; assume data is in native format.
    int endian  = 0;        // set =0 for little-endian machines(PC's,DEC,etc.).
    int ns  = 2901;
    int dt  = 2000;
    Coor_S_R_Scale  = 100.0;      // set value of global variable.
    */

    // parameters for tau-phr transform.
    int ntrMin_PWD   = 21;
    int ntau    = ns;
    int itaus   = 20;
    int itaue   = ns;
    int nphr    = 101;
    float   dphr, phrmin, phrmax;
    phrmax  = 800.;
    phrmax  = phrmax*1.0E-6 ;
    phrmin  = 0.0 ;
    dphr    = (phrmax-phrmin)/(nphr-1);

    // parameters for parallel-calculation.
    int nline_loc   = 10;
    int firstCMPLineNo_PWD  = lineMin;
    int lastCMPLineNo_PWD   = lineMax;

    /*
    // parallel calculation of multiple CMP Lines with MPI+OpenMP.
    verbose = 1;
    parallelProcessing(
            nproc, myid,
            firstCMPLineNo_PWD, lastCMPLineNo_PWD, nline_loc,
            &MyProjCMPInfo, nfile, ntr_files,
            MycmpFileIndex_2d, MyCMPLineInfo, CMPFileNameList,
            cmpfileOffsetBeg, conv, endian, ns, dt,
            ntrMin_PWD, ntau, itaus, itaue,
            nphr, phrmax, phrmin, dphr, dir,
            verbose);
    //goto    jobDone;
    */
    char    OPWdataFilename[FILE_NAME_MAX_LENGTH]="./";
    verbose = 1;
    outputOPWdata(
            myid, firstCMPLineNo_PWD, lastCMPLineNo_PWD, nline_loc,
            &MyProjCMPInfo, MycmpFileIndex_2d, MyCMPLineInfo,
            dir, nphr, ntau,
            OPWdataFilename,
            verbose);

    // merge and re-ordering temp-data.
    // the final output: opw4d[nhpr][nline][ncd][ntau].

    /*
    // parameters for serial-calculation.
    float   ***taup3d   = alloc3float(ntau, nphr, ncdp_proj);
    memset( (void*)&taup3d[0][0][0], 0, sizeof(float)*1L*ntau*nphr*ncdp_proj );
    int lineNo  = 2210;
    fflush(stdout);
    fprintf(stdout, " Processing Line-No: %d\n", lineNo);
    processingOneCMPLine(
    lineNo, &MyProjCMPInfo, nfile, ntr_files,
    MycmpFileIndex_2d, MyCMPLineInfo, CMPFileNameList,
    // specific parameters of the CMP files (might be adjusted).
    cmpfileOffsetBeg, conv, endian, ns, dt,
    // parameters for tau-phr transform.
    ntrMin_PWD, ntau, itaus, itaue,
    nphr, phrmax, phrmin, dphr, taup3d,
    verbose);
    char    taupCMPLineFilename[FILE_NAME_MAX_LENGTH]="./taup3dLine.dat";
    write_3d_float_wb(taup3d, ncdp_proj, nphr, ntau, taupCMPLineFilename);
    free3float(taup3d);
    */

    // free memory.
    free2char(CMPFileNameList);
    free2char(CMPIndexNameList);
    free2cmpFileIndex(MycmpFileIndex_2d);
    free(ntr_files);
    free(MyCMPLineInfo);

jobDone:
    MPI_Finalize();

    return mystat;
}

// The master node will output the 4-D OPW data.
int outputOPWdata(
        int myid,
        int firstCMPLineNo_PWD, // the first CMP Line-No. for tau-phr transform.
        int lastCMPLineNo_PWD,  // the last CMP Line-No. for tau-phr transform.
        int nline_loc,          // the total line number for serially processing in each MPI process.
        // parameters of the project.
        ProjCMPInfo *MyProjCMPInfo,
        cmpFileIndex    **MycmpFileIndex_2d,
        CMPLineInfo *MyCMPLineInfo,
        char    *dir,
        // parameters for tau-phr transform.
        int nphr,
        int ntau,
        char    *OPWdataFilename,
        int verbose
        )
{
    double  stime, etime;
    stime   = MPI_Wtime();

    if ( 0 != myid )
        goto    outputOPWdataWait;

    int cdpMin_proj = MyProjCMPInfo->cdpMin_proj;
    int cdpMax_proj = MyProjCMPInfo->cdpMax_proj;
    int ncdp_proj   = cdpMax_proj - cdpMin_proj + 1;

    int ntask;      // number of task for parallel-computing.
    ntask = (int)( (lastCMPLineNo_PWD-firstCMPLineNo_PWD+1)/nline_loc ) + 1;

    for ( int itask = 1 ; itask <= ntask ; ++ itask )
    {
        int firstLineNo = firstCMPLineNo_PWD + (itask-1)*nline_loc;
        int lastLineNo  = firstLineNo + nline_loc - 1;
        if ( ntask == itask )
            lastLineNo  = lastCMPLineNo_PWD;
        fprintf(stdout, "myid=%d Processing Line-No.[%d,%d]\n", myid, firstLineNo, lastLineNo);

        char    taupSectionFilename[FILE_NAME_MAX_LENGTH]="";
        sprintf(taupSectionFilename,"%s/tauphrdata_nphr_nline_ncdp_ntau_lineNo%d_lineNo%d.dat", dir, firstLineNo, lastLineNo);
        fprintf(stdout, "Now, reading file: %s\n", taupSectionFilename);

        float   bufSizeGB   = 4.0*ntau*ncdp_proj*nline_loc*nphr/1024./1024./1024.;
        fprintf(stderr, "Allocating space: sizeof(taup4d) is %f (GB)\n", bufSizeGB);
        /*
        float   ****taup4d   = alloc4float(ntau, ncdp_proj, nline_loc, nphr);
        if ( NULL == taup4d )
        {
            fprintf(stderr, "Error!\n");
            fprintf(stderr, "Failed to allocate 4-D space for buffer: taup4d[nphr][nline_loc][ncdp_proj][ntau]\n");
            break;
        }
        memset( (void*)&taup4d[0][0][0][0], 0, sizeof(float)*1L*ntau*ncdp_proj*nline_loc*nphr);
        read_1d_float_rb_long ( &taup4d[0][0][0][0], 1L*ntau*ncdp_proj*nline_loc*nphr, taupSectionFilename);

        for ( int lineNo = firstLineNo; lineNo <= lastLineNo; ++ lineNo )
        {
            int iline   = lineNo - firstLineNo;
            for ( int cdpNo = cdpMin_proj ; cdpNo <= cdpMax_proj; ++cdpNo )
            {
                int icdp_proj   = cdpNo - cdpMin_proj;
                for ( int iphr = 0 ; iphr < nphr; ++iphr )
                {
                    for( int it  = 0 ; it  < ntau ; ++ it  )
                    {
                        //taup4d[iphr][iline][icdp_proj][it]    = taup3d[icdp_proj][iphr][it];
                    }
                }
            }
        }
        // free memory.
        free4float(taup4d);
        */
    }

outputOPWdataWait:
    MPI_Barrier( MPI_COMM_WORLD );

    etime   = MPI_Wtime();
    if ( 0 == myid && 1 == verbose )
        printf("\t****Run Time of Func: %s() is %f\n",  __FUNCTION__, etime-stime);

    return 0;
}
