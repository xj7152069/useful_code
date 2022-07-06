//include
#include "fbopw3d.h"

int main( int argc , char *argv[] )
{
    int mystat  = 0;
    char    parFile[FILE_NAME_MAX_LENGTH]="";
    strcpy(parFile, argv[1]);

    int DisPlayCMPInfo  = 1;    // =1. display the CMP file information on the screen.
    int verbose = 1;
    int nfile   = -1;
    char    **CMPFileNameList   = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    **CMPIndexNameList  = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    readpar(parFile, &nfile, CMPFileNameList, CMPIndexNameList, verbose);

    /*
    // for BN3D data.
    int cmpfileOffsetBeg= 3600; // =0, SU-format; =3600, SEGY-format.
    int conv    = 1;        // = 0 ; assume data is in native format.
    int endian  = 0;        // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max = 100000000;
    int ns  = 6001;
    int dt  = 2000;
    */

    // for Luojia3D data.
    int cmpfileOffsetBeg= 0; // =0, SU-format; =3600, SEGY-format.
    int conv    = 0;        // = 0 ; assume data is in native format.
    int endian  = 0;        // set =0 for little-endian machines(PC's,DEC,etc.).
    long    ntr_cmpfile_max = 100000000;
    int ns  = 2901;
    int dt  = 2000;

    cmpFileIndex    *MycmpFileIndex = calloc(ntr_cmpfile_max, sizeof(cmpFileIndex));

    // Read CMP-index files.
    fprintf(stdout, " \n------------------------------------\n");
    cmpFileIndex    **MycmpFileIndex_2d = alloc2cmpFileIndex(ntr_cmpfile_max, nfile);
    long *ntr_files  = calloc(nfile, sizeof(long));
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr_cmpfile;
        read_cmpFileIndex(CMPIndexNameList[ifile], &ntr_cmpfile, MycmpFileIndex_2d[ifile]);
        ntr_files[ifile]    = ntr_cmpfile;
        fprintf(stdout, " total trace number is %ld\n", ntr_cmpfile);
    }

    // get the minimum and maximum line number in the project.
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
    int nfoldMin    = ntrMin;
    int nfoldMax    = ntrMax;
    fprintf(stdout, " lineMin = %d\n", lineMin);
    fprintf(stdout, " lineMax = %d\n", lineMax);
    fprintf(stdout, " nfoldMin= %d\n", nfoldMin);
    fprintf(stdout, " nfoldMax= %d\n", nfoldMax);

    int nline   = lineMax - lineMin + 1;
    CMPLineInfo *MyCMPLineInfo  = calloc(nline, sizeof(CMPLineInfo));

    // get the minimum and maximum cdp number in each CMP line.
    int cdpMin_proj = MycmpFileIndex_2d[0][0].cdp;
    int cdpMax_proj = MycmpFileIndex_2d[0][0].cdp;
#pragma omp parallel for shared(cdpMin_proj, cdpMax_proj)
    for ( int iline = lineMin; iline <= lineMax; ++ iline )
    {
        int cdpMin_loc  = 100000000;
        int cdpMax_loc  = -1;
        int ntrace_line = 0;
        for ( int ifile = 0 ; ifile < nfile; ++ifile )
        {
            long    ntr = ntr_files[ifile];
            for ( long itr = 0 ; itr < ntr; ++itr )
            {
                int ilineNo   = MycmpFileIndex_2d[ifile][itr].line;
                int icdpNo    = MycmpFileIndex_2d[ifile][itr].cdp;
                if( iline == ilineNo )
                {
                    if (icdpNo < cdpMin_loc)    cdpMin_loc  = icdpNo;
                    if (icdpNo > cdpMax_loc)    cdpMax_loc  = icdpNo;
                    ++ ntrace_line;
                }
            }
        }

        MyCMPLineInfo[iline].lineNo         = iline;
        MyCMPLineInfo[iline].ncdp_line      = cdpMax_loc - cdpMin_loc + 1;
        MyCMPLineInfo[iline].ntrace_line    = ntrace_line;
        MyCMPLineInfo[iline].cdpMin_line    = cdpMin_loc;
        MyCMPLineInfo[iline].cdpMax_line    = cdpMax_loc;

        //fprintf(stdout, " lineNo=%d cdpMin=%d cdpMax=%d\n", iline, cdpMin_loc, cdpMax_loc);
        // update the global minimum and maximum cdp number in the project.
#pragma omp critical
        {
            if (cdpMin_loc < cdpMin_proj)    cdpMin_proj  = cdpMin_loc;
            if (cdpMax_loc > cdpMax_proj)    cdpMax_proj  = cdpMax_loc;
        }
    }

    // display project information.
    for ( int iline = lineMin; iline <= lineMax; ++ iline )
    {
        int cdpMin_loc  = MyCMPLineInfo[iline].cdpMin_line;
        int cdpMax_loc  = MyCMPLineInfo[iline].cdpMax_line;
        int lineNo      = MyCMPLineInfo[iline].lineNo;
        int ncdp_line   = MyCMPLineInfo[iline].ncdp_line;
        int ntrace_line = MyCMPLineInfo[iline].ntrace_line;
        fprintf(stdout, " lineNo=%d ncdp=%d ntrace=%d cdpMin=%d cdpMax=%d\n",
                lineNo, ncdp_line, ntrace_line, cdpMin_loc, cdpMax_loc);
    }

    ProjCMPInfo MyProjCMPInfo;
    MyProjCMPInfo.lineMin   = lineMin;
    MyProjCMPInfo.lineMax   = lineMax;
    MyProjCMPInfo.cdpMin_proj   = cdpMin_proj;
    MyProjCMPInfo.cdpMax_proj   = cdpMax_proj;
    fprintf(stdout, "\n\t ----- display project information. -----\n");
    fprintf(stdout, "nline=%d lineMin=%d lineMax=%d\n", lineMax-lineMin+1, lineMin, lineMax);
    fprintf(stdout, "ncdp=%d cdpMin=%d cdpMax=%d\n", cdpMax_proj-cdpMin_proj+1, cdpMin_proj, cdpMax_proj);

    // find a CMP-gather using a given cdp+line.
    int cdpNo=3509, lineNo=2430;
    cmpFileIndex    *MycmpFileIndex_single   = calloc(nfoldMax, sizeof(cmpFileIndex));
    fbsegy_std  *Myfbsegy_std   = calloc(nfoldMax, sizeof(fbsegy_std));
    int itr_single  = 0;
    int ntrCMP  = 0;
    CMPInfo MyCMPInfo;
    MyCMPInfo.itr_cdmpfile  = (long *)calloc(nfoldMax, sizeof(long));

    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr = ntr_files[ifile];
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int iline   = MycmpFileIndex_2d[ifile][itr].line;
            int icdp    = MycmpFileIndex_2d[ifile][itr].cdp;
            if( lineNo == iline && cdpNo == icdp )
            {
                MycmpFileIndex_single[itr_single].line  = iline;
                MycmpFileIndex_single[itr_single].cdp   = icdp;
                MycmpFileIndex_single[itr_single].sx   = MycmpFileIndex_2d[ifile][itr].sx;
                MycmpFileIndex_single[itr_single].sy   = MycmpFileIndex_2d[ifile][itr].sy;
                MycmpFileIndex_single[itr_single].gx   = MycmpFileIndex_2d[ifile][itr].gx;
                MycmpFileIndex_single[itr_single].gy   = MycmpFileIndex_2d[ifile][itr].gy;

                Myfbsegy_std[itr_single].cdp    = cdpNo;
                Myfbsegy_std[itr_single].line   = lineNo;
                Myfbsegy_std[itr_single].cdpt   = itr_single+1;
                Myfbsegy_std[itr_single].tracl  = itr_single+1;
                Myfbsegy_std[itr_single].tracr  = itr_single+1;
                Myfbsegy_std[itr_single].sx = MycmpFileIndex_2d[ifile][itr].sx;
                Myfbsegy_std[itr_single].gx = MycmpFileIndex_2d[ifile][itr].gx;
                Myfbsegy_std[itr_single].sy = MycmpFileIndex_2d[ifile][itr].sy;
                Myfbsegy_std[itr_single].gy = MycmpFileIndex_2d[ifile][itr].gy;
                Myfbsegy_std[itr_single].ns = (unsigned short)ns;
                Myfbsegy_std[itr_single].dt = (unsigned short)dt;
                Myfbsegy_std[itr_single].f1 = 0;
                Myfbsegy_std[itr_single].d1 = 1.*dt/1000000;
                Myfbsegy_std[itr_single].f2 = 1;
                Myfbsegy_std[itr_single].d2 = 1.;

                MyCMPInfo.ifile = ifile;
                MyCMPInfo.itr_cdmpfile[itr_single] = itr;

                ++ itr_single;
            }
        }
    }
    ntrCMP  = itr_single;
    MyCMPInfo.ntr_cdp   = ntrCMP;
    printf("ntrCMP=%d filename is %s\n", ntrCMP, CMPFileNameList[MyCMPInfo.ifile]);

    // open the CMP file.
    open_cmp_gather(CMPFileNameList[MyCMPInfo.ifile]);

    float   **cmpgather = alloc2float(ns, ntrCMP);
    fbsegy_std  Myfbsegy_std_tmp;
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
    {
        printf("itr=%d ntr=%d\n", itr, ntrCMP);
        long    jtr = MyCMPInfo.itr_cdmpfile[itr];
        long    fileOffset  = cmpfileOffsetBeg + 4L*(ns+60)*jtr;
	    fseek( cmpGatherFilefp, fileOffset, 0 );

	    fread( &Myfbsegy_std_tmp, sizeof(fbsegy_std), 1, cmpGatherFilefp );
		if ( ns != fread( &cmpgather[itr][0], 4, ns, cmpGatherFilefp ) )
        {
            fprintf(stderr, "Read file: %s error!\n", CMPFileNameList[MyCMPInfo.ifile]);
            break;
        }

	    if ( 0 == endian && 1 == conv )	// swapbytes of trace.
        {
            for ( int it = 0 ; it < ns; ++it )
		        swap_float_4(&cmpgather[itr][it]);
        }
    }
    //write_2d_float_wb(cmpgather, ntrCMP, ns, "./cdpfile.dat");

    char    cmpFile[FILE_NAME_MAX_LENGTH]="./cdpfile.su";
	FILE *fpcmp = NULL;
	if( (fpcmp=fopen(cmpFile,"wb")) == NULL )
	{
		printf("can not open file:%s to write.\n" , cmpFile );
		return  -1;
	}
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
	{
        printf("---itr=%d ntr=%d\n", itr, ntrCMP);
		fwrite( &Myfbsegy_std[itr], sizeof(fbsegy_std),  1, fpcmp );
		fwrite( &cmpgather[itr][0], sizeof(float),  ns, fpcmp );
	}
    fclose(fpcmp);

    // close the CMP file.
    close_cmp_gather();
    /*
     */

    // display information.
    fprintf(stdout, " \n------------------------------------\n");
    fprintf(stdout, " lineNo = %d\n", lineNo);
    fprintf(stdout, " cdpNo  = %d\n", cdpNo);
    fprintf(stdout, " ntrCMP = %d\n", ntrCMP);
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
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
    /*
    */

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
    free(MycmpFileIndex);
    free2cmpFileIndex(MycmpFileIndex_2d);
    free(ntr_files);

    return mystat;
}
