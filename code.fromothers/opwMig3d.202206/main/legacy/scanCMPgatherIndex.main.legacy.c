/*=====================================================================*
*                                                                      *
*                 Scan index for multiple CMP Files.                   *
*                                                                      *
*======================================================================*
*   Summary:                                                           *
*      The program will scan the CMP-index files.                      *
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

//System
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<limits.h>
#include<float.h>
#include<math.h>

// Global file-pointers.
FILE	*cmpGatherFilefp;	// file pointer for operating the cmp-gather file.

//include
#include "fbsegy_shengli.h"
#include "cmpFileIndex.h"
#include "fballoc.h"
#include "fballoc.c"
#include "fbswapbyte.c"
#include "cmpFileIndex.c"

//Definition of Macros.
#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef FILE_NAME_MAX_LENGTH
#define FILE_NAME_MAX_LENGTH 4096
#endif

#ifndef FILE_NUMBER_MAX
#define FILE_NUMBER_MAX	256
#endif

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

cmpFileIndex    **alloc2cmpFileIndex(size_t n1, size_t n2)
{
        return	(cmpFileIndex **)alloc2(n1, n2, sizeof(cmpFileIndex));
}

void	free2cmpFileIndex(cmpFileIndex **p)
{
        free2((void**)p);
}


char	**alloc2char(size_t n1, size_t n2);
void	free2char(char **p);
void	str_replace(char *str_src, int n, char * str_copy);
void	suffix_replace(char *str_source, char *str_find, char *str_to_copy);
int	readpar(
		char	*parFile,
		int	*fileNumber,
		char	**CMPFileNameList,
		char	**CMPIndexNameList,
		int	verbose
	       );

int	main( int argc , char *argv[] )
{
	int	mystat	= 0;
	int	verbose	= 1;
	int	lineMin, lineMax;
	char	parFile[FILE_NAME_MAX_LENGTH]="";
	strcpy(parFile, argv[1]);

	int	nfile	= -1;
	char	**CMPFileNameList	= alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
	char	**CMPIndexNameList	= alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
	readpar(parFile, &nfile, CMPFileNameList, CMPIndexNameList, verbose);

	// File-index for multiple CMP files.
	//cmpFilesInfo	*MycmpFilesInfo	= calloc(nfile, sizeof(cmpFilesInfo));

	// for BN3D data.
	int	cmpfileOffsetBeg= 3600;	// =0, SU-format; =3600, SEGY-format.
	int	conv	= 1;		// = 0 ; assume data is in native format.
	int	endian	= 0;		// set =0 for little-endian machines(PC's,DEC,etc.).
	long	ntr_cmpfile_max	= 100000000;
	int	DisPlayCMPInfo	= 1;	// =1. display the CMP file information on the screen.

	cmpFileIndex	*MycmpFileIndex	= calloc(ntr_cmpfile_max, sizeof(cmpFileIndex));

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
    lineMin  = 100000000;
    lineMax  = -1;
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr = ntr_files[ifile];
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int ilineNo   = MycmpFileIndex_2d[ifile][itr].line;
            if (ilineNo < lineMin)    lineMin  = ilineNo;
            if (ilineNo > lineMax)    lineMax  = ilineNo;
        }
    }
    fprintf(stdout, " lineMin = %d\n", lineMin);
    fprintf(stdout, " lineMax = %d\n", lineMax);

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

    /*
    // find a CMP-gather using a given cdp+line.
    int cdpNo=3200, lineNo=2421;
    cmpFileIndex    *MycmpFileIndex_single   = calloc(nfoldMax, sizeof(cmpFileIndex));
    int itr_single  = 0;
    for ( int ifile = 0 ; ifile < nfile; ++ifile )
    {
        long    ntr = ntr_files[ifile];
        for ( long itr = 0 ; itr < ntr; ++itr )
        {
            int iline   = MycmpFileIndex_2d[ifile][itr].line;
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
		//fprintf(stdout, " nfoldMax= %d\n", nfoldMax);
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

char	**alloc2char(size_t n1, size_t n2)
{
        return	(char**)alloc2(n1, n2, sizeof(char));
}

void	free2char(char **p)
{
        free2((void**)p);
}

void	suffix_replace(char *str_source, char *str_find, char *str_to_copy)
{
	//char str_source[50] = "the book the source the end!\n";
	//char str_find[10] = "the";
	//char str_to_copy[10] = "+++";
	char *p;

	//printf("%s", str_source);

	p = strstr(str_source, str_find);
	while(p){

		str_replace(p, strlen(str_find), str_to_copy);
		p = p + strlen(str_to_copy);
		p = strstr(p, str_find);
	}

	//printf("%s", str_source);

}

void str_replace(char *str_src, int n, char * str_copy)
{
	int lenofstr;
	char *tmp;
	lenofstr = strlen(str_copy);
	//string to copy is shorter than string to find
	if(lenofstr < n)
	{
		tmp = str_src+n;
		while(*tmp)
		{
			*(tmp-(n-lenofstr)) = *tmp; //n-lenofstr, length of moving
			tmp++;
		}
		*(tmp-(n-lenofstr)) = *tmp; //move '\0'
	}
	else if(lenofstr > n) //string to copy longer than string to find
	{
		tmp = str_src;
		while(*tmp) tmp++;
		while( tmp>=(str_src+n) )
		{
			*(tmp+(lenofstr-n)) = *tmp;
			tmp--;
		}
	}
	strncpy(str_src, str_copy, lenofstr);
}

int	readpar(
		char	*parFile,
		int	*fileNumber,
		char	**CMPFileNameList,
		char	**CMPIndexNameList,
		int	verbose
	       )
{
	char	tmp[FILE_NAME_MAX_LENGTH]="";
	char	Suffix_in[FILE_NAME_MAX_LENGTH]="";
	char	Suffix_out[FILE_NAME_MAX_LENGTH]="";

	FILE	*fp;
	if((fp=fopen(parFile,"r"))==NULL)
	{
		printf("Can not open the paramter file: %s to read!\n", parFile);
		return 1;
	}

	fscanf(fp,"%s%s", tmp, Suffix_in);
	fscanf(fp,"%s%s", tmp, Suffix_out);
	fscanf(fp,"%s%d", tmp, fileNumber);
	int	nfile	= *fileNumber;
	for (int i = 0 ; i < nfile; ++i)
	{
		fscanf(fp,"%s", CMPFileNameList[i]);

		strcpy(CMPIndexNameList[i], CMPFileNameList[i]);

		suffix_replace(CMPIndexNameList[i], Suffix_in, Suffix_out);
	}

	fclose(fp);

	if( verbose == 1 )
	{
		fprintf(stdout, "the input CMP-file suffix is: %s\n", Suffix_in);
		fprintf(stdout, "the output Index-file suffix is: %s\n", Suffix_out);
		fprintf(stdout, "number of CMP files is %d\n", nfile);

		fprintf(stdout, "Input Filename list:\n");
		for (int i = 0 ; i < nfile; ++i)
			fprintf(stdout, "the input CMP filename is: %s\n", CMPFileNameList[i]);
		fprintf(stdout, "\n");

		fprintf(stdout, "\nOutput Filename list:\n");
		for (int i = 0 ; i < nfile; ++i)
			fprintf(stdout, "the output Index filename is: %s\n", CMPIndexNameList[i]);
		fprintf(stdout, "\n");
	}

	return 0;
}

