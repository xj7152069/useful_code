//Sys.
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>

//Definition of Macros.
#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef FILE_NAME_MAX_LENGTH
#define FILE_NAME_MAX_LENGTH 256
#endif

#ifndef DEBUG_READ_WRITE_MODE
#define DEBUG_READ_WRITE_MODE
#endif

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y)) 
#endif

#include"/home/fb/pw3dlib/include/fbsegy.h"

int	open_cmp_gather(char *cmpGatherFilename)
{

	FILE	*cmpGatherFilefp=NULL;

	// 1. open the cmp gather file.
	int	stat	= 0;
	if ( NULL == (cmpGatherFilefp = fopen(cmpGatherFilename,"rb")) )
	{
		// open file failed.
		stat	= +1;
	}

	// 2. read the segy-header of the first trace and get some parameters.
	int	fileOffsetBeg	= 3600;
	fileOffsetBeg	= 0;//su format.
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );

	fbsegy	segy_hdr	= {0} ;
	fread( &segy_hdr, sizeof(fbsegy) , 1 , cmpGatherFilefp );

	int	ns	= segy_hdr.ns;
	int	dtus	= segy_hdr.dt;
	int	cdp1	= segy_hdr.cdp;
	int	line1	= segy_hdr.line;
printf("ns=%d dt=%d(us)\n", ns, dtus);
printf("cdp1=%d line1=%d\n", cdp1, line1);

	fclose(cmpGatherFilefp);

	return stat;
}

int main( int argc , char *argv[] )
{
	char	cmpfn[FILE_NAME_MAX_LENGTH]="";
	strcpy(cmpfn, argv[1]);

	open_cmp_gather(cmpfn);

}
/*
	//int	cdp1	= segy_hdr.cdpno;
	//int	line1	= segy_hdr.lineno;

	int	cdpmin=cdp1, cdpmax=cdp1;
	int	linemin=line1, linemax=line1;
	int	ntr	= 0;
	int	ntrmax	= 0;

	// 2. read the segy-header of each trace.
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
	for( long itr = 0 ; itr < ntr_cmpfile ; ++ itr )
	{
		if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , cmpGatherFilefp ) )
		{
			// read header failed.
			stat	= -1;
			break;
		}

		int	cdp	= segy_hdr.cdp;
		int	line	= segy_hdr.line;
		//int	cdp	= segy_hdr.cdpno;
		//int	line	= segy_hdr.lineno;
		if (cdp < cdpmin)	cdpmin	= cdp;
		if (cdp > cdpmax)	cdpmax	= cdp;
		if (line < linemin)	linemin	= line;
		if (line > linemax)	linemax	= line;
		if (cdp != cdp1 || line != line1)
			ntr	= 0;
		++ ntr;
		cdp1	= cdp;
		line1	= line;
		if (ntr > ntrmax)	ntrmax	= ntr;

		// get the cdp&line number.
		cmpfileidx[itr].cdp	= cdp;
		cmpfileidx[itr].line	= line;

		// Skip a trace.
		fseek( cmpGatherFilefp, (long)ns*4L, SEEK_CUR );
	}

	*ncdp_first	= cdpmin;
	*ncdp_final 	= cdpmax;
	*nline_first	= linemin;
	*nline_final	= linemax;
	*ntr_gather_max	= ntrmax;

	return stat;
}

int	read_cmp_gather(
// Input parameters.
	int	fileOffsetBeg,		// =0, SU-format; =3600, SEGY-format.
	cmpGatherIndex *cmpfileidx,
	long	ntr_cmpfile,		// total trace number in the CMP file.
	int	ns,
	int	iline,
	int	icdp,
	int	ntr_gather_max,		// maximum trace number within a CMP gather.
// Output parameters.
	float	**tcmp,			// CMP gather.
	float	*offset,		// abs-offset.
	int	*ntr_gather_get		// trace number within current CMP gather.
	)
{
	int	stat	= 0;

	long	*traceNumIndex	= calloc(ntr_gather_max, sizeof(long));
	int	ntrget	= 0;

	zero2float( tcmp,   ns, ntr_gather_max);
	zero1float( offset, ntr_gather_max );

	// 1. search the trace number of the given line & cdp number.
	for( long itr = 0 ; itr < ntr_cmpfile ; ++ itr )
	{
		if ( iline == cmpfileidx[itr].line && icdp == cmpfileidx[itr].cdp )
		{
			traceNumIndex[ntrget]	= itr;
			++ ntrget;
		}
	}
	*ntr_gather_get	= ntrget;

	// 2. read the cmp gather by trace-number index.
	int	traceSize	= 4*ns + 240;	// size of SEGY-format trace.
	fbsegy	segy_hdr	= {0} ;
	for( int itr = 0 ; itr < ntrget ; ++ itr )
	{
		long	traceNum	= traceNumIndex[itr];
		long	fileOffset	= fileOffsetBeg + traceNum * traceSize;
		fseek( cmpGatherFilefp, fileOffset, 0 );

		if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , cmpGatherFilefp ) )
		{
			// read header failed.
			stat	= -1;
			break;
		}

		int	sx = segy_hdr.sx;
		int	sy = segy_hdr.sy;
		int	gx = segy_hdr.gx;
		int	gy = segy_hdr.gy;
		offset[itr] = sqrt( (gx-sx)*(gx-sx)+(gy-sy)*(gy-sy) );
		//printf("%8d\t%8d\t%8d\t%8d\n", (int)sx, (int)sy, (int)gx, (int)gy);

		if ( ns != fread( &tcmp[itr][0], 4, ns , cmpGatherFilefp ) )
		{
			// read trace failed.
			stat	= -2;
			break;
		}
	}

	// Free buffers.
	free(traceNumIndex);

	return stat;
}
*/
