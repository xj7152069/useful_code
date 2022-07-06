
#include "fbCommon.h"
#include "fbswapbyte.h"
#include "fbsegy_shengli.h"
#include "cmpFileIndex.h"

// Global file-pointers.
FILE *cmpGatherFilefp = NULL;   // file pointer for operating the cmp-gather file.

int	open_cmp_gather(char *cmpGatherFilename)
{
	int	stat	= 0;

	// 1. open the cmp gather file.
	if ( NULL == (cmpGatherFilefp = fopen(cmpGatherFilename,"rb")) )
	{
		// open file failed.
		stat	= +1;
	}

	return stat;
}

int	close_cmp_gather()
{
	int	stat	= 0;

	fclose(cmpGatherFilefp);

	return stat;
}

int	init_cmp_gather_shengli(
// Input parameters.
	int	fileOffsetBeg,		// =0, SU-format; =3600, SEGY-format.
	int	conv,			// convert data to native format.
					// = 0 ; assume data is in native format.
	int	endian,			// set =0 for little-endian machines(PC's,DEC,etc.).
	long	ntr_cmpfile_max,
// Output parameters.
	int	*ncdp_first,
	int	*ncdp_final,
	int	*nline_first,
	int	*nline_final,
	int	*ntr_gather_max,
	long	*ntr_cmpfile_get,
	float	*Offset_X_Max,
	float	*Offset_Y_Max,
	float	*Offset_Max,
	float	*Offset_Min,
	cmpFileIndex *MycmpFileIndex,
	int	verbose			// =1, print debug information.
	)
{
	int	stat	= 0;
	unsigned short	ns, dtus;
	int	kw181, kw185, kw189, kw193, kw201;
	int	lineno, cdpno, mx, my, ntr;
	int	sx, sy, gx, gy;
	int	kw197, byte197;

	// 1. read the segy-header of the first trace and get some parameters.
	fbsegy_SL	segy_hdr	= {0} ;
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
	fread( &segy_hdr, sizeof(fbsegy_SL) , 1 , cmpGatherFilefp );
	if ( 0 == endian && 1 == conv )	// swapbytes of header keywords.
	{
		printf("conv=%d endian=%d \n", conv, endian);
		swap_u_short_2(&segy_hdr.ns);
		swap_u_short_2(&segy_hdr.dt);
		swap_int_4(&segy_hdr.lineno);
		swap_int_4(&segy_hdr.cdpno);
		swap_int_4(&segy_hdr.mx);
		swap_int_4(&segy_hdr.my);
		swap_int_4(&segy_hdr.ntr);
		swap_int_4(&segy_hdr.byte197);
	}

	ns	= segy_hdr.ns;
	dtus	= segy_hdr.dt;
	kw181	= segy_hdr.lineno;
	kw185	= segy_hdr.cdpno;
	kw189	= segy_hdr.mx;
	kw193	= segy_hdr.my;
	kw201	= segy_hdr.ntr;
	kw197	= segy_hdr.byte197;
	if ( 1 == verbose )
printf("ns=%d dtus=%d kw181=%d kw185=%d kw189=%d kw193=%d kw197=%d kw201=%d\n", ns, dtus, kw181, kw185, kw189, kw193, kw197, kw201);

	int	line1	= lineno;
	int	cdp1	= cdpno;
	int	MAX_INT	= (unsigned)(-1)>>1;
	int	cdpmin=MAX_INT, cdpmax=-1;
	int	linemin=MAX_INT, linemax=-1;
	int	ntrmin	= MAX_INT;
	int	ntrmax	= -1;
	float	offX_Max = 0;
	float	offY_Max = 0;
	float	off_Max = 0;
	float	off_Min	=1E10;
	if ( 1 == verbose )
printf("ns=%d dt=%d(us) line1=%d cdp1=%d MAX_INT=%d\n", ns, dtus, line1, cdp1, MAX_INT);

	// 2. read the segy-header of each trace.
	int	flag_eof = 0;
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
	long	itr	= 0;
	int	ntrace_cmp	= 0;
	while ( 1 )
	{
		fread( &segy_hdr, sizeof(fbsegy_SL) , 1 , cmpGatherFilefp ) ;
		flag_eof = feof(cmpGatherFilefp);
		if (flag_eof)
		{
			// read header failed.
			break;
		}
		//if( 0 == itr%10000 )	printf("itr=%ld\n", itr);

		if ( 0 == endian && 1 == conv )	// swapbytes of header keywords.
		{
			swap_u_short_2(&segy_hdr.ns);
			swap_u_short_2(&segy_hdr.dt);
			swap_int_4(&segy_hdr.lineno);
			swap_int_4(&segy_hdr.cdpno);
			swap_int_4(&segy_hdr.mx);
			swap_int_4(&segy_hdr.my);
			swap_int_4(&segy_hdr.ntr);
			swap_int_4(&segy_hdr.byte197);

			swap_int_4(&segy_hdr.sx);
			swap_int_4(&segy_hdr.sy);
			swap_int_4(&segy_hdr.gx);
			swap_int_4(&segy_hdr.gy);
		}

		ns	= segy_hdr.ns;
		dtus	= segy_hdr.dt;
		lineno	= segy_hdr.lineno;
		cdpno	= segy_hdr.cdpno;
		mx	= segy_hdr.mx;
		my	= segy_hdr.my;
		ntr	= segy_hdr.ntr;
		kw197	= segy_hdr.byte197;
		sx	= segy_hdr.sx;
		sy	= segy_hdr.sy;
		gx	= segy_hdr.gx;
		gy	= segy_hdr.gy;
/*
		float	mxf	= 0.5*(sx+gx);
		float	myf	= 0.5*(sy+gy);
if ( 1 == verbose )
printf("itr=%ld lineno=%d cdpno=%d mx=%d mxf=%f my=%d myf=%f kw197=%d ntr=%d\n", itr, lineno, cdpno, mx, my, mxf, myf, kw197, ntr);
*/

		int	line	= lineno;
		int	cdp	= cdpno;
		if (cdp < cdpmin)	cdpmin	= cdp;
		if (cdp > cdpmax)	cdpmax	= cdp;
		if (line < linemin)	linemin	= line;
		if (line > linemax)	linemax	= line;
		if ( cdp != cdp1 )	// a new CMP gather.
		{
			if (ntrace_cmp < ntrmin)	ntrmin	= ntrace_cmp;
			if (ntrace_cmp > ntrmax)	ntrmax	= ntrace_cmp;
			cdp1	= cdpno;
			ntrace_cmp	= 0;
			if ( line != line1 )	// a new CMP Line.
			{
				line1	= lineno;
				cdp1	= cdpno;
			}
		}

		// get the cdp&line number.
		MycmpFileIndex[itr].line= line;
		MycmpFileIndex[itr].cdp	= cdp;
		MycmpFileIndex[itr].sx	= sx;
		MycmpFileIndex[itr].sy	= sy;
		MycmpFileIndex[itr].gx	= gx;
		MycmpFileIndex[itr].gy	= gy;
		float	offset_x	= fabs(gx - sx);
		float	offset_y	= fabs(gy - sy);
		float	offset_xy	= sqrt(offset_x*offset_x + offset_y*offset_y);
		MycmpFileIndex[itr].offset	= (int)offset_xy;
		if (offset_x > offX_Max)	offX_Max	= offset_x;
		if (offset_y > offY_Max)	offY_Max	= offset_y;
		if (offset_xy > off_Max)	off_Max		= offset_xy;
		if (offset_xy < off_Min)	off_Min		= offset_xy;

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
			stat	= -1;
			break;
		}
	}
printf("ntrmax=%d\n", ntrmax);

	*ncdp_first	= cdpmin;
	*ncdp_final 	= cdpmax;
	*nline_first	= linemin;
	*nline_final	= linemax;
	*ntr_gather_max	= ntrmax;
	*ntr_cmpfile_get= itr;
	*Offset_X_Max	= offX_Max;
	*Offset_Y_Max	= offY_Max;
	*Offset_Max	= off_Max;
	*Offset_Min	= off_Min;

	return stat;
}

int	write_cmpFileIndex(
	char	*CMPIndexName,
	long	ntr_cmpfile,
	cmpFileIndex	*MycmpFileIndex )
{
	int	stat	= 0;

	// 1. open the cmp-index file.
	FILE	*fp	= NULL;
	if ( NULL == (fp = fopen(CMPIndexName,"wb")) )
	{
		// open file failed.
		stat	= +1;
	}
	else
		fprintf(stderr, "Now, writing file: %s\n", CMPIndexName);

	fwrite( &ntr_cmpfile, sizeof(long), 1, fp );

	fwrite( &MycmpFileIndex[0] , sizeof(cmpFileIndex) , ntr_cmpfile , fp ) ;

	fclose(fp);

	return stat;
}

int	read_cmpFileIndex(
	char	*CMPIndexName,
	long	*ntr_cmpfile,
	cmpFileIndex	*MycmpFileIndex )
{
	int	stat	= 0;

	// 1. open the cmp-index file.
	FILE	*fp	= NULL;
	if ( NULL == (fp = fopen(CMPIndexName,"rb")) )
	{
		// open file failed.
		stat	= +1;
	}
	else
		fprintf(stderr, "Now, reading file: %s\n", CMPIndexName);

	long ntr ;
	fread( &ntr, sizeof(long), 1, fp );

	fread( &MycmpFileIndex[0] , sizeof(cmpFileIndex) , ntr, fp ) ;

	fclose(fp);

	*ntr_cmpfile    = ntr ;

	return stat;
}
