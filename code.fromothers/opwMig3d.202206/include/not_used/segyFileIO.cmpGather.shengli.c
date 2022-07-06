/*=================================================================*
*==================================================================*
*       Subroutines for cmp Gather file IO.                        *
*                                                                  *
*==================================================================*
*==================================================================*
*         Author:       Feng Bo                                    *
*         WPI, Tongji University, China.                           *
*                                                                  *
*         Version:      v1.                                        *
*                          2016-04                                 *
*                                                                  *
*    (1) use the global file pointer: FILE    *cmpGatherFilefp     *
*        for cmp-gather file opening and reading.                  *
*==================================================================*/

/* TYPEDEFS */
typedef	struct	{
	int	cdp;
	int	line;
} cmpGatherIndex;

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
	long	ntr_cmpfile,		// total trace number in the CMP file.
	int	conv,			// convert data to native format.
					// = 0 ; assume data is in native format.
	int	endian,			// set =0 for little-endian machines(PC's,DEC,etc.).
// Output parameters.
	int	*lt,
	float	*dt,			// Warning: the unit of dt should be ms.
	float	*dmx,			// cdp-x interval (m).
	float	*dmy,			// cdp-y interval (m).
	int	*ncdp_first,
	int	*ncdp_final,
	int	*nline_first,
	int	*nline_final,
	int	*ntr_gather_max,
	cmpGatherIndex *cmpfileidx,
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
	fbsegy	segy_hdr	= {0} ;
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
	fread( &segy_hdr, sizeof(fbsegy) , 1 , cmpGatherFilefp );
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
	//int	ep	= segy_hdr.ep;
	//int	fldr	= segy_hdr.fldr;

	int	MAX_INT	= (unsigned)(-1)>>1;
	int	cdpmin=MAX_INT, cdpmax=-1;
	int	linemin=MAX_INT, linemax=-1;
	int	ntrmin	= MAX_INT;
	int	ntrmax	= -1;
	if ( 1 == verbose )
printf("ns=%d dt=%d(us) line1=%d cdp1=%d MAX_INT=%d\n", ns, dtus, line1, cdp1, MAX_INT);

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
		float	mxf	= 0.5*(sx+gx);
		float	myf	= 0.5*(sy+gy);
if ( 1 == verbose )
printf("itr=%ld lineno=%d cdpno=%d mx=%d mxf=%f my=%d myf=%f kw197=%d ntr=%d\n", itr, lineno, cdpno, mx, my, mxf, myf, kw197, ntr);

		int	line	= lineno;
		int	cdp	= cdpno;
		if (cdp < cdpmin)	cdpmin	= cdp;
		if (cdp > cdpmax)	cdpmax	= cdp;
		if (line < linemin)	linemin	= line;
		if (line > linemax)	linemax	= line;
		if (ntr < ntrmin)	ntrmin	= ntr;
		if (ntr > ntrmax)	ntrmax	= ntr;
		/*
		if (cdp != cdp1 || line != line1)
			ntr	= 0;
		++ ntr;
		cdp1	= cdp;
		line1	= line;
		if (ntr > ntrmax)	ntrmax	= ntr;
		 */

		// get the cdp&line number.
		cmpfileidx[itr].cdp	= cdp;
		cmpfileidx[itr].line	= line;

		// Skip a trace.
		fseek( cmpGatherFilefp, (long)ns*4L, SEEK_CUR );
	}
printf("ntrmax=%d\n", ntrmax);

	float	dcdp	= 12.5;
	float	dline	= 12.5;
	*dmx	= dcdp;
	*dmy	= dline;
	*lt	= ns;
	*dt	= (float)dtus/1000.0;
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
printf("ntrget=%d\n", ntrget);

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
		printf("%d\t%d\t%d\t%d\n", (int)sx, (int)sy, (int)gx, (int)gy);

		if ( ns != fread( &tcmp[itr][0], 4, ns , cmpGatherFilefp ) )
		{
			// read trace failed.
			stat	= -2;
			break;
		}
		// need swap-byte if necessary.
	}

	// Free buffers.
	free(traceNumIndex);

	return stat;
}
