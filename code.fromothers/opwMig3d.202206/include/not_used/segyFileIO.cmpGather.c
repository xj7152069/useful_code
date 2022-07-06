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

/*
int	init_cmp_gather(
// Input parameters.
	int	fileOffsetBeg,		// =0, SU-format; =3600, SEGY-format.
	long	ntr_cmpfile,		// total trace number in the CMP file.
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
	cmpGatherIndex *cmpfileidx
	)
{
	int	stat	= 0;

	// 1. read the segy-header of the first trace and get some parameters.
	fbsegy	segy_hdr	= {0} ;
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
	fread( &segy_hdr, sizeof(fbsegy) , 1 , cmpGatherFilefp );

	int	ns	= segy_hdr.ns;
	int	dtus	= segy_hdr.dt;
	float	dcdp	= segy_hdr.dcdp;
	float	dline	= segy_hdr.dline;
	*lt	= ns;
	*dt	= (float)dtus/1000.0;
	*dmx	= dcdp;
	*dmy	= dline;
printf("ns=%d dt=%d(us)\n", ns, dtus);

	int	cdp1	= segy_hdr.cdp;
	int	line1	= segy_hdr.line;
printf("cdp1=%d line1=%d\n", cdp1, line1);
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
printf("ntrmax=%d\n", ntrmax);

	*ncdp_first	= cdpmin;
	*ncdp_final 	= cdpmax;
	*nline_first	= linemin;
	*nline_final	= linemax;
	*ntr_gather_max	= ntrmax;

	return stat;
}
*/

int	init_cmp_gather_shengli(
// Input parameters.
	int	fileOffsetBeg,		// =0, SU-format; =3600, SEGY-format.
	long	ntr_cmpfile,		// total trace number in the CMP file.
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
	cmpGatherIndex *cmpfileidx
	)
{
	int	stat	= 0;

	// 1. read the segy-header of the first trace and get some parameters.
	fbsegy	segy_hdr	= {0} ;
	fseek( cmpGatherFilefp, fileOffsetBeg, 0 );
	fread( &segy_hdr, sizeof(fbsegy) , 1 , cmpGatherFilefp );
	// need swap-header.

	int	ns	= segy_hdr.ns;
	int	dtus	= segy_hdr.dt;
	float	dcdp	= 25.0;
	float	dline	= 25.0;
	//float	dcdp	= segy_hdr.dcdp;
	//float	dline	= segy_hdr.dline;
	*lt	= ns;
	*dt	= (float)dtus/1000.0;
	*dmx	= dcdp;
	*dmy	= dline;
printf("ns=%d dt=%d(us)\n", ns, dtus);

	int	kw181	= segy_hdr.lineno;
	int	kw185	= segy_hdr.cdpno;
	int	kw189	= segy_hdr.mx;
	int	kw193	= segy_hdr.my;
	int	kw201	= segy_hdr.ntr;
printf("kw181=%d kw185=%d kw189=%d kw193=%d kw201=%d\n", kw181, kw185, kw189, kw193, kw201);
	return	0;

	int	line1	= kw181;
	int	cdp1	= kw185;
/*
	int	cdp1	= segy_hdr.cdp;
	int	line1	= segy_hdr.lineno;
	int	ep	= segy_hdr.ep;
	int	fldr	= segy_hdr.fldr;
printf("cdp1=%d line1=%d ep=%d fldr=%d \n", cdp1, line1, ep, fldr);
 */

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

		int	cdp	= segy_hdr.cdpno;
		int	line	= segy_hdr.lineno;
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
printf("ntrmax=%d\n", ntrmax);

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
