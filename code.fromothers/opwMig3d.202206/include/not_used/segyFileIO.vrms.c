/*=================================================================*
*==================================================================*
*       Subroutines for velocity file IO.                          *
*                                                                  *
*==================================================================*
*==================================================================*
*         Author:       Feng Bo                                    *
*         WPI, Tongji University, China.                           *
*                                                                  *
*         Version:      v1.                                        *
*                          2016-04                                 *
*                                                                  *
*    (1) use the global file pointer: FILE    *vrmsFilefp          *
*        for stacking velocity file opening and reading.           *
*==================================================================*/

/* TYPEDEFS */
typedef	struct	{
	int	cdp;
	int	line;
} vrmsGatherIndex;

int	open_stk_velocity(char *stkVelocityFilename)
{
	int	stat	= 0;

	// 1. open the velocity file.
	if ( NULL == (vrmsFilefp = fopen(stkVelocityFilename,"rb")) )
	{
		// open file failed.
		stat	= +1;
	}

	return stat;
}

int	close_stk_velocity()
{
	int	stat	= 0;

	fclose(vrmsFilefp);

	return stat;
}

int	init_stk_velocity(
// Input parameters.
	int	fileOffsetBeg,		// =0, SU-format; =3600, SEGY-format.
	long	ntr_vrmsfile,		// total trace number in the velocity file.
// Output parameters.
	int	*nvtau,
	float	*dvtau,			// Warning: the unit of dt should be ms.
	float	*dvx,			// cdp-x interval (m).
	float	*dvy,			// cdp-y interval (m).
	int	*nvcdp_first,
	int	*nvcdp_final,
	int	*nvline_first,
	int	*nvline_final,
	vrmsGatherIndex *vrmsfileidx
	)
{
	int	stat	= 0;

	// 1. read the segy-header of the first trace and get some parameters.
	fbsegy	segy_hdr	= {0} ;
	fseek( vrmsFilefp, fileOffsetBeg, 0 );
	fread( &segy_hdr, sizeof(fbsegy) , 1 , vrmsFilefp );

	int	ns	= segy_hdr.ns;
	int	dtus	= segy_hdr.dt;
	//float	dcdp	= segy_hdr.dcdp;
	//float	dline	= segy_hdr.dline;
	float	dcdp	= 25.0;
	float	dline	= 25.0;
	*nvtau	= ns;
	*dvtau	= (float)dtus/1000.0;
	*dvx	= dcdp;
	*dvy	= dline;

	int	cdp1	= segy_hdr.cdpno;
	int	line1	= segy_hdr.lineno;
	int	cdpmin=cdp1, cdpmax=cdp1;
	int	linemin=line1, linemax=line1;

	// 2. read the segy-header of each trace.
	fseek( vrmsFilefp, fileOffsetBeg, 0 );
	for( long itr = 0 ; itr < ntr_vrmsfile ; ++ itr )
	{
		if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , vrmsFilefp ) )
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

		// get the cdp&line number.
		vrmsfileidx[itr].cdp	= cdp;
		vrmsfileidx[itr].line	= line;

		// Skip a trace.
		fseek( vrmsFilefp, (long)ns*4L, SEEK_CUR );
	}

	*nvcdp_first	= cdpmin;
	*nvcdp_final 	= cdpmax;
	*nvline_first	= linemin;
	*nvline_final	= linemax;

	return stat;
}

int	read_stk_velocity(
// Input parameters.
	int	fileOffsetBeg,		// =0, SU-format; =3600, SEGY-format.
	vrmsGatherIndex *vrmsfileidx,
	long	ntr_vrmsfile,		// total trace number in the velocity file.
	int	nvtau,			// trace length of velocity file.
	int	iline,
	int	ncdp_first,		// the first cdp number.
	int	ncdp_final,		// the last  cdp number.
	int	ns,
// Output parameters.
	float	**vrms			// stacking velocity: vrms[nmx][ns]. (ns<=nvtau)
	)
{
	int	stat	= 0;

	int	nmx	= ncdp_final - ncdp_first + 1;
	long	*traceNumIndex	= calloc(nmx, sizeof(long));
	int	ntrget	= 0;

	// 1. search the trace number of the given line & cdp number.
	for( long itr = 0 ; itr < ntr_vrmsfile ; ++ itr )
	{
		if ( iline == vrmsfileidx[itr].line &&
			vrmsfileidx[itr].cdp >= ncdp_first &&
			vrmsfileidx[itr].cdp <= ncdp_final )
		{
			traceNumIndex[ntrget]	= itr;
			++ ntrget;
		}
	}
	if ( ntrget != nmx)
	{
		stat	= 1;	// read velocity error.
		return	stat;
	}

	// 2. read the stacking velocity by trace-number index.
	int	traceSize	= 4*nvtau + 240;	// size of SEGY-format trace.
	fbsegy	segy_hdr	= {0} ;
	for( int itr = 0 ; itr < ntrget ; ++ itr )
	{
		long	traceNum	= traceNumIndex[itr];
		long	fileOffset	= fileOffsetBeg + traceNum * traceSize;
		fseek( vrmsFilefp, fileOffset, 0 );

		if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , vrmsFilefp ) )
		{
			// read header failed.
			stat	= -1;
			break;
		}

		if ( ns != fread( &vrms[itr][0], 4, ns , vrmsFilefp ) )
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
