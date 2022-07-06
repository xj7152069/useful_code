/*=================================================================*
*==================================================================*
*       Subroutines for plane-wave Gather file IO.                 *
*                                                                  *
*==================================================================*
*==================================================================*
*         Author:       Feng Bo                                    *
*         WPI, Tongji University, China.                           *
*                                                                  *
*         Version:      v1.                                        *
*                          2016-04                                 *
*                                                                  *
*    (1) use the global file pointer: FILE    *pwdGatherFilefp     *
*        for pwd-gather file opening and reading.                  *
*==================================================================*/

int	open_pwd_gather(char *pwdGatherFilename, char *mode)
{
	int	stat	= 0;

	// 1. open the plane-wave gather file.
	if ( NULL == (pwdGatherFilefp = fopen(pwdGatherFilename, mode)) )
	{
		// open file failed.
		stat	= +1;
	}

	return stat;
}

int	close_pwd_gather()
{
	int	stat	= 0;

	fclose(pwdGatherFilefp);

	return stat;
}

int	write_pwd_gather(
// Input parameters.
	float	***pdat3d,		// pwd gather: pdat3d[nmx][nphr][ntau].
	int	ntau,
	int	nphr,
	int	nmx,
	int	nline_first,		// the first line number.
	int	iline			// the current line number.
	)
{
	int	stat	= 0;

	long	pwdLineSize	= 4L*ntau*nphr*nmx;
	long	fileOffset	= (iline - nline_first) * pwdLineSize ;
	fseek( pwdGatherFilefp, fileOffset, 0 );

	float	**pdat2d	= alloc2float(ntau, nmx);
	for( int iphr = 0 ; iphr < nphr ; ++ iphr )
	{
		for( int imx = 0 ; imx < nmx  ; ++ imx )
		for( int it  = 0 ; it  < ntau ; ++ it  )
			pdat2d[imx][it]	= pdat3d[imx][iphr][it];
	
		int	count	= nmx*ntau;
		if ( count != fwrite( &pdat2d[0][0], sizeof(float) , count, pwdGatherFilefp ) )
		{
			// write plane-wave gather failed.
			stat	= -1;
			break;
		}
	}

	// Free buffers.
	free2float(pdat2d);

	return stat;
}
