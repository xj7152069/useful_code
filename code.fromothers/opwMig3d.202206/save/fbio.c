
#define _Debug_ScanAcquisition_2d_ 1

typedef struct {
//	Contains the shot&receiver coordinates.
	float	sx ;
	float	sz ;
	int	ishot_num ;	// the current shot-point number.
	int	ntr ;		// the trace number in current shot-gather.
	float	*RecvInfo ;
	float	*SeisRecord ;
} ShotGatherInfo ;

int	ScanAcquisition_2d( char *fn, ShotGatherInfo *MyShotGatherInfo,
			int *num_shots, int *num_receivers, int storageflag )
{
/*	discription.    */
//	char *fn	The input acquisition filename.
//	ShotGatherInfo *MyShotGatherInfo	A struct-array which contains the source&receiver coordinates.
//	int *num_shots	Return the shot number in the acquisition file.
//	int *num_receivers	Return the receiver number in the acquisition file.
//	int *storageflag	==0, only return the value of num_shots and num_receivers.
//				!=0, store the source&receiver coordinates;

	FILE	*fp = NULL ;
	if( NULL == (fp=fopen(fn,"r")) )
	{
		fprintf(stderr, "Can not open file: %s to read.\n" , fn ) ;
		return 1 ;
	}

	float	coordinate_z, coordinate_x, tmp1, tmp2 ;
	int	flag ;

	int	num_shots_find	= 0 ;
	int	num_traces_max	= 0 ;
	int	num_traces1	= 0 ;
	int	num_traces_min	= 1000000 ;
	int	num_traces2	= 1000000 ;
	int	num_traces	= 0 ;
	long	num_lines	= 0 ;

	while ( fscanf(fp, "%f%f%f%f%d", &coordinate_z, &coordinate_x, &tmp1, &tmp2, &flag) != EOF )
	{
		if ( 0 == flag )
		{
			if ( 0 != storageflag )
			{
				MyShotGatherInfo[num_shots_find].sz = coordinate_z ;
				MyShotGatherInfo[num_shots_find].sx = coordinate_x ;
			}

			++ num_shots_find ;
			++ num_lines ;

			// Find the maximum and minimum receiver number within shot-gathers.
			if ( num_traces1 > num_traces_max )
				num_traces_max = num_traces1 ;

			if ( num_traces2 < num_traces_min )
				num_traces_min = num_traces2 ;

			num_traces1 = 0 ;
			num_traces2 = 0 ;
			num_traces  = 0 ;
		}
		else if ( 1 == flag )
		{
			if ( 0 != storageflag )
			{
				MyShotGatherInfo[num_shots_find-1].RecvInfo[num_traces*2]   = coordinate_z ;
				MyShotGatherInfo[num_shots_find-1].RecvInfo[num_traces*2+1] = coordinate_x ;
				MyShotGatherInfo[num_shots_find-1].ntr = num_traces + 1 ;
			}

			++ num_traces1 ;
			++ num_traces2 ;
			++ num_traces  ;
			++ num_lines ;
		}
		else
			return 1 ;
	}

	fclose(fp);
	// Find the maximum and minimum receiver number within shot-gathers.
	if ( num_traces1 > num_traces_max )
		num_traces_max = num_traces1 ;

	if ( num_traces2 < num_traces_min )
		num_traces_min = num_traces2 ;

#ifdef _Debug_ScanAcquisition_2d_
	if ( 0 != storageflag )
	{
		for ( int ishot = 0 ; ishot < num_shots_find ; ++ ishot )
			printf("\t\t Shot[%6d] has %6d receivers.\n", ishot+1, MyShotGatherInfo[ishot].ntr );
	}
#endif
	*num_shots	= num_shots_find ;
	*num_receivers	= num_traces_max ;

	return 0 ;
}

int	ScanAcquisitionSystem_2d( char *fn, ShotGatherInfo *MyShotGatherInfo,
			int *Num_shots, int *Num_receivers,  int *Num_traces,
			int storageflag, int verbose )
{
/*	discription.    */
//	char *fn	The input acquisition filename.
//	ShotGatherInfo *MyShotGatherInfo	A struct-array which contains the source&receiver coordinates.
//	int *num_shots	Return the shot number in the acquisition file.
//	int *num_receivers	Return the receiver number in the acquisition file.
//	int *storageflag	==0, only return the value of num_shots and num_receivers.
//				!=0, store the source&receiver coordinates;
//	Modified by FengBo, 2013.07.03, from ScanAcquisition_2d().

	FILE	*fp = NULL ;
	if( NULL == (fp=fopen(fn,"r")) )
	{
		fprintf(stderr, "Can not open file: %s to read.\n" , fn ) ;
		return 1 ;
	}

	float	coordinate_z, coordinate_x, tmp1, tmp2 ;
	int	flag ;

	int	num_shots_find	= 0 ;
	int	num_traces_max	= 0 ;
	int	num_traces1	= 0 ;
	int	num_traces_min	= 1000000 ;
	int	num_traces2	= 1000000 ;
	int	num_traces	= 0 ;
	int	num_lines	= 0 ;
	int	num_traces_total= 0 ;

	while ( fscanf(fp, "%f%f%f%f%d", &coordinate_z, &coordinate_x, &tmp1, &tmp2, &flag) != EOF )
	{
		if ( 0 == flag )
		{
			if ( 0 != storageflag )
			{
				MyShotGatherInfo[num_shots_find].sz = coordinate_z ;
				MyShotGatherInfo[num_shots_find].sx = coordinate_x ;
			}

			++ num_shots_find ;
			++ num_lines ;

			// Find the maximum and minimum receiver number within shot-gathers.
			if ( num_traces1 > num_traces_max )
				num_traces_max = num_traces1 ;

			if ( num_traces2 < num_traces_min )
				num_traces_min = num_traces2 ;

			num_traces1 = 0 ;
			num_traces2 = 0 ;
			num_traces  = 0 ;
		}
		else if ( 1 == flag )
		{
			if ( 0 != storageflag )
			{
				MyShotGatherInfo[num_shots_find-1].RecvInfo[num_traces*2]   = coordinate_z ;
				MyShotGatherInfo[num_shots_find-1].RecvInfo[num_traces*2+1] = coordinate_x ;
				MyShotGatherInfo[num_shots_find-1].ntr = num_traces + 1 ;
			}

			++ num_traces1 ;
			++ num_traces2 ;
			++ num_traces  ;
			++ num_lines ;
		}
		else
			return 1 ;
	}

	fclose(fp);

	// Find the maximum and minimum receiver number within shot-gathers.
	if ( num_traces1 > num_traces_max )
		num_traces_max = num_traces1 ;

	if ( num_traces2 < num_traces_min )
		num_traces_min = num_traces2 ;

	// the total trace number in the acquisition system.
	num_traces_total	= num_lines - num_shots_find ;

	if ( 1 == verbose )
	{
		fprintf(stderr, "\t Line number   scanned: %8d\n", num_lines);
		fprintf(stderr, "\t Shot number   scanned: %8d\n", num_shots_find);
		fprintf(stderr, "\t Max-receiver   number: %8d\n", num_traces_max);
		fprintf(stderr, "\t Min-receiver   number: %8d\n", num_traces_min);
		fprintf(stderr, "\t total-receiver number: %8d\n", num_traces_total);
	}

	*Num_shots	= num_shots_find ;
	*Num_receivers	= num_traces_max ;
	*Num_traces	= num_traces_total ;

	return 0 ;
}

int	ScanAcquisition_shotcoor_2d( char *fn, float *shotx, float *shotz,
		long *nline, int *nshot )
{
/*	discription.    */
//	char *fn	The input acquisition filename.
//	float *shotx	The shot X coordinates.
//	float *shotz	The shot Z coordinates.
//	long *nline	The line number of input file.
//	int *nshot	The shot number of input file.

	FILE	*fp = NULL ;
	if( NULL == (fp=fopen(fn,"r")) )
	{
		fprintf(stderr, "Can not open file: %s to read.\n" , fn ) ;
		return 1 ;
	}

	float	coordinate_z, coordinate_x, tmp1, tmp2 ;
	int	flag ;

	long	num_lines	= 0 ;
	int	num_shots_find	= 0 ;

	while ( fscanf(fp, "%f%f%f%f%d", &coordinate_z, &coordinate_x, &tmp1, &tmp2, &flag) != EOF )
	{
		if ( 0 == flag )
		{
			shotx[num_shots_find]	= coordinate_x ;
			shotz[num_shots_find]	= coordinate_z ;

			++ num_shots_find ;
			++ num_lines ;
		}
		else if ( 1 == flag )
		{
			++ num_lines ;
		}
		else
			return 1 ;
	}

	fclose(fp);

#ifdef _Debug_ScanAcquisition_shotcoor_2d_
	for ( int ishot = 0 ; ishot < num_shots_find ; ++ ishot )
		printf("\t\t sx = %f\tsz = %f\n", shotx[ishot], shotz[ishot]);
	printf("\t Line number scanned: %8ld\n", num_lines);
	printf("\t Shot number scanned: %8d\n", num_shots_find);
#endif
	*nline	= num_lines ;
	*nshot	= num_shots_find ;

	return 0 ;
}

typedef struct {
	int shot_number ;
	int ntrace ;
	long fileoffset ;
} shotfile_info ;


typedef struct {
	int	sx;
//	int	sz;
	int	gx;
//	int	gz;
	long	itr ;
} shotfile_index ;

int	scan_shotgather_2d( char *fn, int *shot_mark, int *ns_input, int *ns_output,
		int *nshot_num, int *ntrace_num, int*ntr_max, shotfile_info *shotinfo )
{
/*	Description:
 *	char *fn: The filename of the input shot-gather. [Input]
 *	int *shot_mark:	A flag which can indicate the shot-number. [Input]
 *	int *ns_input:	The user-input sampling-number. [Input]
 *	int *ns_output:	The sampling-number read from the SU-Header. [Output]
 *	int *nshot_num:	The real shot number of the the input shot-gather. [Output]
 *	int *ntrace_num:The real trace number of the the input shot-gather. [Output]
 */

	FILE	*fp = NULL ;
	if ( NULL == (fp=fopen(fn,"rb")) )
	{
		printf(" Open file: %s error.\n" , fn ) ;
		exit(1) ;
	}
	
	//	Check the input parameters.
	int	ns_in = (*ns_input) ;
	if ( 3 != (*shot_mark) && 5 != (*shot_mark) )
	{
		printf(" \t**** shot_mark error! (shot_mark should be 3 or 5).\n");
		printf(" shot_mark = %d \n", *shot_mark);
		return 1 ;
	}

	fbsegy	segy_hdr={0} ;
	if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , fp ) )
	{
		printf(" Read file: %s error. \n" , fn ) ;
		exit(1) ;
	}

	int	ns, ep_cur, ep_next ;

	if ( 3 == *shot_mark )
		ep_cur	= segy_hdr.fldr ;
	else if ( 5 == *shot_mark )
		ep_cur	= segy_hdr.ep ;

	ns	= segy_hdr.ns ;
	if ( ns_in != ns )
	{
		printf(" \t**** Input ns is %6d\n", ns_in);
		printf(" \t**** ns in file %s is %6d\n", fn, ns);
		printf(" \t**** Use the input ns. (default)");
		ns	= ns_in ;
	}
	*ns_output	= ns ;

	int	shot_num = 1 ;
	int	ntrace_max = -1 ;
	long	trace_num = 0L ;
	int	itracf = 0 ;

	rewind( fp );
	while (1)
	{
		if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , fp ) )
			break ;

		++ trace_num ;
		++ itracf ;

		if ( 3 == *shot_mark )
			ep_next	= segy_hdr.fldr ;
		else if ( 5 == *shot_mark )
			ep_next	= segy_hdr.ep ;

		if ( ep_cur != ep_next )
		{
			-- itracf ;
			//printf(" Current shot number %d ntr %d \n", ep_cur, itracf);
			if ( itracf > ntrace_max)
				ntrace_max = itracf ;

			shotinfo[shot_num].shot_number	= ep_cur ;
			shotinfo[shot_num].ntrace	= itracf ;
			shotinfo[shot_num+1].fileoffset	= ftell(fp) - 240L ;

			ep_cur	= ep_next ;

			++ shot_num ;
			itracf = 1 ;
		}

		fseek( fp, (long)ns*4L, SEEK_CUR );	// Skip a trace.
		if ( 0 != feof(fp) )
			break ;
	}

	// Set the info of the first and the last shot.
	shotinfo[1].fileoffset		= 0L;
	shotinfo[shot_num].shot_number	= ep_cur ;
	shotinfo[shot_num].ntrace	= itracf ;

	if ( itracf > ntrace_max)
		ntrace_max = itracf ;

	*nshot_num	= shot_num ;
	*ntrace_num	= trace_num ;
	*ntr_max	= ntrace_max ;

/*
	for ( int ishot = 1 ; ishot <= shot_num ; ++ ishot )
	{
		printf(" shotinfo[%d].shot_number=%d\n", ishot, shotinfo[ishot].shot_number);
		printf(" shotinfo[%d].ntrace=%d\n", ishot, shotinfo[ishot].ntrace);
		printf(" shotinfo[%d].fileoffset=%ld\n", ishot, shotinfo[ishot].fileoffset);
	}
*/

	fclose(fp);

	return 0 ;
}

int	index_shotgather_2d( char *fn, long ntrace_num, int ns,
		shotfile_index *shotindex, short *scalco )
{
/*	Description:
 *	char *fn: The filename of the input shot-gather. [Input]
 *	long ntrace_num:The real trace number of the the input shot-gather. [Input]
 */

	FILE	*fp = NULL ;
	if ( NULL == (fp=fopen(fn,"rb")) )
	{
		printf(" Open file: %s error.\n" , fn ) ;
		exit(1) ;
	}
	
	fbsegy	segy_hdr={0} ;
	long	itr=-1;
	for ( itr = 0 ; itr < ntrace_num ; ++ itr )
	{
		if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , fp ) )
			break ;
		else
		{
			shotindex[itr].sx	= segy_hdr.sx ;
			//shotindex[itr].sz	= segy_hdr.sz ;
			shotindex[itr].gx	= segy_hdr.gx ;
			//shotindex[itr].gz	= segy_hdr.gz ;
			shotindex[itr].itr	= itr ;
		}

		fseek( fp, (long)ns*4L, SEEK_CUR );	// Skip a trace.
		if ( 0 != feof(fp) )
			break ;
	}

/*
	if ( (itr+1) != ntrace_num )
	{
		fprintf(stderr, "Read file error.\n itr=%ld ntr=%ld\n", itr, ntrace_num);
		return 1 ;
	}
*/

	rewind( fp );
	if ( 1 != fread( &segy_hdr, sizeof(fbsegy) , 1 , fp ) )
		return 1 ;

	*scalco	= segy_hdr.scalco ;

	//printf("****scalco=%d\n", segy_hdr.scalco);

	fclose(fp);

	return 0 ;
}

int	get_su_header_from_gather( char *fn, int ntrace_num, int ns, fbsegy *segy_hdr )
{
/*	Description:
 *	char *fn: The filename of the input shot-gather. [Input]
 *	long ntrace_num:The real trace number of the the input shot-gather. [Input]
 *	fbsegy *segy_hdr: The SU-Header of each seismic trace. [Output]
 */

	FILE	*fp = NULL ;
	if ( NULL == (fp=fopen(fn,"rb")) )
	{
		printf(" Open file: %s error.\n" , fn ) ;
		exit(1) ;
	}
	
	long	itr=-1;
	for ( itr = 0 ; itr < ntrace_num ; ++ itr )
	{
		if ( 1 != fread( &segy_hdr[itr], sizeof(fbsegy) , 1 , fp ) )
			break ;

		fseek( fp, (long)ns*4L, SEEK_CUR );	// Skip a trace.

		if ( 0 != feof(fp) )
			break ;
	}

/*
	if ( (itr+1) != ntrace_num )
	{
		fprintf(stderr, "Read file error.\n itr=%ld ntr=%ld\n", itr, ntrace_num);
		return 1 ;
	}
*/

	fclose(fp);

	return 0 ;
}

int	read_shotgather_2d( char *fn, int *shot_mark,
		int shot_idx, int nshot_num, int ns,
		fbsegy *segy_hdr, float **shotgather,
		int *ntr_gather, shotfile_info *shotinfo )
{
	FILE	*fp = NULL ;
	if ( NULL == (fp=fopen(fn,"rb")) )
	{
		printf(" Open file: %s error.\n" , fn ) ;
		exit(1) ;
	}
	
	//	Check the input parameters.
	if ( 3 != (*shot_mark) && 5 != (*shot_mark) )
	{
		printf(" \t**** shot_mark error! (shot_mark should be 3 or 5).\n");
		printf(" shot_mark = %d \n", *shot_mark);
		return 1 ;
	}

	int	ishot_idx = -1 ;
	for ( int is = 1 ; is <= nshot_num ; ++ is )
	{
		if ( shot_idx == shotinfo[is].shot_number )
		{
			ishot_idx = is ;
			break ;
		}
	}
	//printf(" \t\t\t ishot = %d ishot_idx = %d\n", shot_idx, ishot_idx);
	if ( -1 == ishot_idx )
	{
		printf("\t **** IN FUNCTION: %s () \n", __FUNCTION__ );
		printf(" \t**** shot_idx error!\n");
		printf(" \t shot_idx = %d\n", shot_idx);
		exit(1) ;
	}
	//printf("ishot_idx = %d\n", ishot_idx);
	//printf("shotinfo[ishot_idx].fileoffset = %ld\n", shotinfo[ishot_idx].fileoffset);

	// Set the location of the file-pointer.
	fseek( fp, shotinfo[ishot_idx].fileoffset, 0 );

	*ntr_gather	= shotinfo[ishot_idx].ntrace ;

	for ( int itr = 0 ; itr < shotinfo[ishot_idx].ntrace ; ++ itr )
	{
		if ( 1 != fread( &segy_hdr[itr], sizeof(fbsegy) , 1 , fp ) )
		{
			printf(" Read SU-header error.\n") ;
			exit(1) ;
		}

		if ( ns != fread( &shotgather[itr][0], sizeof(float) , ns , fp ) )
		{
			printf(" Read SU-header error.\n") ;
			exit(1) ;
		}
	}

/*
	// Mute.
	for ( int itr = 0 ; itr < shotinfo[ishot_idx].ntrace ; ++ itr )
	{
		int	imutes = segy_hdr[itr].muts ;
		int	imutee = segy_hdr[itr].mute ;
		if ( imutes < 0 || imutes > ns )
			imutes = 0 ;
		if ( imutee < 0 || imutee > ns )
			imutes = ns ;
		for ( int it = imutes ; it < imutee ; ++ it )
			shotgather[itr][it] = 0 ;

		//fprintf(stderr, "itr=%d imutes=%d imutee=%d\n", itr, imutes, imutee);
	}
*/

	fclose(fp);

	return 0 ;
}
