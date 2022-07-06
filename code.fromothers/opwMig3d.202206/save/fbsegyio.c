/*====================================================================*/
/*===== Information Table of vpar & ivpar in gdbfcinitcutslicedata_().
*      ivpar[0] = 4th_num      ! Sample number of the 4th dim.
*      ivpar[1] = 3rd_num      ! Sample number of the 3rd dim.
*      ivpar[2] = third1       ! First sample number of the 3rd dim.
*      ivpar[3] = third2       ! Last  sample number of the 3rd dim.
*      ivpar[4] = dthird       ! Sampling   interval of the 3rd dim.
*      ivpar[5] = 2nd_num      ! Sample number of the 2nd dim.
*      ivpar[6] = second1      ! First sample number of the 2nd dim.
*      ivpar[7] = second2      ! Last  sample number of the 2nd dim.
*      ivpar[8] = dsecond      ! Sampling   interval of the 2nd dim.
*      ivpar[9] = 1th_num      ! Sample number of the 1th dim.
*      ivpar[10]= first1       ! First sample number of the 1th dim.
*      ivpar[11]= first2       ! Last  sample number of the 1th dim.
*      ivpar[12]= dfirst       ! Sampling   interval of the 1th dim.
*      ivpar[13]= nslice       ! The slice number of each file.*
*      ivpar[14]= format       ! Format(=1 float; =2 short; =3 long; =4 complex).
*      ivpar[15]= nf           ! Total file number.*
*      ivpar[16]= nfi          ! ith   file number.*
*      ivpar[17]= domain       ! File format.
*      ivpar[18]= ptype        ! Project format(=2 2D; =3 3D).

*      vpar[0]  = forth1       ! First value of the 4th dim.
*      vpar[1]  = ddforth      ! Real grid interval of the 4th dim.
*      vpar[2]  = ddthird      ! Real grid interval of the 3rd dim.
*      vpar[3]  = ddsecond     ! Real grid interval of the 2nd dim.
*      vpar[4]  = ddfirst      ! Real grid interval of the 1th dim.
*      vpar[5]  = vmin         ! Minimum value of data file.*
*      vpar[6]  = vmax         ! Maximum value of data file.*
*      vpar[7]  = x0           ! Origin coordinate of X.
*      vpar[8]  = y0           ! Origin coordinate of Y.
*/
int	fbioInitVel3dZXY()
{
	int	stat	= 0;

	return	stat;
}
//===== Initializing the Velocity Field File.      =======================//
/*
	gdbfcinitcutslicedata_( fn2, ivpar, vpar, &ifn2 );
	if( ivpar[5] == -99999 )
	{
		printf("  gdbfcinitcutslicedata_ %s error.\n", fn1);
		printf("  Program Exit; return 1; \n");
		return 1;
	}
	int	nvline_first, nvline_final, nvline_step, nvline_number, nvmy;
	int	nvcdp_first,  nvcdp_final,  nvcdp_step,  nvcdp_number,  nvmx;
	int	nvtau ;
	int	ivflag ;
	float	dvx , dvy, dvtau ;

	nvline_number = ivpar[1] ;
	nvline_first  = ivpar[2] ;
	nvline_final  = ivpar[3] ;
	nvline_step   = ivpar[4] ;
	nvcdp_number  = ivpar[5] ;
	nvcdp_first   = ivpar[6] ;
	nvcdp_final   = ivpar[7] ;
	nvcdp_step    = ivpar[8] ;
	nvtau         = ivpar[9] ;
	dvy           =  vpar[2] ;    // Warning: the unit of dt should be ms.
	dvx           =  vpar[3] ;
	dvtau         =  vpar[4] ;

	if( myid == 0 ) printf(" Initialize file %s done.\n ifn2=%d\n",fn2, ifn2 );

	nvcdp_number  = nvcdp_final  - nvcdp_first  + 1 ;
	nvline_number = nvline_final - nvline_first + 1 ;
	nvmx          = nvcdp_number  ;
	nvmy          = nvline_number ;

	if( myid == 0 )
	{
		printf("                                                      \n");
		printf(" ********  Parameters of 3d Velocity Field.  ******** \n");
		printf(" nvline_first   = %d\n", nvline_first  );
		printf(" nvline_final   = %d\n", nvline_final  );
		printf(" nvline_step    = %d\n", nvline_step   );
		printf(" nvline_number  = %d\n", nvmy          );
		printf(" nvcdp_first    = %d\n", nvcdp_first   );
		printf(" nvcdp_final    = %d\n", nvcdp_final   );
		printf(" nvcdp_step     = %d\n", nvcdp_step    );
		printf(" nvcdp_number   = %d\n", nvmx          );
		printf(" nvtau          = %d\n", nvtau         );
		printf(" dvx(m)         = %f\n", dvx           );
		printf(" dvy(m)         = %f\n", dvy           );
		printf(" dvtau(ms)      = %f\n", dvtau         );
		printf("                                                      \n");
	}
*/
