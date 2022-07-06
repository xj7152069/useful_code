/*=====================================================================*
*                                                                      *
*                 3D planewave decompostion                            *
*                                                                      *
*======================================================================*
*       Summary:                                                       *
*         The program performs 3d Plane-Wave Decompositon,             *
*             under iCluster sys.                                      *
*       Detail:                                                        *
*         Input data:                                                  *
*            1. the CMP gather.                                        *
*            2. the RMS velocity field.                                *
*         Output data:                                                 *
*            1. the planewave data with the following data format      *
*               pdat4d(ns,nphr,nmx,nmy).                               *
*======================================================================*
*       Author: Feng Bo                                                *
*         Date  : 2009.02-2010.07                                      *
*         Current Version: v16                                         *
*======================================================================*
*       History:                                                       *
*         Stable Version: v14.                                         *
*         Date  : 2009.08                                              *
* Name:    3d_Phx_planewave_decomposition_for_iCluster_v13.c           *
* Package:  pw3dlib.update.090820.tar                                  *
*======================================================================*
*         Logs.                                                        *
*======================================================================*
*         Updated for V15, 2010-07-15                                  *
*           Modified regular-file writting, from data() to file().     *
*           Modified: 2010-07, with OpenMP-support.                    *
*           Using Multi-threads CDP grouping computing.                *
*         Updated for V16, 2010-07-17                                  *
*           Adding Nan-value scaning and zeroing features.             *
*         Updated for V16, 2016-04-23                                  *
*           write new IO functions for SEGY and high-dim data,         *
*           (1) segyFileIO.cmpGather.c                                 *
*           (2) segyFileIO.vrms.c                                      *
*           (3) regular4dFileIO.pwdGather.c                            *
*=====================================================================*/


//Sys.
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<math.h>

//OpenMP
#include<omp.h>

//MPI
#include "mpi.h"

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

// Define the global file pointers.
FILE	*cmpGatherFilefp;	// file pointer for reading the cmp-gather file.
FILE	*pwdGatherFilefp;	// file pointer for writing the pwd-gather file.
FILE	*vrmsFilefp;		// file pointer for reading the stacking velocity file.

//iCluster Headers.
//#include "gDBFCFuncIo.h"

//include
#include "fbsegy.h"
#include "fballoc.h"
#include "fbrw.h"
#include "fballoc.c"
#include "fbrw.c"
#include "sinc_interpolation.c"
#include "segyFileIO.cmpGather.c"
#include "segyFileIO.vrms.c"
#include "regular4dFileIO.pwdGather.c"

//src
#include "tauph_transform.c"

/*	Function Prototypes.	*/
int	readpar(char *parfn, char *fn1, char *fn2, char *fn3, int *ntau, int *nphr, float *phrmax, int *diy,
		int *nline_first_do, int *nline_final_do, int *ncdp_first_do, int *ncdp_final_do);
void	echopar(char *fn1, char *fn2, char *fn3,
		int ntau, int nline_first, int nline_final, int ncdp_first, int ncdp_final);
/*
int	read_vrms_of_current_line_tmp( float **velin, float dt_in, int nt_in,
		float **vrms, float dt_out, int ns, 
		int ncdp_first, int ncdp_final);
int	read_vrms_of_current_line( float **vrms, int iline, int ncdp_first, int ncdp_final, int ns, int ifn2);
int	read_cmp_gather_of_current_line(float **tcmp, float *offset,
		int iline, int icdp, int ntr, int *ntr_gather, int ns, int ifn1);
int	read_cmp_gather_of_current_line_new(float **tcmp, float *offset,
		int iline, int icdp, int ntr, int *ntr_gather, int ns, int ifn1);
int	write_planewave_data4d( float **pdat, int ns, int ntau, int nphr, int icdp, int iline, int ifn3);
int	write_planewave_data4d_group( float ***pdat3d, int ntau, int nphr, int cdp1, int cdp2, int iline, int ifn3);
float	sinc_interpolation( float *trace , int lt , float dt , float time_s_r );
int	init_planewave_file(char *fn3, int *ifn3, int *icpar, float *cpar,
	int ntau, float dt, 
	int nmx,  float dmx,  int ncdp_first,  int ncdp_final,  int ncdp_step,  float x0min,
	int nmy,  float dmy,  int nline_first, int nline_final, int nline_step, float y0min,
	int nphr, float dphr, int iphr1, int iphr2, int iphrstep, float phrmin);
*/

int main( int argc , char *argv[] )
{
	char	fn1[FILE_NAME_MAX_LENGTH] ;       // the input 3d CMP gather.
	char	fn2[FILE_NAME_MAX_LENGTH] ;       // the input RMS velocity field.
	char	fn3[FILE_NAME_MAX_LENGTH] ;       // the output planewave data.
	char	parfn[FILE_NAME_MAX_LENGTH] ;
	int 	ifn1=1, ifn2=2, ifn3=3;

	int	ispar[20]={0}, ivpar[25]={0}, icpar[25]={0};
	float	spar[20]={0}, vpar[20]={0}, cpar[20]={0};

	int	nline_first, nline_final, nline_step, nline_number, nmy ;
	int	ncdp_first,  ncdp_final,  ncdp_step,  ncdp_number,  nmx ;
	int	ns, ntr_all, ntr, ntau ;
	float	dt, dmx, dhx, dmy, amp, offset_max, offset_min ;
    
	int	nphr ;
	float	dphr, phrmin, phrmax;

	int	diy, nline_first_do, nline_final_do, ncdp_first_do,  ncdp_final_do ;
	int	myid, np;
	int	mystat	= 0;

	MPI_Status Status ;
	MPI_Init (&argc , &argv ) ;
	MPI_Comm_rank(MPI_COMM_WORLD , &myid ) ;
	MPI_Comm_size(MPI_COMM_WORLD , &np ) ;

	if( argc != 2 )
	{
		printf(" Wrong Parameters!\n");
		return 1 ;
	}

	// parameters init.
	int	cmpfileOffsetBeg	= 0;	// =0, SU-format; =3600, SEGY-format.
	int	vrmsfileOffsetBeg	= 0;	// =0, SU-format; =3600, SEGY-format.
	int	ntr_cmpfile	= 2304024;	// trace number of salt3d cmp gather.
	int	ntr_vrmsfile	= 10000;

	// for test only.
	cmpfileOffsetBeg= 0;	// default: cmp-gather is SU-format.
	//ntr_cmpfile	= 2024;	// trace number of salt3d cmp gather.

	nline_first_do = -1 ;
	nline_final_do = -1 ;
	ncdp_first_do  = -1 ;
	ncdp_final_do  = -1 ;

	strcpy(parfn,argv[1]);
	readpar(parfn, fn1, fn2, fn3, &ntau, &nphr, &phrmax, &diy,
		&nline_first_do, &nline_final_do, &ncdp_first_do, &ncdp_final_do);

/*====================================================================*/
/*===== Information Table of spar & ispar in gdbfcinitreadsegy_().
*      ispar[0] = ns           ! Sample length
*      ispar[1] = ntr_all      ! Trace number in a file.
*      ispar[2] = ntr_gather   ! Trace number in a gather.
*      ispar[3] = line1        ! First line number.
*      ispar[4] = dmy          ! Line step.
*      ispar[5] = nline        ! Total line number.
*      ispar[6] = shot1/cdp1   ! First shot/cdp number.
*      ispar[7] = dshot/dcdp   ! Shot/cdp step.
*      ispar[8] = nshot/ncdp   ! Total shot/cdp number.
*      ispar[9] = format       ! Format.
*      spar[0]  = dt           ! Sample rate(ms).
*      spar[1]  = ddcdp        ! Cdp step(m).
*      spar[2]  = ddline       ! Line step(m).
*      spar[3]  = 1.0/Amp_max  ! 1.0/Amp_max
*      spar[4]  = offset1      ! Minimum offset.
*      spar[5]  = offset2      ! Maximum offset.
*/
/*===== Initializing the CMP gather file.      =======================*/
/*
	gdbfcinitreadsegy_( fn1, ispar, spar, &ifn1);
	if( ispar[0] == 0 )
	{
		printf("  gdbfcinitreadsegy %s error.\n", fn1);
		printf("  Program Exit; return 1; \n");
		return 1;
	}

	ns            = ispar[0];
	ntr           = ispar[2];
	nline_first   = ispar[3];
	nline_step    = ispar[4];
	nline_number  = ispar[5];
	ncdp_first    = ispar[6];
	ncdp_step     = ispar[7];
	ncdp_number   = ispar[8];
	dt            =  spar[0];     // Warning: the unit of dt should be ms. //
	dmx           =  spar[1];
	dmy           =  spar[2];
	offset_min    =  spar[4];
	offset_max    =  spar[5];

	if( myid == 0 ) printf(" Initialize file %s done.\n ifn1=%d\n",fn1, ifn1 );
*/

	// open and init. cmp gather file.
	int	ostat	= open_cmp_gather( fn1 );
	if ( 0 != ostat )
	{
		mystat	= 11;	// open cmp gather file failed.
		goto	mainProgStop;
	}

	cmpGatherIndex	*cmpfileidx	= calloc(ntr_cmpfile, sizeof(cmpGatherIndex));
	int	ntr_gather_max	= 1;
	int	istat	= init_cmp_gather( cmpfileOffsetBeg, ntr_cmpfile,
		&ns, &dt, &dmx, &dmy,
		&ncdp_first, &ncdp_final,
		&nline_first, &nline_final,
		&ntr_gather_max, cmpfileidx);
	if ( 0 != istat )
	{
		mystat	= 12;	// read cmp gather file failed.
		goto	mainProgStop;
	}
	if( myid == 0 ) printf(" Initialize file %s done.\n", fn1);
	nline_step	= 1;
	ncdp_step	= 1;
	ntr		= ntr_gather_max;

	if( myid == 0 )
	{
		printf("                                                  \n");
		printf(" ********  Parameters of 3d CMP Gather.  ******** \n");
		printf(" nline_first   = %d\n", nline_first  );
		printf(" nline_final   = %d\n", nline_final  );
		printf(" nline_step    = %d\n", nline_step   );
		printf(" nline_number  = %d\n", nmy          );
		printf(" ncdp_first    = %d\n", ncdp_first   );
		printf(" ncdp_final    = %d\n", ncdp_final   );
		printf(" ncdp_step     = %d\n", ncdp_step    );
		printf(" ncdp_number   = %d\n", nmx          );
		printf(" ns            = %d\n", ns           );
		printf(" ntr           = %d\n", ntr          );
		printf(" dt(ms)        = %f\n", dt           );
		printf(" dmx(m)        = %f\n", dmx          );
		printf(" dmy(m)        = %f\n", dmy          );
		printf(" offset_min(m) = %f\n", offset_min   );
		printf(" offset_max(m) = %f\n", offset_max   );
		printf("                                                  \n");
	}

	if( diy == 1 )
	{
		nline_first = nline_first_do ;
		nline_final = nline_final_do ;
		ncdp_first  = ncdp_first_do  ;
		ncdp_final  = ncdp_final_do  ;

		ncdp_number   = ncdp_final  - ncdp_first  + 1 ;
		nline_number  = nline_final - nline_first + 1 ;
		nmx           = ncdp_number  ;
		nmy           = nline_number ;
	}
	else
	{
		ncdp_number   = ncdp_final  - ncdp_first  + 1 ;
		nline_number  = nline_final - nline_first + 1 ;
		nmx          = ncdp_number ;
		nmy          = nline_number ;
	}

	if( myid == 0 )
		echopar(fn1, fn2, fn3, ntau, nline_first, nline_final, ncdp_first, ncdp_final);

/*====================================================================*/

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
/*===== Initializing the Velocity Field File.      =======================*/
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
	dvy           =  vpar[2] ;    // Warning: the unit of dt should be ms. //
	dvx           =  vpar[3] ;
	dvtau         =  vpar[4] ;

	if( myid == 0 ) printf(" Initialize file %s done.\n ifn2=%d\n",fn2, ifn2 );
*/
	int	nvline_first, nvline_final, nvline_step, nvline_number, nvmy;
	int	nvcdp_first,  nvcdp_final,  nvcdp_step,  nvcdp_number,  nvmx;
	int	nvtau ;
	float	dvx , dvy, dvtau ;
/*
	// open and init. stacking velocity file.
	ostat	= open_stk_velocity( fn2 );
	if ( 0 != ostat )
	{
		mystat	= 21;	// open stacking velocity file failed.
		goto	mainProgStop;
	}
*/

	vrmsGatherIndex	*vrmsfileidx	= calloc(ntr_vrmsfile, sizeof(vrmsGatherIndex));
/*
	istat	= init_stk_velocity( vrmsfileOffsetBeg, ntr_vrmsfile,
		&nvtau, &dvtau, &dvx, &dvy,
		&nvcdp_first, &nvcdp_final,
		&nvline_first, &nvline_final,
		vrmsfileidx);
	if ( 0 != istat )
	{
		mystat	= 22;	// read stacking velocity file failed.
		goto	mainProgStop;
	}
	if( myid == 0 ) printf(" Initialize file %s done.\n", fn2);
*/
	nvline_step    = 1;
	nvcdp_step     = 1;

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
/*====================================================================*/

/*===== Check the parameters of the input data file.==================*/
/*
	if( ntau > ns )
	{
		printf(" Tau-domain sampling number(ntau)  is %d\n", ntau);
		printf(" Seismic trace sampling number(ns) is %d\n", ns  );
		printf(" ntau <= ns is accepted.\n");
		return 1;
	}

	if( ns > nvtau )
	{
		printf(" Velocity Field Sampling Number(nvtau) Lowers Than Trace Sampling Number(ns).\n");
		printf(" Velocity Field Sampling Number(nvtau) is %d\n", nvtau);
		printf(" Trace Sampling Number(ns) is %d\n", ns);
		return 1;
	}

	if( ncdp_first < nvcdp_first || ncdp_final > nvcdp_final )
	{
		printf(" Velocity Field Not Match(InLine).\n");
		printf(" ncdp_first = %d\t, nvcdp_first = %d\n", ncdp_first, nvcdp_first);
		printf(" ncdp_final = %d\t, nvcdp_final = %d\n", ncdp_final, nvcdp_final);
		return 1;
	}

	if( nline_first < nvline_first || nline_final > nvline_final )
	{
		printf(" Velocity Field Not Match(CrossLine).\n");
		printf(" nline_first = %d\t, nvline_first = %d\n", nline_first, nvline_first);
		printf(" nline_final = %d\t, nvline_final = %d\n", nline_final, nvline_final);
		return 1;
	}

	// check velocity field sampling.//
	if( (int)(dmx+0.5) != (int)(dvx+0.5) || (int)(dmy+0.5) != (int)(dvy+0.5) )
	{
		printf(" Velocity Field Sampling Interval Not Match.\n");
		printf(" Velocity Field: dvx = %f\t, dvy = %f\n", dvx, dvy);
		printf(" CMP     Gather: dmx = %f\t, dmy = %f\n", dmx, dmy);
		return 1;
	}
*/
 
/*====================================================================*/
/*===== Check the parameters of 3d Offset Plane-Wave decompositon.====*/
	phrmax	= phrmax*1.0E-6 ;
	phrmin	= 0.0 ;
	dphr	= (phrmax-phrmin)/(nphr-1) ;
	if( myid == 0 )
	{
		printf("                                                      \n");
		printf(" ********  Parameters of 5d PlaneWave Data.  ******** \n");
		printf(" nphr    = %d\n", nphr        );
		printf(" dphr    = %6.4f(micro-s/m)\n", dphr*1.0E6   );
		printf(" phrmin  = %6.4f(micro-s/m)\n", phrmin*1.0E6 );
		printf(" phrmax  = %6.4f(micro-s/m)\n", phrmax*1.0E6 );
		printf(" nmx     = %d\n", nmx         );
		printf(" nmy     = %d\n", nmy         );
		printf("                                                      \n");
	}

/*======Initialize the 5d PlaneWade Data File.========================*/
/*
	int	ipwflag = 0 ;
	ipwflag = init_planewave_file(fn3, &ifn3, icpar, cpar, ntau, dt,
		nmx, dmx, ncdp_first, ncdp_final, ncdp_step, 0,
		nmy, dmy, nline_first, nline_final, nline_step, 0,
		nphr, dphr, 1, nphr, 1, phrmin);
	if( ipwflag != 0 )
	{
		printf("  gdbfcinitwriteregulardata_ %s error.\n", fn3);
		printf("  Program Exit; return 1; \n");
		return 1;
	}
	if( myid == 0 ) printf(" Initialize file %s done.\n ifn3=%d\n",fn3, ifn3);
*/
/*====================================================================*/
/*
	ostat	= open_pwd_gather( fn3, "wb" );
	if ( 0 != ostat )
	{
		mystat	= 31;	// open pwd gather file failed.
		goto	mainProgStop;
	}
	if( myid == 0 ) printf(" Initialize file %s done.\n", fn3);
*/

	int	cdp1, cdp2 , ncdp, ncdpg ;
	int	ngroup, igroup, dgroup ;
	long	memory_max , memory_need ;

	ncdp		= nmx ;
	memory_max	= 1024*1024*1024 ;
	memory_need	= (long)ntau*nphr*ncdp*sizeof(float) ;
	if( memory_need <= memory_max )
	{
		ngroup	= 1 ;
		dgroup	= ncdp ;
	}
	else
	{
		ngroup	= (int)(memory_need/memory_max) + 1 ;
		dgroup	= ncdp/ngroup + 1 ;
	}

/*======3D Array Allocation according to the given maximum memory.====*/
	float	***pdat3d	= NULL ;	//pdat3d[dgroup][nphr][ntau].
	float	***tcmp		= NULL ;	//tcmp[dgroup][ntr][ns].
	float	**pdat2d	= NULL ;	//pdat2d[nphr][ntau].
	float	**vrms		= NULL ;	//vrms[nmx][ns].
	float	**offset	= NULL ;	//offset[dgroup][ntr].

	pdat3d	= alloc3float( ntau, nphr, dgroup);
	tcmp	= alloc3float( ns, ntr, dgroup );
	pdat2d	= alloc2float( ntau, nphr);
	vrms	= alloc2float( ns, nmx );
	offset	= alloc2float( ntr, dgroup );

	if( pdat3d == NULL || pdat2d == NULL || tcmp == NULL || vrms == NULL || offset == NULL  )
	{
		printf("   Array Allocation Failed.\n");
		return 1 ;
	}

	memset( (void*)&pdat3d[0][0][0], 0, sizeof(float)*1l*ntau*nphr*dgroup );
	memset( (void*)&tcmp[0][0][0],   0, sizeof(float)*1l*ns*ntr*dgroup );
	memset( (void*)&pdat2d[0][0],    0, sizeof(float)*1l*ntau*nphr );
	memset( (void*)&vrms[0][0],      0, sizeof(float)*1l*ns*nmx );
	memset( (void*)&offset[0][0],    0, sizeof(float)*1l*ntr*dgroup );

/*====================================================================*/
	if( myid == 0 )
	{
		printf(" ncdp_first  = %4d \n", ncdp_first);
		printf(" ncdp_final  = %4d \n", ncdp_final);
		printf(" nphr        = %4d \n", nphr      );
		printf(" ncdp        = %4d \n", ncdp      );
		printf(" ntau        = %4d \n", ntau      );
		printf(" ns          = %4d \n", ns        );
		printf(" ntr         = %4d \n", ntr       );
		printf(" ngroup      = %4d \n", ngroup    );
		printf(" dgroup      = %4d \n", dgroup    );
		printf(" memory_max  = %6.2f (MB) \n",  memory_max/1024.0/1024.0);
		printf(" memory_need = %6.2f (MB) \n", memory_need/1024.0/1024.0);

		printf("                                            \n");
		printf("********************************************\n");
		printf("*                                          *\n");
		printf("*   3d PlaneWave Decomposition Begin!      *\n");
		printf("*                                          *\n");
		printf("********************************************\n");
		printf("                                            \n");
	}

/*===== Main Process of 3d PlaneWave Decomposition. ==================*/
	int	iline, icdp, imx, idx_g, ntr_gather, ilen ;
	int	itaus, itaue;
	int	it, iphr, itr ;
	int	NumSend, isend, itag, iline_do, iline_done ;
	char	processname[25];

	int	s_table[4], r_table[4];
	int	itask, ntask, itask_do, itask_done ;

	ntask	= nmy*ngroup ;

	if ( np == 1 )
//	Sequential Computing.
	{
		printf("*************************** Sequential Computing.\n");
		for( iline = nline_first ; iline <= nline_final ; ++ iline )
		{
			//read_vrms_of_current_line( vrms, iline, ncdp_first, ncdp_final, ns, ifn2);
			for ( int imx = 0 ; imx < nmx ; ++ imx)
			for ( int it  = 0 ; it  < ns  ; ++ it )
				vrms[imx][it]	= 2000;
		/*
			int	rvstat	= read_stk_velocity(vrmsfileOffsetBeg, vrmsfileidx,
				ntr_vrmsfile, nvtau, iline,
				ncdp_first, ncdp_final, ns, vrms);

			if ( 0 != rvstat )
			{
				mystat	= 22;	// read stacking velocity file failed.
				goto	mainProgStop;
			}
		*/

			printf("*************************** myid: %d\t, Line: %d\n", myid, iline);
			for ( igroup = 0 ; igroup < ngroup ; ++ igroup )
			{
				cdp1	= ncdp_first + igroup*dgroup ;
				cdp2	= cdp1 + dgroup - 1 ;
				if( cdp2 >= ncdp_final )
					cdp2 = ncdp_final ;
				printf(" igroup = %4d \t", igroup);
				printf(" cdp1   = %4d \t", cdp1  );
				printf(" cdp2   = %4d \n", cdp2  );

				zero3float(pdat3d, ntau, nphr, dgroup);
				zero3float(tcmp,   ns,   ntr,  dgroup);
				zero2float(offset, ntr, dgroup);

				itaus	= 20 ;
				itaue	= ntau ;
				//printf(" itaus   = %4d \n", itaus  );
				//printf(" itaue   = %4d \n", itaue  );
//OpenMP Multi-Threads Mode.
#pragma omp parallel for private(ntr_gather, imx, idx_g)
				for( icdp = cdp1 ; icdp <= cdp2 ; ++ icdp )
				{
					imx	= icdp - ncdp_first ;
					idx_g	= icdp - cdp1 ;
					printf(" icdp = %4d imx = %4d idx_g %4d \n", icdp, imx, idx_g);

					ntr_gather = 0 ;
					//read_cmp_gather_of_current_line_new(tcmp[idx_g], offset[idx_g],
					//	iline, icdp, ntr, &ntr_gather, ns, ifn1);
					//if( ntr_gather <= 0 )	continue ;
					int	rcstat	= read_cmp_gather( cmpfileOffsetBeg,
						cmpfileidx, ntr_cmpfile, ns,
						iline, icdp, ntr_gather_max,
						tcmp[idx_g], offset[idx_g], &ntr_gather);	
write_2d_float_wb(tcmp[idx_g], ntr_gather, ns, "cmp_check.dat");
					//continue ;
	
					if( 0 != rcstat )	continue ;
					if( ntr_gather <= 3 )	continue ;

					if( myid == 0 && icdp%10 == 0 ) 
						printf(" CDP: %d\t, Ntr: %d\t, Line: %d\n", icdp, ntr_gather, iline);

//					FengBo add code, scaning for nan-value and zeroing it. 2010-07-19.
					for( itr = 0 ; itr < ntr_gather ; ++ itr )
					{
						for( it  = 0 ; it  < ns ; ++ it  )
						{
							if( isnan(tcmp[idx_g][itr][it]) )
								tcmp[idx_g][itr][it] = 0.0 ;
						}
					}

					tauph_transform_3d_to_2d_robust(ns, ntau, nphr, ntr_gather, itaus, itaue,
						tcmp[idx_g], pdat3d[idx_g], vrms[imx], offset[idx_g], dt, phrmin, dphr);

				}
				if( nmx == 1 )
				{
					write_2d_float_wb(pdat3d[0], nphr, ntau, "./taupdat3d.dat");
					write_2d_float_wb(tcmp[0], ntr, ns, "./cmp3d.dat");
					write_1d_float_wb(vrms[0], ns, "./vrms.dat");
				}
				//write_3d_float_wb(pdat3d, dgroup, nphr, ntau, "./taupdat3d.dat");
				//write_3d_float_ab(pdat3d, dgroup, nphr, ntau, "taupdat3d.dat", igroup);
				//write_planewave_data4d_group(pdat3d, ntau, nphr, cdp1, cdp2, iline, ifn3);

			}// end loop for cdp group.

/*
			// write the plane-wave gather file.
			int	wstat	= write_pwd_gather(pdat3d, ntau, nphr, nmx, nline_first, iline);
			if ( 0 != wstat )
			{
				mystat	= 33;	// write pwd gather file failed.
				goto	mainProgStop;
			}
*/

		}// end loop for line.

//		Close the 3d cmp gather file(fn1).
		//gdbfcreadsegyfree_(&ifn1);
		fclose(cmpGatherFilefp);
		printf(" Close the CMP Gather: %s done.\n", fn1);

//		Close the 3d stacking velicty file(fn2).
		//ivflag = 0 ;
		//gdbfcclosecutslicedata_(&ifn2, &ivflag);
		//fclose(vrmsFilefp);
		printf(" Close the Velocity File: %s done.\n", fn2);

//		Close the 3d offset plane-wave file(fn3).
		//ipwflag = 0 ;
		//gdbfccloseregulardata_(&ifn3, &ipwflag);
		//fclose(pwdGatherFilefp);
		printf(" Close the Plane-Wave File: %s done.\n", fn3);
	}
/*
	else
//	Parallel Computing.
	{
		if ( myid == 0 )
//		Master Node.
		{
			printf("\tTask number: %6d\t Thread number: %4d \n", ntask, np);
			NumSend = 1 ;
			for ( itask = 1 ; itask <= MIN(np-1,ntask) ; ++ itask )
			{
				igroup = itask%ngroup ;
				if ( igroup == 0 )	igroup = ngroup ;
				iline = (itask-igroup)/ngroup + 1 ;

				iline = nline_first + iline - 1 ;
				cdp1  = ncdp_first + (igroup-1)*dgroup ;
				cdp2  = cdp1 + dgroup - 1 ;
				if( cdp2 >= ncdp_final )
					cdp2 = ncdp_final ;
				s_table[0] = iline ;
				s_table[1] = cdp1  ;
				s_table[2] = cdp2  ;
				s_table[3] = itask ;

				MPI_Send(s_table, 4, MPI_INT, itask, itask, MPI_COMM_WORLD);
				NumSend = NumSend + 1 ;
				printf("\t**** send2pid: %3d\t, Line: %4d\t cdp1: %4d\t cdpn: %4d\n", itask, iline, cdp1, cdp2);
			}
			for ( itask = 1 ; itask <= ntask ; ++ itask )
			{
				MPI_Recv(&itask_done, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
				isend = Status.MPI_SOURCE ;
				if ( NumSend <= ntask )
				{
					igroup = NumSend%ngroup ;
					if ( igroup == 0 )	igroup = ngroup ;
					iline = (NumSend-igroup)/ngroup + 1 ;

					iline = nline_first + iline - 1 ;
					cdp1  = ncdp_first + (igroup-1)*dgroup ;
					cdp2  = cdp1 + dgroup - 1 ;
					if( cdp2 >= ncdp_final )
						cdp2 = ncdp_final ;
					s_table[0] = iline ;
					s_table[1] = cdp1  ;
					s_table[2] = cdp2  ;
					s_table[3] = NumSend ;

					MPI_Send(s_table, 4, MPI_INT, isend, NumSend, MPI_COMM_WORLD);
					NumSend = NumSend + 1 ;
					printf("\t**** send2id: %3d\t, Line: %4d\t cdp1: %4d\t cdpn: %4d\n", isend, iline, cdp1, cdp2);
				}
				else
					MPI_Send(0, 0, MPI_INT, isend, 0, MPI_COMM_WORLD);
			}
		}
		else
//		Working Node.
		{
			if( myid <= ntask )
			while(1)
			{
				MPI_Recv(r_table, 4, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
				itag = Status.MPI_TAG;
				if( itag == 0 )	break ;

				iline = r_table[0] ;
				cdp1  = r_table[1] ;
				cdp2  = r_table[2] ;
				itask = r_table[3] ;

				zero3float(pdat3d, ntau, nphr, dgroup);
				zero3float(tcmp,   ns,   ntr,  dgroup);
				zero2float(offset, ntr, dgroup);

				//read_vrms_of_current_line( vrms, iline, ncdp_first, ncdp_final, ns, ifn2);
//FengBo modified for Tianjin data.
			read_vrms_of_current_line_tmp(velin, dt_in, nt_in, vrms, dt, ns,
				ncdp_first, ncdp_final);

				itaus	= 50 ;
				itaue	= ntau ;
//OpenMP Multi-Threads Mode.
#pragma omp parallel for private(ntr_gather, imx, idx_g)
				for( icdp = cdp1 ; icdp <= cdp2 ; ++ icdp )
				{
					imx	= icdp - ncdp_first ;
					idx_g	= icdp - cdp1 ;
					ntr_gather = 0 ;
					read_cmp_gather_of_current_line_new(tcmp[idx_g], offset[idx_g],
						iline, icdp, ntr, &ntr_gather, ns, ifn1);
					if( ntr_gather <= 0 )	continue ;
	
					tauph_transform_3d_to_2d_robust(ns, ntau, nphr, ntr_gather, itaus, itaue,
						tcmp[idx_g], pdat3d[idx_g], vrms[imx], offset[idx_g], dt, phrmin, dphr);
				}

				write_planewave_data4d_group(pdat3d, ntau, nphr, cdp1, cdp2, iline, ifn3);
				MPI_Send(&itask, 1, MPI_INT, 0, itag, MPI_COMM_WORLD);
			}
		}

//		Close the 3d cmp gather file(fn1).
		gdbfcreadsegyfree_(&ifn1);
		if( myid == 0 ) printf(" Close the CMP Gather: %s done.\n", fn1);

//		Close the 3d velicty field file(fn2).
		ivflag = 0 ;
		gdbfcclosecutslicedata_(&ifn2, &ivflag);
		if( myid == 0 ) printf(" Close the Velocity File: %s done.\n", fn2);

//		Slave nodes close fn3; Master node Build Index File of fn3.
		MPI_Barrier(MPI_COMM_WORLD);
		ipwflag = 0 ;
		gdbfccloseregulardata_(&ifn3, &ipwflag);
		if( myid == 0 ) printf(" Close the PlaneWave File: %s done.\n", fn3);
	}
*/

/*===== Free Arrays and Close the Files.  ============================*/
	free(cmpfileidx);
	free(vrmsfileidx);

	free3float(pdat3d);
	free3float(tcmp);
	free2float(pdat2d);
	free2float(vrms);
	free2float(offset);

/*===== 3d PlaneWave Decomposition Finished! =========================*/
	if( myid == 0 )
	{
		printf(" phrmax   = %6.4f(micro-s/m)\n", phrmax*1.0E6 );
		printf(" phrmin   = %6.4f(micro-s/m)\n", phrmin*1.0E6 );
		printf(" dphr     = %6.4f(micro-s/m)\n", dphr*1.0E6   );
		printf(" nphr     = %4d\n", nphr  );
		printf(" nmy      = %4d\n", nmy   );
		printf(" nmx      = %4d\n", nmx   );
		printf(" ntau     = %4d\n", ntau  );
		printf("                                           \n");
		printf("*******************************************\n");
		printf("*                                         *\n");
		printf("*   3d PlaneWave Decomposition Finished!  *\n");
		printf("*                                         *\n");
		printf("*******************************************\n");
		printf("                                           \n");
	}

mainProgStop:
	fprintf(stdout, "myid = %d  mystat = %d\n", myid, mystat);

	MPI_Finalize();

	return 0 ;
}

int	readpar(char *parfn, char *fn1, char *fn2, char *fn3, int *ntau, int *nphr, float *phrmax, int *diy,
		int *nline_first_do, int *nline_final_do, int *ncdp_first_do, int *ncdp_final_do)
{
	FILE *fp;
	if((fp=fopen(parfn,"r"))==NULL)
	{
		printf("Can Not Open Paramter file: %s For Reading!\n", parfn);
		return 1;
	}
	fscanf(fp,"%s",fn1);
	fscanf(fp,"%s",fn2);
	fscanf(fp,"%s",fn3);
	fscanf(fp,"%d",ntau);
	fscanf(fp,"%d",nphr);
	fscanf(fp,"%f",phrmax);
	fscanf(fp,"%d",diy);

	if( *diy == 1 )
	{
		fscanf(fp,"%d",nline_first_do);
		fscanf(fp,"%d",nline_final_do);
		fscanf(fp,"%d",ncdp_first_do);
		fscanf(fp,"%d",ncdp_final_do);
	}
	fclose(fp);

	return 0;
}

void	echopar(char *fn1, char *fn2, char *fn3,
		int ntau, int nline_first, int nline_final, int ncdp_first, int ncdp_final)
{
	printf(" The CMP     gather filename is :%s\n", fn1 );
	printf(" The Velocity Field filename is :%s\n", fn2 );
	printf(" The PlaneWade Data filename is :%s\n", fn3 );
	printf(" Tau-domain sampling number is %d\n",  ntau );
	printf(" The First Line Number is %d\n", nline_first);
	printf(" The Last  Line Number is %d\n", nline_final);
	printf(" The First cdp  Number is %d\n", ncdp_first );
	printf(" The Last  cdp  Number is %d\n", ncdp_final );
}

/*
int	read_vrms_of_current_line( float **vrms, int iline, int ncdp_first, int ncdp_final, int ns, int ifn2)
{
	int ivflag, ivcdp, n4, its, ite, ix ;
	n4  = 1 ;
	its = 1 ;
	ite = ns;
	for( ivcdp = ncdp_first ; ivcdp <= ncdp_final ; ++ivcdp )
	{
		ivflag = 1 ;
		ix     = ivcdp - ncdp_first ;
		gdbfcreadcutslicedata_(&ifn2, &n4, &iline, &ivcdp, &its, &ite, &ivflag, &vrms[ix][0]);
		if( ivflag != 0 )
		{
			printf(" gdbfcreadcutslicedata_ error.\n File Number: %d\n", ifn2);
			return 1;
		}
	}
	return 0;
}

int	read_vrms_of_current_line_tmp( float **velin, float dt_in, int nt_in,
		float **vrms, float dt_out, int ns, 
		int ncdp_first, int ncdp_final)
{
printf("nt_in=%d",nt_in);
printf("dt_in=%f",dt_in);
printf("ns=%d",ns);
printf("dt_out=%f",dt_out);
	// dt_in and dt_out is "ms".
	int	ivcdp1	= 5363;
	int	ivcdp2	= 9243;
	for( int ivcdp = ncdp_first ; ivcdp <= ncdp_final ; ++ivcdp )
	{
		int	ixout   = ivcdp - ncdp_first ;
		int	ixin	= ncdp_first - ivcdp1;
		if (ixin < 0 )
			ixin	= 0;

		float	tmax	= nt_in * dt_in;
		for (int it = 0; it < ns; ++ it)
		{
			float	tt	= dt_out*it;
			if ( tt > tmax )
				vrms[ixout][it] = velin[ixin][nt_in-1];
			else
				vrms[ixout][it]	= sinc_interpolation(velin[ixin], dt_in, nt_in, tt);
		}
	}

	return 0;
}

int	read_cmp_gather_of_current_line(float **tcmp, float *offset,
		int iline, int icdp, int ntr, int *ntr_gather, int ns, int ifn1)
{
	int   **headpar=NULL ;
	int	itr, ntr1;
	float	sx, gx, sy, gy ;

	headpar = alloc2int( 4, ntr );
	zero2int(   headpar, 4, ntr);

	//zero2float( tcmp,   ns, ntr);
	//zero1float( offset, ntr );

	gdbfcreadsegygatherdatafloat_(&ifn1, &iline, &icdp, &ns, &ntr1, (int*)&headpar[0][0], (float*)&tcmp[0][0]);
	*ntr_gather = ntr1 ;
	if( ntr1 <= 0 ) return 1;

	for( itr = 0 ; itr < ntr1 ; ++ itr )
	{
		sx = headpar[itr][0];
		sy = headpar[itr][1];
		gx = headpar[itr][2];
		gy = headpar[itr][3];
		offset[itr] = sqrt( (gx-sx)*(gx-sx)+(gy-sy)*(gy-sy) );
//		printf("%8d\t%8d\t%8d\t%8d\n", (int)sx, (int)sy, (int)gx, (int)gy);
	}

	free2int(headpar);
	return 0;
}

int	read_cmp_gather_of_current_line_new(float **tcmp, float *offset,
		int iline, int icdp, int ntr, int *ntr_gather, int ns, int ifn1)
{
	int	itr, ntr1;
	int	ithead1[20]={0}, iflag1, its, ite, itrace ;
	float	sx, gx, sy, gy ;

	//zero2float( tcmp,   ns, ntr);
	//zero1float( offset, ntr );

	its    = 1  ;
	ite    = ns ;
	ntr1   = 0  ;
	iflag1 = 0  ;
	for( itr = 0 ; itr < ntr ; ++ itr )
	{
		itrace = itr + 1 ;
		gdbfcreadsegytrace_(&ifn1, &iline, &icdp, &itrace, &its, &ite, &iflag1, ithead1, (void*)&tcmp[itr][0]);
		if( iflag1 != 0 )
		{
			*ntr_gather = ntr1 ;
			return 1 ;
		}
		sx = ithead1[4] ;
		sy = ithead1[5] ;
		gx = ithead1[7] ;
		gy = ithead1[8] ;
		offset[itr] = sqrt( (gx-sx)*(gx-sx)+(gy-sy)*(gy-sy) );
		ntr1 ++ ;
//		printf("%8d\t%8d\t%8d\t%8d\n", (int)sx, (int)sy, (int)gx, (int)gy);
	}
	*ntr_gather = ntr1 ;

	return 0;
}

int	write_planewave_data4d( float **pdat, int ns, int ntau, int nphr, int icdp, int iline, int ifn3)
{
	// pdat4d[nphr][nmy][nmx][ns]. //
	int ipar240[20]={0} ;
	int iflag, iphr ;

	ipar240[0] = 1   ;        //! first grid number of 1th Dimension.
	//ipar240[1] = ns  ;        //! last  grid number of 1th Dimension.
	ipar240[1] = ntau  ;      //! last  grid number of 1th Dimension.
	ipar240[2] = icdp  ;      //! location of 2nd Dimension.
	ipar240[3] = iline ;      //! location of 3rd Dimension.

	for( iphr = 0 ; iphr < nphr ; ++ iphr )
	{
		ipar240[4] = iphr+1 ;      //! location of 4th Dimension.
		ipar240[5] = 1    ;        //! location of 5th Dimension.
		iflag = 1 ;
		gdbfcwriteregulardata_( &ifn3, ipar240, &iflag, (void *)&pdat[iphr][0]);
		if( iflag != 0 )
		{
			printf(" gdbfcwriteregulardata_ error.\n File Number: %d\n", ifn3);
			printf("ns = %d \n", ns);
			printf("icdp = %d \n", icdp);
			printf("iline = %d \n", iline);
			printf("iphr = %d \n", iphr);
			printf("iflag = %d \n", iflag);
			exit(-1);
		}
	}
	return 0;
}

int	write_planewave_data4d_group( float ***pdat3d, int ntau, int nphr, int cdp1, int cdp2, int iline, int ifn3)
{
	// pdat4d[nphr][nmy][nmx][ns]. //
	int ipar240[20]={0} ;
	int iflag, iphr, icdp, imx ;

//	FengBo updated IO functions, using gdbfcinitwriteregularfile_()
//		instead of gdbfcinitwriteregulardata_(). 2010-07-14.
// Input Paramater:
//  int *fno   File Mark
//  int *n4    Write File Point
//  int *flag Write Flag Mark
//       flag = 1 Write 1th Data
//            n4[0] = Start 1th Value
//            n4[1] = End   1th Value
//            n4[2] = Start 2th Value
//            n4[3] = Start 3th Value
//            n4[4] = Start 4th Value
//            n4[5] = Start 5th Value
//       flag = 2 Write 2th Data
//            n4[0] = Start 1th Value
//            n4[1] = Start 2th Value
//            n4[2] = End   2th Value
//            n4[3] = Start 3th Value
//            n4[4] = Start 4th Value
//            n4[5] = Start 5th Value
//       flag = 3 Write 2th Data
//            n4[0] = Start 1th Value
//            n4[1] = Start 2th Value
//            n4[2] = Start 3th Value
//            n4[3] = End   3th Value
//            n4[4] = Start 4th Value
//            n4[5] = Start 5th Value
//            .
//            .
//  Reaturn: flag = 0 Right
//               != 0 Error

	for( iphr = 0 ; iphr < nphr ; ++ iphr )
	{
		for( icdp = cdp1 ; icdp <= cdp2 ; ++ icdp )
		{
			iflag = 1 ;
			ipar240[0] = 1     ;       //! first grid number of 1th Dimension.
			ipar240[1] = ntau  ;       //! last  grid number of 1th Dimension.
			ipar240[2] = icdp  ;       //! location of 2nd Dimension.
			ipar240[3] = iline ;       //! location of 3rd Dimension.
			ipar240[4] = iphr+1 ;      //! location of 4th Dimension.
			ipar240[5] = 1    ;        //! location of 5th Dimension.
			imx = icdp - cdp1 ;
			//gdbfcwriteregulardata_( &ifn3, ipar240, &iflag, (void *)&pdat3d[imx][iphr][0]);
			gdbfcwriteregularfile_( &ifn3, ipar240, &iflag, (void *)&pdat3d[imx][iphr][0]);
			if( iflag != 0 )
			{
				printf(" gdbfcwriteregulardata_ error.\n File Number: %d\n", ifn3);
				printf("ntau  = %d \n", ntau );
				printf("icdp  = %d \n", icdp );
				printf("iline = %d \n", iline);
				printf("iphr  = %d \n", iphr );
				printf("iflag = %d \n", iflag);
				exit(-1);
			}
		}
	}
	return 0;
}

int	init_planewave_file(char *fn3, int *ifn3, int *icpar, float *cpar,
	int ntau, float dt, 
	int nmx,  float dmx,  int ncdp_first,  int ncdp_final,  int ncdp_step,  float x0min,
	int nmy,  float dmy,  int nline_first, int nline_final, int nline_step, float y0min,
	int nphr, float dphr, int iphr1, int iphr2, int iphrstep, float phrmin)
{

	icpar[0]  = 5 ;            //! Dimension Number
	icpar[1]  = 1 ;            //! Domain 1:time 2:depth
	icpar[2]  = 1 ;            //! Data Format 1:float 2:int 3:short...
	icpar[3]  = 1 ;            //! File Type 1:Seismic 2:Velocity 3:Other
	cpar[0]   = 0.0 ;          //! Value Minimum ( Write File Don't use)
	cpar[1]   = 0.0 ;          //! Value Maximum ( Write File Don't use)

	icpar[4]  = ntau ;         //! 1th Dimension grid number.
	icpar[5]  = 1  ;           //! first grid number of 1th Dimension.
	icpar[6]  = ntau ;         //! last  grid number of 1th Dimension.
	icpar[7]  = 1  ;           //! 1th Dimension grid step.
	cpar[2]   = 0.0;           //! origin coordinate of 1th Dimension grid.
	cpar[3]   = dt ;           //! sampling rate of 1th Dimension grid(ms).

	icpar[8]  = nmx ;          //! 2nd Dimension grid number.
	icpar[9]  = ncdp_first ;   //! first grid number of 2nd Dimension.
	icpar[10] = ncdp_final ;   //! last  grid number of 2nd Dimension.
	icpar[11] = ncdp_step  ;   //! 1th Dimension grid step. 
	cpar[4]   = x0min      ;   //! origin coordinate of 2nd Dimension grid.
	cpar[5]   = dmx ;          //! sampling rate of 2nd Dimension grid.

	icpar[12] = nmy ;          //! 3rd Dimension grid number.
	icpar[13] = nline_first ;  //! first grid number of 3rd Dimension.
	icpar[14] = nline_final ;  //! last  grid number of 3rd Dimension.
	icpar[15] = nline_step  ;  //! 3rd Dimension grid step.
	cpar[6]   = y0min       ;  //! origin coordinate of 3rd Dimension grid.
	cpar[7]   = dmy ;          //! sampling rate of 3rd Dimension grid.

	icpar[16] = nphr  ;        //! 4th Dimension grid number.
	icpar[17] = iphr1 ;        //! first grid number of 4th Dimension.
	icpar[18] = iphr2 ;        //! last  grid number of 4th Dimension.
	icpar[19] = iphrstep ;     //! 4th Dimension grid step.
	cpar[8]   = phrmin ;       //! origin coordinate of 4th Dimension grid.
	cpar[9]   = dphr ;         //! sampling rate of 4th Dimension grid.

	icpar[20] = 1   ;          //! 5th Dimension grid number.
	icpar[21] = 1   ;          //! first grid number of 5th Dimension.
	icpar[22] = 1   ;          //! last  grid number of 5th Dimension.
	icpar[23] = 1   ;          //! 5th Dimension grid step.
	cpar[10]  = 0.0 ;          //! origin coordinate of 5th Dimension grid.
	cpar[11]  = 1.0 ;          //! sampling rate of 5th Dimension grid.

	int	ipwflag = 0 ;
	gdbfcinitwriteregularfile_(fn3, icpar, cpar, &ipwflag, ifn3) ;
	//gdbfcinitwriteregulardata_(fn3, icpar, cpar, &ipwflag, ifn3) ;

	return ipwflag ;
}
*/
