#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define LengthMax 1024

int	main( int argc , char *argv[] )
{
	int	ntr, nshot ;
	float	sx0 , sy0, gx0, gy0, dsx, dsy, dgx, dgy, sx, sy, gx, gy, sz, gz ;
	int	itrx, itry;
	int ntrx, ntry;
	int nshotx, nshoty;
	int ishotx, ishoty;
	int ishot;
	char	fn[LengthMax]="";

	if ( 2 != argc )
	{
		fprintf(stderr, "\tParameter error!\n");
		return 1;
	}

	strcpy(fn, argv[1]);

	FILE	*fp = NULL ;
	if( NULL == (fp = fopen(fn, "w")) )
	{
		fprintf(stderr, "\tOpen file %s error!\n", fn);
		return 1;
	}

	// Initialization.
	// Multi-shot version.
//	nshot	= 3;
	nshotx	= 1;
	ntrx = 601;

	sx0	= 3000.0;
	sz	= 0.0 ;
	dsx	= 100.0;

	gx0	= 0.0 ;
	gz	= 0.0 ;
	dgx	= 10.0;


		for(ishotx =0 ; ishotx <nshotx ; ishotx ++)
		{

			sx	= sx0 + (float)ishotx * dsx ;
			ishot= ishoty*nshotx +ishotx+1;

			fprintf(fp, "%f\t%f\t%f\t%f\t%d\n", sz, sx, 0.0, 0.0, 0);


				for(itrx =0; itrx < ntrx ;itrx++)
				{
					gx = gx0  +  (float)itrx * dgx ;
					fprintf(fp, "%f\t%f\t%f\t%f\t%d\n", gz, gx, 0.0, 0.0, 1);
				}
		
		}



	fclose(fp);

	return 0;
}
