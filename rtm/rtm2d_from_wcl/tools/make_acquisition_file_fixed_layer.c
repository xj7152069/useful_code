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
	int iflag;
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
//	nshotx	= 41;
	nshotx	= 61;
	ntrx = 1201;

//	sx0	= 1000.0;
	sx0	= 0.0;
	sy  = 0.0;
	sz	= 0.0;
	dsx	= 100.0;

	gx0	= 0.0 ;
	gy  = 0.0 ;
	gz	= 0.0 ;
	dgx	= 10.0;

	iflag = 2;


	//检波器固定，炮点移动
	if(iflag == 1)
	{
		for(ishotx =0 ; ishotx <nshotx ; ishotx ++)
		{

			sx	= sx0 + (float)ishotx * dsx ;
			ishot= ishotx+1;

			fprintf(fp, "%d\t%f\t%f\t%f\n", ishot, sx, sy, sz);

				for(itrx =0; itrx < ntrx ;itrx++)
				{
					gx = gx0  +  (float)itrx * dgx ;
					fprintf(fp, "%d\t%f\t%f\t%f\n", -1, gx, gy, gz);
				}
		
		}
	}

	//检波器位置随着炮点位置移动
	//双边观测
	if(iflag == 2)
	{
		for(ishotx =0 ; ishotx <nshotx ; ishotx ++)
		{

			sx	= sx0 + (float)ishotx * dsx ;
			ishot= ishotx+1;


			fprintf(fp, "%d\t%f\t%f\t%f\n", ishot, sx, sy, sz);


			for(itrx =0; itrx < ntrx ;itrx++)
			{
				gx = sx + dgx*(itrx- ntrx/2);
				fprintf(fp, "%d\t%f\t%f\t%f\n", -1, gx, gy, gz);
			}
		
		}
	}







	fclose(fp);

	return 0;
}
