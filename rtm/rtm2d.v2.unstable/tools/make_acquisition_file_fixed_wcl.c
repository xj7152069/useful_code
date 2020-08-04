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
	nshot	= 3;
	ntrx = 100;
	ntry = 100;

	sx0	= 490.0;
	sy0	= 490.0;
	sz	= 0.0 ;
	dsx	= 100.0;
	dsy	= 100.0;

	gx0	= -10.0 ;
	gy0	= -10.0 ;
	gz	= 0.0 ;
	dgx	= 10.0;
	dgy	= 10.0;


	for ( int ishot = 0 ; ishot < nshot ; ++ ishot )
	{
		sx	= sx0 + (float)ishot * dsx ;
		sy	= sy0 + (float)ishot * dsy ;
		fprintf(fp, "%d\t%f\t%f\t%f\n", ishot, sx, sy, sz);


		for(itry =0; itry <ntry ; itry++)
		{
			gy = gy0 +  (float)itry * dgy ;

			for(itrx =0; itrx < ntrx ;itrx++)
			{
				gx = gx0  +  (float)itrx * dgx ;
				fprintf(fp, "%d\t%f\t%f\t%f\n", -1, gx, gy, gz);
			}
		}



	}

	fclose(fp);

	return 0;
}
