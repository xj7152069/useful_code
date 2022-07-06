
#include "fbCommon.h"

void read_2d_float_rb_su ( float **trace , fbsegy_std *segy2d, int n1 , int n2 , char *fn )
{
	FILE *fp = NULL;
	int i , j ;
	if( (fp=fopen(fn,"rb")) == NULL )
	{
		printf("can not open file:%s to read.\n" , fn ) ;
		exit(1) ;
	}
	for( i = 0 ; i < n1 ; ++ i )
	{
		if ( 1  != fread( &segy2d[i] ,  sizeof(fbsegy_std) ,  1 , fp ) )
		{
			fprintf(stderr, " Read SU-header error.\n") ;
			exit(1) ;
		}
		if ( n2 != fread( &trace[i][0] , sizeof(float) ,  n2 , fp ) )
		{
			fprintf(stderr, " Read SU-data error.\n") ;
			exit(1) ;
		}
		/*
		int	gx	= segy2d[i].gx;
		int	sx	= segy2d[i].sx;
		int	scalco	= segy2d[i].scalco;
		printf("i=%d sx=%d gx=%d scalco=%d\n", i, sx, gx, scalco);
		*/
	}
	fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
	fprintf(stderr, "Succeed in reading file: %s\n",fn);
#endif
}

void write_2d_float_wb_suhdr ( float **trace , fbsegy_std *segy2d, int n1 , int n2 , char *fn )
{
	FILE *fp = NULL;
	int i , j ;
	if( (fp=fopen(fn,"wb")) == NULL )
	{
		printf("can not open file:%s to write.\n" , fn ) ;
		exit(1) ;
	}
	for( i = 0 ; i < n1 ; ++ i )
	{
		fwrite( &segy2d[i] ,   sizeof(fbsegy_std) ,  1 , fp ) ;
		fwrite( &trace[i][0] , sizeof(float) ,  n2 , fp ) ;
	}
	fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
	fprintf(stderr, "Succeed in writing file: %s\n",fn);
#endif
}

void write_3d_float_wb_suhdr ( float ***trace , fbsegy_std **segy2d, int n1 , int n2 , int n3, char *fn )
{
	FILE *fp = NULL;
	int i , j ;
	if( (fp=fopen(fn,"wb")) == NULL )
	{
		printf("can not open file:%s to write.\n" , fn ) ;
		exit(1) ;
	}
	for( i = 0 ; i < n1 ; ++ i )
	for( j = 0 ; j < n2 ; ++ j )
	{
		fwrite( &segy2d[i][j] ,   sizeof(fbsegy_std) ,  1 , fp ) ;
		fwrite( &trace[i][j][0] , sizeof(float) ,  n3 , fp ) ;
	}
	fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
	fprintf(stderr, "Succeed in writing file: %s\n",fn);
#endif
}
