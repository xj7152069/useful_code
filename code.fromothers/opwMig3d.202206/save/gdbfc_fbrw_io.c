#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int	gdbfcwrite3float_( char prefix[], float *data, int *n1f, int *n2f, int *n3f, int *idf )
{
	int   n1, n2, n3 , id;
	n1 = *n1f ;
	n2 = *n2f ;
	n3 = *n3f ;
	id = *idf ;

	long  filesize = (long)n1*n2*n3 ;
	char  fn[1000]="" ;
	char  suffix[8]="" ;
	char  idx[4]="000" ;

	if( id >= 1 && id < 10 )
	{
		char id00x[2];
		sprintf( id00x , "%1d" , id ) ;
		idx[2] = id00x[0] ;
	}
	else if ( id >= 10 && id < 100 )
	{
		char id0xx[3];
		sprintf( id0xx , "%2d" , id ) ;
		idx[1] = id0xx[0] ;
		idx[2] = id0xx[1] ;
	}
	else if ( id >= 100 && id < 1000 )
	{
		char idxxx[4];
		sprintf( idxxx , "%2d" , id ) ;
		idx[0] = idxxx[0] ;
		idx[1] = idxxx[1] ;
		idx[2] = idxxx[2] ;
	}
	else
	{
		strcpy ( idx , "bad" ) ;
	}
	strcpy ( suffix , idx ) ;
	strcat ( suffix , ".dat" ) ;
	strcat ( fn, prefix ) ;
	strcat ( fn, suffix ) ;
	printf(" The out filename is %s\n", fn );

	FILE  *fp = NULL;
	if( (fp=fopen(fn,"wb")) == NULL )
	{
		printf("\tcan not open file:%s to write.\n" , fn ) ;
		exit(0) ;
	}
	fwrite( &data[0] , sizeof(float) , filesize, fp ) ;
	fclose(fp);

#ifdef DEBUG_READ_WRITE_MODE
	printf("Succeed in writing file: %s\n",fn);
#endif

	return 0 ;
}

int	gdbfcwritecomplex_( char prefix[], float *data, long *numberf, int *idf,
	int *iofff, int *flag, int *rwflag )
{
	int   id   = *idf ;
	int   ioff = *iofff ; 

	long  filesize = *numberf;
	char  fn[1000]="" ;
	char  suffix[8]="" ;
	char  idx[4]="000" ;

	if( id >= 1 && id < 10 )
	{
		char id00x[2];
		sprintf( id00x , "%1d" , id ) ;
		idx[2] = id00x[0] ;
	}
	else if ( id >= 10 && id < 100 )
	{
		char id0xx[3];
		sprintf( id0xx , "%2d" , id ) ;
		idx[1] = id0xx[0] ;
		idx[2] = id0xx[1] ;
	}
	else if ( id >= 100 && id < 1000 )
	{
		char idxxx[4];
		sprintf( idxxx , "%2d" , id ) ;
		idx[0] = idxxx[0] ;
		idx[1] = idxxx[1] ;
		idx[2] = idxxx[2] ;
	}
	else
	{
		strcpy ( idx , "bad" ) ;
	}
	strcpy ( suffix , idx ) ;
	strcat ( suffix , ".dat" ) ;
	strcat ( fn, prefix ) ;
	strcat ( fn, suffix ) ;
	printf(" The out filename is %s\n", fn );

	FILE  *fp = NULL;
	if( *rwflag == 0 )
	{
		if( (fp=fopen(fn,"wb")) == NULL )
		{
			printf("\tcan not open file:%s to write.\n" , fn ) ;
			*flag = 1 ;
			return 1 ;
		}
	}
	else
	{
		if( (fp=fopen(fn,"rb+")) == NULL )
		{
			printf("\tcan not open file:%s to write.\n" , fn ) ;
			*flag = 1 ;
			return 1 ;
		}
	}
	long	fileoffset = filesize*(ioff-1)*sizeof(float)*2;
	fseek(fp, fileoffset, 0);

	*flag = 0 ;
	fwrite( &data[0] , sizeof(float)*2 , filesize, fp ) ;
	fclose(fp);

#ifdef DEBUG_READ_WRITE_MODE
	printf("Succeed in writing file: %s\n",fn);
#endif

//	printf("\tnumber = %d\n", (int)filesize);
	return 0 ;
}

int	gdbfcreadcomplex_( char prefix[], float *data, long *numberf, int *idf, int *iofff, int *flag )
{
	int   id   = *idf ;
	int   ioff = *iofff ; 

	long  filesize = *numberf;
	char  fn[1000]="" ;
	char  suffix[8]="" ;
	char  idx[4]="000" ;

	if( id >= 1 && id < 10 )
	{
		char id00x[2];
		sprintf( id00x , "%1d" , id ) ;
		idx[2] = id00x[0] ;
	}
	else if ( id >= 10 && id < 100 )
	{
		char id0xx[3];
		sprintf( id0xx , "%2d" , id ) ;
		idx[1] = id0xx[0] ;
		idx[2] = id0xx[1] ;
	}
	else if ( id >= 100 && id < 1000 )
	{
		char idxxx[4];
		sprintf( idxxx , "%2d" , id ) ;
		idx[0] = idxxx[0] ;
		idx[1] = idxxx[1] ;
		idx[2] = idxxx[2] ;
	}
	else
	{
		strcpy ( idx , "bad" ) ;
	}
	strcpy ( suffix , idx ) ;
	strcat ( suffix , ".dat" ) ;
	strcat ( fn, prefix ) ;
	strcat ( fn, suffix ) ;
	printf(" The out filename is %s\n", fn );

	FILE  *fp = NULL;
	if( (fp=fopen(fn,"rb")) == NULL )
	{
		printf("\tcan not open file:%s to read.\n" , fn ) ;
		*flag = 1 ;
		return 1 ;
	}
	long	fileoffset = filesize*(ioff-1)*sizeof(float)*2;
	fseek(fp, fileoffset, 0);

	*flag = 0 ;
	fread( &data[0] , sizeof(float)*2 , filesize, fp ) ;
	fclose(fp);

#ifdef DEBUG_READ_WRITE_MODE
	printf("Succeed in reading file: %s\n",fn);
#endif

//	printf("\tnumber = %d\n", (int)filesize);
	return 0 ;
}

