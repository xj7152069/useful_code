
#include "fbCommon.h"

/*  create empty file. */
void create_empty_file( char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"w")) == NULL )
    {
        printf("can not create empty file:%s\n" , fn ) ;
        exit(0) ;
    }
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in creating empty file: %s\n" , fn ) ;
#endif
}


/* function for 1D case */
void read_1d_int_r  ( int *data , int n , char *fn )
{
    FILE *fp = NULL;
    int i;
    if( (fp=fopen(fn,"r")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n ; ++ i )
        fscanf( fp , "%d" , &data[i] ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif
}


void read_1d_int_rb ( int *data , int n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    fread( data , sizeof(int) , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n",fn);
#endif
}


void write_1d_int_w  ( int *data , int n , char *fn )
{
    FILE *fp = NULL;
    int i;
    if( (fp=fopen(fn,"w")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n ; ++ i )
        fprintf(fp,"%d\n",data[i]);
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_1d_int_wb ( int *data , int n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    fwrite( data , sizeof(int) , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_1d_int_a  ( int *data , int n , char *fn , int flag )
{
    FILE *fp = NULL;
    int i;
    if( flag == 0 )
        fp=fopen(fn,"w") ;
    else
        fp=fopen(fn,"a") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n ; ++ i )
        fprintf(fp,"%d\n",data[i]);
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


void write_1d_int_ab ( int *data , int n , char *fn , int flag )
{
    FILE *fp = NULL;
    if( flag == 0 )
        fp=fopen(fn,"wb") ;
    else
        fp=fopen(fn,"ab") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    fwrite( data , sizeof(int) , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


/* function for 1d case */
void read_1d_float_r  ( float *data , int n , char *fn )
{
    FILE *fp = NULL;
    int i;
    if( (fp=fopen(fn,"r")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n ; ++ i )
        fscanf( fp , "%f" , &data[i] ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif
}


void read_1d_float_rb ( float *data , int n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    fread( data , FLOAT_SIZE_BYTES , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n",fn);
#endif
}


void write_1d_float_w  ( float *data , int n , char *fn )
{
    FILE *fp = NULL;
    int i;
    if( (fp=fopen(fn,"w")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n ; ++ i )
        fprintf(fp,"%f\n",data[i]);
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_1d_float_wb ( float *data , int n , char *fn )
{
    FILE *fp = NULL;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    fwrite( data , FLOAT_SIZE_BYTES , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_1d_float_a  ( float *data , int n , char *fn , int flag )
{
    FILE *fp = NULL;
    int i;
    if( flag == 0 )
        fp=fopen(fn,"w") ;
    else
        fp=fopen(fn,"a") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n ; ++ i )
        fprintf(fp,"%f\n",data[i]);
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


void write_1d_float_ab ( float *data , int n , char *fn , int flag )
{
    FILE *fp = NULL;
    if( flag == 0 )
        fp=fopen(fn,"wb") ;
    else
        fp=fopen(fn,"ab") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    fwrite( data , FLOAT_SIZE_BYTES , n , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


/* function for 2d case */
void read_2d_int_r ( int **trace , int n1 ,int n2 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"r")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    for( j = 0 ; j < n2 ; ++ j )
        fscanf( fp , "%d" , &trace[i][j] ) ;
    fclose(fp) ;
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif

}


void read_2d_int_rb ( int **trace , int n1 , int n2 , char *fn )
{
    FILE *fp = NULL;
    int i;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
        fread( trace[i] , sizeof(int) , n2 , fp ) ;
    fclose(fp) ;
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif

}


void write_2d_int_w ( int **trace , int n1 , int n2 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"w")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    {
    for( j = 0 ; j < n2 ; ++ j )
        fprintf(fp,"%d\t",trace[i][j]);
    fprintf(fp,"\n");
    }
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_2d_int_wb ( int **trace , int n1 , int n2 , char *fn )
{
    FILE *fp = NULL;
    int i ;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
        fwrite( trace[i] , sizeof(int) , n2 , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_2d_int_a ( int **trace , int n1 , int n2 , char *fn , int flag )
{
    FILE *fp = NULL;
    int i , j ;
    if( flag == 0 )
        fp=fopen(fn,"w") ;
    else
        fp=fopen(fn,"a") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    {
    for( j = 0 ; j < n2 ; ++ j )
        fprintf(fp,"%d\t",trace[i][j]);
    fprintf(fp,"\n");
    }
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


void write_2d_int_ab ( int **trace , int n1 , int n2 , char *fn , int flag )
{
    FILE *fp = NULL;
    int i ;
    if( flag == 0 )
        fp=fopen(fn,"wb") ;
    else
        fp=fopen(fn,"ab") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
        fwrite( trace[i] , sizeof(int) , n2 , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


/* function for 2d case */
void read_2d_float_r ( float **trace , int n1 ,int n2 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"r")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    for( j = 0 ; j < n2 ; ++ j )
        fscanf( fp , "%f" , &trace[i][j] ) ;
    fclose(fp) ;
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif

}


void read_2d_float_rb ( float **trace , int n1 , int n2 , char *fn )
{
    FILE *fp = NULL;
    int i;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
        fread( trace[i] , FLOAT_SIZE_BYTES , n2 , fp ) ;
    fclose(fp) ;
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif

}


void write_2d_float_w ( float **trace , int n1 , int n2 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"w")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    {
    for( j = 0 ; j < n2 ; ++ j )
        fprintf(fp,"%f\t",trace[i][j]);
    fprintf(fp,"\n");
    }
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_2d_float_wb ( float **trace , int n1 , int n2 , char *fn )
{
    FILE *fp = NULL;
    int i ;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file:%s to write\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
        fwrite( trace[i] , FLOAT_SIZE_BYTES , n2 , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}


void write_2d_float_a ( float **trace , int n1 , int n2 , char *fn , int flag )
{
    FILE *fp = NULL;
    int i , j ;
    if( flag == 0 )
        fp=fopen(fn,"w") ;
    else
        fp=fopen(fn,"a") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    {
    for( j = 0 ; j < n2 ; ++ j )
        fprintf(fp,"%f\t",trace[i][j]);
    fprintf(fp,"\n");
    }
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}


void write_2d_float_ab ( float **trace , int n1 , int n2 , char *fn , int flag )
{
    FILE *fp = NULL;
    int i ;
    if( flag == 0 )
        fp=fopen(fn,"wb") ;
    else
        fp=fopen(fn,"ab") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
        fwrite( trace[i] , FLOAT_SIZE_BYTES , n2 , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}

/* function for 3d case */
void read_3d_float_rb ( float ***trace , int n1 , int n2 , int n3 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"rb")) == NULL )
    {
        printf("can not open file:%s to read\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    for( j = 0 ; j < n2 ; ++ j )
        fread( trace[i][j] , FLOAT_SIZE_BYTES , n3 , fp ) ;
    fclose(fp) ;
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in reading file: %s\n" , fn ) ;
#endif

}

void write_3d_float_wb ( float ***trace , int n1 , int n2 , int n3 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file:%s to writing\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    for( j = 0 ; j < n2 ; ++ j )
        fwrite( trace[i][j] , FLOAT_SIZE_BYTES , n3 , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}

void write_3d_float_ab ( float ***trace , int n1 , int n2 , int n3 , char *fn , int flag )
{
    FILE *fp = NULL;
    int i , j ;
    if( flag == 0 )
        fp=fopen(fn,"wb") ;
    else
        fp=fopen(fn,"ab") ;
    if( fp == NULL )
    {
        printf("can not open file:%s to add\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    for( j = 0 ; j < n2 ; ++ j )
        fwrite( trace[i][j] , FLOAT_SIZE_BYTES , n3 , fp ) ;
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in adding file: %s\n",fn);
#endif
}

void write_3d_float_wb_su ( float ***trace , int ***head, int n1 , int n2 , int n3 , char *fn )
{
    FILE *fp = NULL;
    int i , j ;
    if( (fp=fopen(fn,"wb")) == NULL )
    {
        printf("can not open file:%s to writing\n" , fn ) ;
        exit(0) ;
    }
    for( i = 0 ; i < n1 ; ++ i )
    for( j = 0 ; j < n2 ; ++ j )
    {
        fwrite( head[i][j] , sizeof(int) , 60 , fp ) ;
        fwrite( trace[i][j] , FLOAT_SIZE_BYTES , n3 , fp ) ;
    }
    fclose(fp);
#ifdef DEBUG_READ_WRITE_MODE
    printf("Succeed in writing file: %s\n",fn);
#endif
}

