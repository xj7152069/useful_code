#ifndef FB_RW_H
#define FB_RW_H

void read_1d_int_r     ( int *data , int n , char *fn ) ;
void read_1d_int_rb    ( int *data , int n , char *fn ) ;
void write_1d_int_w    ( int *data , int n , char *fn ) ;
void write_1d_int_wb   ( int *data , int n , char *fn ) ;
void write_1d_int_a    ( int *data , int n , char *fn , int flag ) ;
void write_1d_int_ab   ( int *data , int n , char *fn , int flag ) ;

void read_1d_float_r   ( float *data , int n , char *fn ) ;
void read_1d_float_rb  ( float *data , int n , char *fn ) ;
void write_1d_float_w  ( float *data , int n , char *fn ) ;
void write_1d_float_wb ( float *data , int n , char *fn ) ;
void write_1d_float_a  ( float *data , int n , char *fn , int flag ) ;
void write_1d_float_ab ( float *data , int n , char *fn , int flag ) ;

void read_2d_int_r     ( int **trace , int n1 ,int n2 , char *fn ) ;
void read_2d_int_rb    ( int **trace , int n1 ,int n2 , char *fn ) ;
void write_2d_int_w    ( int **trace , int n1 ,int n2 , char *fn ) ;
void write_2d_int_wb   ( int **trace , int n1 ,int n2 , char *fn ) ;
void write_2d_int_a    ( int **trace , int n1 ,int n2 , char *fn , int flag ) ;
void write_2d_int_ab   ( int **trace , int n1 ,int n2 , char *fn , int flag ) ;

void read_2d_float_r   ( float **trace , int n1 ,int n2 , char *fn ) ;
void read_2d_float_rb  ( float **trace , int n1 ,int n2 , char *fn ) ;
void write_2d_float_w  ( float **trace , int n1 ,int n2 , char *fn ) ;
void write_2d_float_wb ( float **trace , int n1 ,int n2 , char *fn ) ;
void write_2d_float_a  ( float **trace , int n1 ,int n2 , char *fn , int flag ) ;
void write_2d_float_ab ( float **trace , int n1 ,int n2 , char *fn , int flag ) ;

void read_3d_float_rb  ( float ***trace , int n1 , int n2 , int n3 , char *fn );
void write_3d_float_wb ( float ***trace , int n1 , int n2 , int n3 , char *fn );
void write_3d_float_ab ( float ***trace , int n1 , int n2 , int n3 , char *fn , int flag );

#endif
