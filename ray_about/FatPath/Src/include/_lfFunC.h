/*************************************************************
 *
 * 简介：		动态内存开辟
 *
 * 作者：		罗飞
 * 完成日期：	2016.06.12
 *
 *
 **/

/*------------------------------------------------------------*/
void *alloc1 (size_t n1 /* fast dimension */, 
			  size_t size /* size of one element */);
/*< output-checking 1 allocation >*/


/*------------------------------------------------------------*/
void *realloc1 (void *v /* previous data */, 
				size_t n1 /* fast dimension */, 
				size_t size /* size of one element */);
/* output-checing re-allocate a 1-d array */


/*------------------------------------------------------------*/
void **alloc2 (size_t n1 /* fast dimension */, 
			   size_t n2 /* slow dimension */, 
			   size_t size /* size of one element */);
/*< output-checking 2 allocation >*/


/*------------------------------------------------------------*/
void ***alloc3 (size_t n1 /* fast dimension */, 
		 		size_t n2 /* slower dimension */, 
				size_t n3 /* slowest dimension */, 
				size_t size /* size of one element */);
/*< output-checking 3 allocation >*/


/*------------------------------------------------------------*/
void ****alloc4 (size_t n1 /* fast dimension */,
				 size_t n2 /* slower dimension */, 
				 size_t n3 /* slower dimension */, 
				 size_t n4 /* slowest dimension */, 
				 size_t size /* size of one element */);
/*< output-checking 4 allocation >*/


/*------------------------------------------------------------*/
void *****alloc5 (size_t n1 /* fast dimension */, 
				  size_t n2 /* slower dimension */, 
				  size_t n3 /* slower dimension */,
				  size_t n4 /* slower dimension */,
				  size_t n5 /* slowest dimension */,
				  size_t size /* size of one element */);
/*< output-checking 5 allocation >*/


/*------------------------------------------------------------*/
void ******alloc6 (size_t n1 /* fast dimension */,
				   size_t n2 /* slower dimension */,
				   size_t n3 /* slower dimension */, 
				   size_t n4 /* slower dimension */,
				   size_t n5 /* slower dimension */,
				   size_t n6 /* slowest dimension */, 
                   size_t size /* size of one element */);
/*< output-checking 6 allocation >*/


/*------------------------------------------------------------*/
void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free4 (void ****p);
void free5 (void *****p);
void free6 (void ******p);
/*< free 1~6D allocation >*/


/*------------------------------------------------------------*/
int *alloc1int (size_t n1 /* fast dimension */);
/*< int 1-D allocation >*/


/*------------------------------------------------------------*/
int *realloc1int (int *v /* previous data */, 
				  size_t n1 /* fast dimension */);
/*< re-allocate a int 1-D array >*/


/*------------------------------------------------------------*/
int **alloc2int (size_t n1 /* fast dimension */,
				 size_t n2 /* slow dimension */);
/*< int 2-D allocation >*/


/*------------------------------------------------------------*/
int ***alloc3int (size_t n1 /* fast dimension */,
				  size_t n2 /* slower dimension */, 
				  size_t n3 /* slowest dimension */);
/*< int 3-D allocation >*/


/*------------------------------------------------------------*/
void free1int (int *p);
void free2int (int **p);
void free3int (int ***p);
/*< free int 1~3D allocation >*/


/*------------------------------------------------------------*/
float *alloc1float (size_t n1 /* fast dimension */);
/*< float 1-D allocation >*/


/*------------------------------------------------------------*/
float *realloc1float (float *v /* previous data */, 
					  size_t n1 /* fast dimension */);
/*< re-allocate a float 1-D array >*/


/*------------------------------------------------------------*/
float **alloc2float (size_t n1 /* fast dimension */, 
					 size_t n2 /* slow dimension */);
/*< float 2-D allocation >*/


/*------------------------------------------------------------*/
float ***alloc3float (size_t n1 /* fast dimension */,
					  size_t n2 /* slower dimension */,
					  size_t n3 /* slowest dimension */);
/*< float 3-D allocation >*/


/*------------------------------------------------------------*/
void free1float (float *p);
void free2float (float **p);
void free3float (float ***p);
/*< free float 1~3D allocation >*/


/*------------------------------------------------------------*/
float ****alloc4float (size_t n1 /* fast dimension */,
					   size_t n2 /* slower dimension */,
					   size_t n3 /* slower dimension */,
					   size_t n4) /* slowest dimension */;
/*< float 4-D allocation >*/


/*------------------------------------------------------------*/
void free4float (float ****p);
/*< free float 4D allocation >*/


/*------------------------------------------------------------*/
float *****alloc5float (size_t n1 /* fast dimension */,
						size_t n2 /* slower dimension */,
						size_t n3 /* slower dimension */,
						size_t n4 /* slower dimension */,
						size_t n5 /* slowest dimension */);
/*< float 5-D allocation >*/


/*------------------------------------------------------------*/
void free5float (float *****p);
/*< free float 5D allocation >*/


/*------------------------------------------------------------*/
float ******alloc6float (size_t n1 /* fast dimension */,
						 size_t n2 /* slower dimension */,
						 size_t n3 /* slower dimension */,
						 size_t n4 /* slower dimension */,
						 size_t n5 /* slower dimension */,
						 size_t n6 /* slowest dimension */);
/*< float 6-D allocation >*/


/*------------------------------------------------------------*/
void free6float (float ******p);
/*< free float 6D allocation >*/


/*------------------------------------------------------------*/
int ****alloc4int (size_t n1 /* fast dimension */,
				   size_t n2 /* slower dimension */,
				   size_t n3 /* slower dimension */,
				   size_t n4 /* slowest dimension */);
/*< int 4-D allocation >*/


/*------------------------------------------------------------*/
void free4int (int ****p);
/*< free int 4D allocation >*/


/*------------------------------------------------------------*/
int *****alloc5int (size_t n1 /* fast dimension */, 
					size_t n2 /* slower dimension */,
					size_t n3 /* slower dimension */,
					size_t n4 /* slower dimension */,
					size_t n5 /* slowest dimension */);
/*< int 5-D allocation >*/


/*------------------------------------------------------------*/
void free5int (int *****p);
/*< free int 5D allocation >*/


/*------------------------------------------------------------*/
unsigned char ******alloc6uchar(size_t n1 /* fast dimension */,
								size_t n2 /* slower dimension */,
								size_t n3 /* slower dimension */,
								size_t n4 /* slower dimension */,
								size_t n5 /* slower dimension */,
								size_t n6 /* slowest dimension */);
/*< unsigned char 6-D allocation >*/


/*------------------------------------------------------------*/
unsigned char *****alloc5uchar(size_t n1 /* fast dimension */,
							   size_t n2 /* slower dimension */,
							   size_t n3 /* slower dimension */,
							   size_t n4 /* slower dimension */,
							   size_t n5 /* slowest dimension */);
/*< unsigned char 5-D allocation >*/


/*------------------------------------------------------------*/
void free6uchar(unsigned char ******p);
void free5uchar(unsigned char *****p);
/*< free unsigned char 6~5D allocation >*/


/*------------------------------------------------------------*/
unsigned short ******alloc6ushort(size_t n1 /* fast dimension */,
								  size_t n2 /* slower dimension */,
								  size_t n3 /* slower dimension */,
								  size_t n4 /* slower dimension */,
								  size_t n5 /* slower dimension */,
								  size_t n6 /* slowest dimension */);
/*< unsigned short 6-D allocation >*/


/*------------------------------------------------------------*/
unsigned short *****alloc5ushort(size_t n1 /* fast dimension */,
								 size_t n2 /* slower dimension */,
								 size_t n3 /* slower dimension */,
								 size_t n4 /* slower dimension */,
								 size_t n5 /* slowest dimension */);
/*< unsigned short 5-D allocation >*/


/*------------------------------------------------------------*/
unsigned short ***alloc3ushort(size_t n1 /* fast dimension */,
							   size_t n2 /* slower dimension */,
							   size_t n3 /* slowest dimension */);
/*< unsigned short 3-D allocation >*/


/*------------------------------------------------------------*/
unsigned short **alloc2ushort(size_t n1 /* fast dimension */,
							  size_t n2 /* slow dimension */);
/*< unsigned short 2-D allocation >*/


/*------------------------------------------------------------*/
void free2ushort(unsigned short **p);
void free3ushort(unsigned short ***p);
void free5ushort(unsigned short *****p);
void free6ushort(unsigned short ******p);
/*< free unsigned short 2~6D allocation >*/


/*------------------------------------------------------------*/
double *alloc1double (size_t n1 /* fast dimension */);
/*< double 1-D allocation >*/


/*------------------------------------------------------------*/
double *realloc1double (double *v, size_t n1);
/*< re-allocate a double 1-D array >*/


/*------------------------------------------------------------*/
double **alloc2double (size_t n1 /* fast dimension */, 
					   size_t n2 /* slow dimension */);
/*< double 2-D allocation >*/


/*------------------------------------------------------------*/
double ***alloc3double (size_t n1 /* fast dimension */,
						size_t n2 /* slower dimension */,
						size_t n3 /* slowest dimension */);
/*< double 3-D allocation >*/


/*------------------------------------------------------------*/
void free1double (double *p);
void free2double (double **p);
void free3double (double ***p);
/*< free double 1~3D allocation >*/


/*------------------------------------------------------------*/
void zero1int(int *p /* Input Data */, 
			  size_t n1 /* fast dimension */);
/*< initialize the 1-d int array with zero >*/


/*------------------------------------------------------------*/
void zero2int(int **p /* Input Data */,
			  size_t n1 /* fast dimension */,
			  size_t n2 /* slow dimension */);
/*< initialize the 2-d int array with zero >*/


/*------------------------------------------------------------*/
void zero3int(int ***p /* Input Data */, 
			  size_t n1 /* fast dimension */,
			  size_t n2 /* slower dimension */,
			  size_t n3 /* slowest dimension */);
/*< initialize the 3-d int array with zero >*/


/*------------------------------------------------------------*/
void zero2ushort(unsigned short **p /* Input Data */, 
				 size_t n1 /* fast dimension */,
				 size_t n2 /* slow dimension */);
/*< initialize the 2-d unsigned short array with zero >*/


/*------------------------------------------------------------*/
void zero3ushort(unsigned short ***p /* Input Data */,
				 size_t n1 /* fast dimension */,
				 size_t n2 /* slower dimension */,
				 size_t n3 /* slowest dimension */);
/*< initialize the 3-d unsigned short array with zero >*/


/*------------------------------------------------------------*/
void zero1float(float *p, size_t n1);
void zero2float(float **p, size_t n1, size_t n2);
void zero3float(float ***p, size_t n1, size_t n2, size_t n3);
void zero4float(float ****p, size_t n1, size_t n2, size_t n3, size_t n4);
/*< initialize the 1~4d float array with zero >*/


/*------------------------------------------------------------*/
void zero1double (double *p, size_t n1);
void zero2double (double **p, size_t n1, size_t n2);
void zero3double (double ***p, size_t n1, size_t n2, size_t n3);
/*< initialize the 1~3d double array with zero >*/


/*------------------------------------------------------------*/
float complex  *alloc1floatcomplex(size_t n1);
void free1floatcomplex(float complex *p);
/*< floatcomplex 1-D allocation >*/


/*------------------------------------------------------------*/
double complex *alloc1doublecomplex(size_t n1);
double complex *realloc1doublecomplex(double complex *v, size_t n1);
void free1doublecomplex(double complex *p);

float  complex **alloc2floatcomplex(size_t n1,size_t n2);
void free2floatcomplex(float complex **p);

double complex **alloc2doublecomplex(size_t n1, size_t n2);
void free2doublecomplex(double complex **p);

float complex ***alloc3floatcomplex(size_t n1,size_t n2, size_t n3);
void free3floatcomplex(float complex ***p);

double complex ***alloc3doublecomplex(size_t n1, size_t n2, size_t n3);
void free3doublecomplex(double complex ***p);

float complex ****alloc4floatcomplex(size_t n1,size_t n2, size_t n3, size_t n4);
void free4floatcomplex(float complex ****p);

double complex ****alloc4doublecomplex(size_t n1, size_t n2, size_t n3, size_t n4);
void free4doublecomplex(double complex ****p);

float complex *****alloc5floatcomplex(size_t n1, size_t n2, size_t n3, size_t n4,size_t n5);
void free5floatcomplex(float complex *****p);

double complex *****alloc5doublecomplex(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5doublecomplex(double complex *****p);
/*< double/float complex  1~6D allocation and free >*/


/*------------------------------------------------------------*/
void zero2floatcomplex(float complex **p,size_t n1,size_t n2);
void zero3floatcomplex(float complex ***p,size_t n1,size_t n2,size_t n3);
void zero4floatcomplex(float complex ****p,size_t n1,size_t n2,size_t n3,size_t n4);
void zero5floatcomplex(float complex *****p,size_t n1,size_t n2,size_t n3,size_t n4,size_t n5);
/*< initialize the 2~5d floatcomplex array with zero >*/


/*------------------------------------------------------------*/
void zero2doublecomplex(double complex **p, size_t n1, size_t n2);
void zero3doublecomplex(double complex ***p, size_t n1, size_t n2, size_t n3);
void zero4doublecomplex(double complex ****p, size_t n1, size_t n2, size_t n3, size_t n4);
void zero5doublecomplex(double complex *****p, size_t n1, size_t n2, size_t n3, size_t n4,size_t n5);
/*< initialize the 2~5d doublecomplex array with zero >*/

/*********************************************************************************************************/
