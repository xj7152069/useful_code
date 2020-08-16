/*
* test_lsqr.h
*
* 08 Sep 1999: First version from James W. Howse <jhowse@lanl.gov>
*/

/*------------------*/
/* Type definitions */
/*------------------*/

typedef struct PROD_DATA {
  dvec  *d_vec;
  dvec  *hy_vec;
  dvec  *hz_vec;
  dvec  *work_vec;
} prod_data;

/*---------------------*/
/* Function prototypes */
/*---------------------*/

void test_lsqr( long, long, long, long, double );
void lstp( long, long, double, dvec *, dvec *, prod_data *,
           double *, double * );
void test_mult( long, dvec *, dvec *, void * );
void lstp_mult( dvec *, dvec *, dvec * );
