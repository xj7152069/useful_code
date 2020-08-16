/*
* test_mult.c
* For testing the matrix-vector product routine for the LSQR test problems.
*
* 08 Sep 1999: First version from James W. Howse <jhowse@lanl.gov>
*/

#include "lsqr.h"

#include "test_lsqr.h"

void test_mult( long  mode,
	        dvec  *x, 
                dvec  *y,
                void  *prod )
/*
*     ------------------------------------------------------------------
*     This is the matrix-vector product routine required by LSQR
*     for a test matrix of the form  A = HY*D*HZ.  The quantities
*     defining D, HY, HZ and the work array W are in the structure data.
*     ------------------------------------------------------------------
*/
{
  long       vec_indx,
             num_rows,
             num_cols;

  prod_data  *data;
  
  data = (prod_data *) prod;  

  num_cols = x->length;
  num_rows = y->length;
/*  
*     Compute  Y = Y + A*X
*/
  if( mode == 0 )
    {
      lstp_mult( data->hz_vec, x, data->work_vec );
      
      for( vec_indx = 0; vec_indx < num_cols; vec_indx++ )
	data->work_vec->elements[vec_indx] *= data->d_vec->elements[vec_indx];
      
      for( vec_indx = num_cols; vec_indx < num_rows; vec_indx++ )
	data->work_vec->elements[vec_indx] = 0.0;
      
      lstp_mult( data->hy_vec, data->work_vec, data->work_vec );

      for( vec_indx = 0; vec_indx < num_rows; vec_indx++ )
	y->elements[vec_indx] += data->work_vec->elements[vec_indx];
    }
/*  
*     Compute  X = X + A^T*Y
*/  
  if( mode == 1 )
    {
      lstp_mult( data->hy_vec, y, data->work_vec );
      
      for( vec_indx = 0; vec_indx < num_cols; vec_indx++ )
	data->work_vec->elements[vec_indx] *= data->d_vec->elements[vec_indx];
      
      lstp_mult( data->hz_vec, data->work_vec, data->work_vec );

      for( vec_indx = 0; vec_indx < num_cols; vec_indx++ )
	x->elements[vec_indx] += data->work_vec->elements[vec_indx];
    }
  
  return;
}


void lstp_mult( dvec  *h,
	        dvec  *x, 
                dvec  *y )
/*
*     ------------------------------------------------------------------
*     HPROD  applies a Householder transformation stored in H to get
*     Y = ( I - 2*HZ*HZ(transpose) ) * X.
*     ------------------------------------------------------------------
*/
{
  long    vec_indx,
          vec_length;

  double  dwork;
  
  vec_length = h->length;
  dwork = 0.0;

  for( vec_indx = 0; vec_indx < vec_length; vec_indx++ )
    dwork += h->elements[vec_indx] * x->elements[vec_indx];
  
  dwork += dwork;
  
  for( vec_indx = 0; vec_indx < vec_length; vec_indx++ )
    y->elements[vec_indx] = x->elements[vec_indx] - dwork * 
h->elements[vec_indx];
  
  return;
}
