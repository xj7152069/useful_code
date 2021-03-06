/*
* test_lsqr.c
* For testing C version of LSQR.
*
* 08 Sep 1999: First version from James W. Howse <jhowse@lanl.gov>
* 02 Sep 2007: "max" redefined here to be "lsqr_max".
*              test_in->rel_mat_err changed from 1e-10 to 1e-15 and
*              test_in->rel_rhs_err changed from 1e-10 to 1e-15
*              to ensure "success" on the ill-conditioned rectangular
*              test cases.
*/

#include "lsqr.h"

#include "test_lsqr.h"

# define lsqr_max(a,b) ( (a) < (b) ? (b) : (a) )

void test_lsqr( long    num_rows,
		long    num_cols,
		long    duplicate_param,
		long    power_param,
		double  damp_param )
/*     
*     ------------------------------------------------------------------
*     This is an example driver routine for running LSQR.
*     It generates a test problem, solves it, and examines the results.
*     Note that subroutine APROD must be declared 'extern' if it is 
*     used only in the call to LSQR.
*     ------------------------------------------------------------------
*/
{
  double       dvec_norm2( dvec * );

  long         col_indx;
  
  double       unorm,
               enorm,
               etol,
               act_mat_cond_num,
               act_resid_norm;
  
  dvec         *act_rhs_vec,
               *act_sol_vec;

  lsqr_input   *test_in;
  lsqr_output  *test_out;
  lsqr_work    *test_work;
  lsqr_func    *test_func;
  prod_data    *test_prod;

  extern void  test_mult( long, dvec *, dvec *, void * );
  extern void  lstp_mult( dvec *, dvec *, dvec * );
  
  alloc_lsqr_mem( &test_in, &test_out, &test_work, &test_func, num_rows, 
                  num_cols );

  test_func->mat_vec_prod = test_mult;
  test_in->lsqr_fp_out = stdout;
/*
*     Allocate the the data structure of type 'prod_data'.
*/
  test_prod = (prod_data *) malloc( sizeof(prod_data) );

  test_prod->d_vec = (dvec *) alloc_dvec( num_cols );
  if (!test_prod->d_vec) lsqr_error("test_prog: vector \'d\' allocation failure "
"in function lsqr_test()", -1);

  test_prod->hy_vec = (dvec *) alloc_dvec( num_rows );
  if (!test_prod->hy_vec) lsqr_error("test_prog: vector \'hy\' allocation "
"failure in function lsqr_test()", -1);

  test_prod->hz_vec = (dvec *) alloc_dvec( num_cols );
  if (!test_prod->hz_vec) lsqr_error("test_prog: vector \'hz\' allocation "
"failure in function lsqr_test()", -1);

  test_prod->work_vec = (dvec *) alloc_dvec( lsqr_max( num_rows, num_cols ) );
  if (!test_prod->work_vec) lsqr_error("test_prog: vector \'work\' allocation "
"failure in function lsqr_test()", -1);
/*
*     Allocate the other needed vectors.
*/
  act_rhs_vec = (dvec *) alloc_dvec( num_rows );
  if (!act_rhs_vec) lsqr_error("test_prog: vector \'b\' allocation failure in "
"function lsqr_test()", -1);

  act_sol_vec = (dvec *) alloc_dvec( num_cols );
  if (!act_sol_vec) lsqr_error("test_prog: vector \'x_ans\' allocation failure "
"in function lsqr_test()", -1);
/*
*     Set the true solution for the test problem.
*/
  for( col_indx = 0; col_indx < num_cols; col_indx++ )
    act_sol_vec->elements[col_indx] = (double)( num_cols - (col_indx + 1) );
/*
*     Call the routine LSTP which generates the least squares test problem.
*/
  lstp( duplicate_param, power_param, damp_param, act_sol_vec, act_rhs_vec, 
        test_prod, &act_mat_cond_num, &act_resid_norm );
/*
*     Copy the right-hand side vector generated for the test problem by LSTP 
*     into the right-hand side vector for LSQR. 
*/
  dvec_copy( act_rhs_vec, test_in->rhs_vec );
/*
*  Set the initial guess for LSQR.
*/
  for( col_indx = 0; col_indx < num_cols; col_indx++ )
    test_in->sol_vec->elements[col_indx] = 0.0;
/*
*     Print information about the test problem.
*/
  fprintf( test_in->lsqr_fp_out,
	   
"--------------------------------------------------------------------\n" );
  fprintf( test_in->lsqr_fp_out,
    "Least-Squares Test Problem      P( %5li %5li %5li "
    "%5li %12.2e )\n\n",
	   num_rows, num_cols, duplicate_param, power_param, damp_param );
  fprintf( test_in->lsqr_fp_out,
    "Condition No. = %12.4e     Residual Function = "
    "%17.9e\n",
	   act_mat_cond_num, act_resid_norm );
  fprintf( test_in->lsqr_fp_out,
	   
"--------------------------------------------------------------------\n\n" );
/*
*     Set the input parameters for LSQR.
*/
  test_in->num_rows = num_rows;
  test_in->num_cols = num_cols;
  test_in->damp_val = damp_param;
  test_in->rel_mat_err = 1.0e-15;  /* Previously 1.0e-10; */
  test_in->rel_rhs_err = 1.0e-15;  /* Previously 1.0e-10; */
  test_in->cond_lim = 10.0 * act_mat_cond_num;
  test_in->max_iter = num_rows + num_cols + 50;
/*
*     Solve the test problem generated by LTSP by calling the routine LSQR.
*/
  lsqr( test_in, test_out, test_work, test_func, test_prod );

/*
*     Examine the results.
*     Print the residual norms RNORM and ARNORM given by LSQR, and then compute 
*     their true values in terms of the solution X obtained by  LSQR.
*     At least one of these norms should be small.
*/
  fprintf( test_in->lsqr_fp_out, "\t\t\tResidual Norm       Residual Norm       "
"Solution Norm\n" );
  fprintf( test_in->lsqr_fp_out, "\t\t\t||A x - b||_2   ||A^T A x - A^T b||_2   "
"||x - x0||_2\n\n" );
  fprintf( test_in->lsqr_fp_out, "Estimated by LSQR  %17.5e   %17.5e   "
"%17.5e\n",
	   test_out->resid_norm, test_out->mat_resid_norm, test_out->sol_norm );
/*
*     Compute  U = A*x - b.
*     This is the negative of the usual residual vector.  It will be close to 
*     zero only if b is a compatible rhs and x is close to a solution.
*/
  dvec_copy( act_rhs_vec, test_in->rhs_vec );
  dvec_scale( (-1.0), test_in->rhs_vec );
  test_mult( 0, test_out->sol_vec, test_in->rhs_vec, test_prod );
/*
*     Compute  v = A^T*u  +  DAMP**2 * x.
*     This will be close to zero in all cases if x is close to a solution.
*/
  dvec_copy( test_out->sol_vec, test_work->bidiag_wrk_vec );
  dvec_scale( ( lsqr_sqr(damp_param) ), test_work->bidiag_wrk_vec );
  test_mult( 1, test_work->bidiag_wrk_vec, test_in->rhs_vec, test_prod );
/*
*     Compute the norms associated with  x, u, v.
*/
  test_out->sol_norm = dvec_norm2( test_out->sol_vec );
  unorm = dvec_norm2( test_in->rhs_vec );
  test_out->resid_norm = sqrt( lsqr_sqr( unorm ) + lsqr_sqr( damp_param ) * 
			       lsqr_sqr( test_out->sol_norm ) ); 
  test_out->mat_resid_norm = dvec_norm2( test_work->bidiag_wrk_vec );
  fprintf( test_in->lsqr_fp_out, "Computed from X    %17.5e   %17.5e   "
"%17.5e\n\n",
	   test_out->resid_norm, test_out->mat_resid_norm, test_out->sol_norm );
/*
*     Print the solution and standard error estimates from LSQR.
*/
  fprintf( test_in->lsqr_fp_out, "Solution X\n" );
  for( col_indx = 0; col_indx < num_cols; col_indx++ )
    {
      fprintf( test_in->lsqr_fp_out, "%6li %14.6g", col_indx, 
	       test_out->sol_vec->elements[col_indx] );
      if( ( (col_indx + 1) % 4 ) == 0 )
	fprintf( test_in->lsqr_fp_out, "\n" );
    }

  fprintf( test_in->lsqr_fp_out, "\n\n" );

  fprintf( test_in->lsqr_fp_out, "Standard Errors SE\n" );
  for( col_indx = 0; col_indx < num_cols; col_indx++ )
    {
      fprintf( test_in->lsqr_fp_out, "%6li %14.6g", col_indx, 
	       test_out->std_err_vec->elements[col_indx] );
      if( ( (col_indx + 1) % 4 ) == 0 )
	fprintf( test_in->lsqr_fp_out, "\n" );
    }
  fprintf( test_in->lsqr_fp_out, "\n\n" );
/*
*     Print information about the accuracy of the LSQR solution.
*/
  for( col_indx = 0; col_indx < num_cols; col_indx++ )
    test_work->srch_dir_vec->elements[col_indx] = 
test_out->sol_vec->elements[col_indx]
      - act_sol_vec->elements[col_indx];

  enorm = dvec_norm2( test_work->srch_dir_vec ) /
                      ( 1.0 + dvec_norm2( act_sol_vec ) );
  etol = 1.0e-5;
  
  if( enorm <= etol )
    fprintf( test_in->lsqr_fp_out, "LSQR appears to be successful.    Relative "
"error in X = %10.2e\n\n",
	     enorm );
  if( enorm > etol )
    fprintf( test_in->lsqr_fp_out, "LSQR appears to have failed.      Relative "
"error in X = %10.2e\n\n",
	     enorm );
/*
*     Free the memory allocated for LSQR.
*/
  free_lsqr_mem( test_in, test_out, test_work, test_func );

  free_dvec( test_prod->d_vec );
  free_dvec( test_prod->hy_vec );
  free_dvec( test_prod->hz_vec );
  free_dvec( test_prod->work_vec );
  free( (prod_data *) (test_prod) );

  free_dvec( act_rhs_vec );
  free_dvec( act_sol_vec );

  return;
}

