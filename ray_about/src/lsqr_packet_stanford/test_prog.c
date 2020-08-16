/*
* test_prog.c
* Tests C version of LSQR.
*
* 08 Sep 1999: First version from James W. Howse <jhowse@lanl.gov>
*/

#include "lsqr.h"

#include "test_lsqr.h"

void main()
{
  double  damp1,
          damp2,
          damp3,
          damp4;
  
  damp1 = 0.1;
  damp2 = 0.01;
  damp3 = 0.001;
  damp4 = 0.0001;
  
  test_lsqr(  1,  1, 1, 1, 0.0 );
  test_lsqr(  2,  1, 1, 1, 0.0 );
  test_lsqr( 40, 40, 4, 4, 0.0 );
  test_lsqr( 40, 40, 4, 4, damp2 );
  test_lsqr( 80, 40, 4, 4, damp2 );

  test_lsqr( 10, 10, 1, 8, 0.0 );
  test_lsqr( 40, 40, 4, 7, 0.0 );
  test_lsqr( 80, 40, 4, 6, 0.0 );
  test_lsqr( 20, 10, 1, 6, 0.0 );
  test_lsqr( 20, 10, 1, 8, 0.0 );

} /*  end main  */
