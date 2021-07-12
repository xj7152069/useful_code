#include "../lib/ShengShen_head.h"
#include "../lib/Parameter_func.h"
#include "../lib/Forward_func.h"
#include "../lib/RickerWavelet_func.h"
#include "../lib/IO_func.h"


int main()
{
    //---------------------------------------------//
    //            计时器--开始                     //
    //---------------------------------------------//
    double start = omp_get_wtime();

    //---------------------------------------------//
    //            正演模拟                         //
    //---------------------------------------------//
    //*各种参数*//
    PARAMETER param;
    parameter_func(&param);
    //**正演**//
    forward_func(&param);

    //---------------------------------------------//
    //            反演成像                         //
    //---------------------------------------------//



    //---------------------------------------------//
    //            计时器--结束                     //
    //---------------------------------------------//
    double end = omp_get_wtime();
    double totaltime=end-start;
    printf("It took %f seconds\n",totaltime);
    return 0;
}
