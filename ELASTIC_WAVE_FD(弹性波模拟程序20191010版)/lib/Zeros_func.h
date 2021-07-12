#ifndef ZEROS_FUNC_H
#define ZEROS_FUNC_H
#include "ShengShen_head.h"

//数组置零
double zeros_func(float **Vx_aft, float **Vxx_aft, float **Vxz_aft,
        float **Vx_now, float **Vxx_now, float **Vxz_now,
        float **Vz_aft, float **Vzx_aft, float **Vzz_aft,
        float **Vz_now, float **Vzx_now, float **Vzz_now,
        float **Txx_aft, float **Txxx_aft, float **Txxz_aft,
        float **Txx_now, float **Txxx_now,float **Txxz_now,
        float **Tzz_aft, float **Tzzx_aft, float **Tzzz_aft,
        float **Tzz_now, float **Tzzx_now, float **Tzzz_now,
        float **Txz_aft, float **Txzx_aft, float **Txzz_aft,
        float **Txz_now, float **Txzx_now, float **Txzz_now,
        float **record,struct PARAMETER *param);



#endif
