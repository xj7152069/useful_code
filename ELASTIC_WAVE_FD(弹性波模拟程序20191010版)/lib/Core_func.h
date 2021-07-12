#ifndef CORE_FUNC_H
#define CORE_FUNC_H
#include "ShengShen_head.h"

//计算Vx,Vz
double core01_func( float **Vx_now,  float **Vx_aft,
					float **Vxx_now, float **Vxx_aft,
					float **Vxz_now, float **Vxz_aft,
					float **Vz_now,  float **Vz_aft,
					float **Vzx_now, float **Vzx_aft,
					float **Vzz_now, float **Vzz_aft,
					float **Txx_now, float **Txz_now,float **Tzz_now,
				    float wavelet, float **DEN,float **absorbx,float **absorbz,
					struct PARAMETER *param);

//计算应力
double core02_func(float **Txx_now, float **Txx_aft,
					float **Txxx_now, float **Txxx_aft,
					float **Txxz_now, float **Txxz_aft,
					float **Tzz_now, float **Tzz_aft,
					float **Tzzx_now, float **Tzzx_aft,
					float **Tzzz_now, float **Tzzz_aft,
					float **Txz_now, float **Txz_aft,
					float **Txzx_now, float **Txzx_aft,
					float **Txzz_now, float **Txzz_aft,
					float **Vx_now, float **Vz_now,
					float **lamda, float **miu,float wavelet,
				    float **DEN,float **absorbx,float **absorbz,struct PARAMETER *param);

#endif
