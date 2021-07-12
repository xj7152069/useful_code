#include "../lib/Zeros_func.h"

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
			  float **record,struct PARAMETER *param)
{
	//*初始化*//
    for(int i=0; i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz;j++)
        {
            Vx_now[i][j] = 0.0;
            Vxx_now[i][j]= 0.0;
            Vxz_now[i][j]= 0.0;
            Vx_aft[i][j] = 0.0;
            Vxx_aft[i][j]= 0.0;
            Vxz_aft[i][j]= 0.0;

            Vz_now[i][j] = 0.0;
            Vzx_now[i][j]= 0.0;
            Vzz_now[i][j]= 0.0;
            Vz_aft[i][j] = 0.0;
            Vzx_aft[i][j]= 0.0;
            Vzz_aft[i][j]= 0.0;

            Txx_now[i][j] = 0.0;
            Txxx_now[i][j]= 0.0;
            Txxz_now[i][j]= 0.0;
            Txx_aft[i][j] = 0.0;
            Txxx_aft[i][j]= 0.0;
            Txxz_aft[i][j]= 0.0;

            Txz_now[i][j] = 0.0;
            Txzx_now[i][j]= 0.0;
            Txzz_now[i][j]= 0.0;
            Txz_aft[i][j] = 0.0;
            Txzx_aft[i][j]= 0.0;
            Txzz_aft[i][j]= 0.0;

            Tzz_now[i][j] = 0.0;
            Tzzx_now[i][j]= 0.0;
            Tzzz_now[i][j]= 0.0;
            Tzz_aft[i][j] = 0.0;
            Tzzx_aft[i][j]= 0.0;
            Tzzz_aft[i][j]= 0.0;
        }
    }


    for(int i=0; i<param->Nr; i++)
    {
        for(int j=0; j<param->Nt;j++)
        {
            record[i][j]=0.0;
        }
    }

return 0.0;
}
