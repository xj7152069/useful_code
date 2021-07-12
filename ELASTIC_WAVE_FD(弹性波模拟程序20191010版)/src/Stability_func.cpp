#include"../lib/Stability_func.h"

double stability_func(struct PARAMETER *param,float **V_p, float **V_s)
{
    float V_pmax = 0;
    float V_smax = 0;

    for(int i=0; i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            if(V_pmax <= V_p[i][j]){V_pmax = V_p[i][j];}
            if(V_smax <= V_s[i][j]){V_smax = V_s[i][j];}
        }
    }

    param->V_pmax=V_pmax;

    float temp=param->dt * sqrt(
            pow(V_pmax*1.0/param->dx,2) + pow(V_smax*1.0/param->dz,2));

    if(temp >= 1.0)
    {
        printf("WARNING!! THE PARAMETER DOES NOT SATISFY THE STABILITY CONDITION!\n");
    }
    else
    {
        printf("GOOD!! THE PARAMETER SATISFY THE STABILITY CONDITION!\n");
    }

    return 0.0;
}

