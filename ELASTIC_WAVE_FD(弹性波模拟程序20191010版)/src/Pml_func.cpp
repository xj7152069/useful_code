#include"../lib/Pml_func.h"

double pml_sub_func(struct PARAMETER *param, float distance);

double pml_func(struct PARAMETER *param, float **absorbx, float **absorbz)
{
    int i,j;

    for(i=0; i<param->Nx; i++)
    {
        for(j=0; j<param->Nz;j++)
        {
            absorbx[i][j]=0.0;//吸收衰减函数初始值为0.0
            absorbz[i][j]=0.0;//吸收衰减函数初始值为0.0
        }
    }

    //左边界//
    for(i=0;i<param->PML;i++)
    {
        for(j=0;j<param->Nz;j++)
        {
            absorbx[i][j] = pml_sub_func(param, param->PML-i);
        }
    }
    //右边界//
    for(i=param->NX+param->PML;i<param->Nx;i++)
    {
        for(j=0;j<param->Nz;j++)
        {

            absorbx[i][j]=pml_sub_func(param, i - param->NX - param->PML);
        }
    }

    //上边界//
    for(i=0;i<param->Nx;i++)
    {
        for(j=0;j<param->PML;j++)
        {
            absorbz[i][j]=pml_sub_func(param, param->PML-j);
        }
    }
    //下边界//
    for(i=0;i<param->Nx;i++)
    {
        for(j=param->NZ+param->PML;j<param->Nz;j++)
        {
            absorbz[i][j]=pml_sub_func(param, j - param->NZ - param->PML);
        }
    }



    /********************PML保存********************/
    if (param->pml_save_flag == 1)
    {
        FILE *fpx;
        if((fpx = fopen ("../file/xPML", "w"))!=NULL)
        {

            for (i=0;i<param->Nx;i++)
            {
                for (j=0;j<param->Nz;j++)
                {
                    fwrite (&absorbx[i][j] , sizeof(float), 1, fpx);

                }
            }
            fclose (fpx);
        }

        FILE *fpz;
        if((fpz = fopen ("../file/zPML", "w"))!=NULL)
        {

            for (i=0;i<param->Nx;i++)
            {
                for (j=0;j<param->Nz;j++)
                {
                    fwrite (&absorbz[i][j] , sizeof(float), 1, fpz);

                }
            }
            fclose (fpz);
        }
    }
    return 0.0;
}

double pml_sub_func(struct PARAMETER *param, float distance)
{

    float temp;
    float V_max=param->V_pmax;

    //Wang
    //temp=- V_max *1.0 / param->PML  * log(0.000001)
    //*(0.25*distance*1.0/param->PML + 0.75*pow((distance)*1.0/param->PML,2));

    //Collino
    //temp=3.0*V_max*1.0/(2.0*param->PML)*log10(1.0/0.001)*
    //pow(distance*1.0/param->PML,2);

    //余弦衰减函数
    temp=500*(1.0-cos(PI*distance*1.0/(2.0*param->PML)));

    return temp;
}



















