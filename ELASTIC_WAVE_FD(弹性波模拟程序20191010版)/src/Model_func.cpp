#include "../lib/Model_func.h"

double model_func(struct PARAMETER *param,float **V, float **V_p, float **V_s,
        float **DEN, float **lamda, float **miu)
{

    /********************模型赋值********************/
    param->acoustic_flag=1;
    for(int i=0;i<param->Nx;i++)
    {
        for(int j=0;j<param->Nz;j++)
        {
            V_p[i][j]   = 2000.0;
            DEN[i][j]   = 2400.0;
            miu[i][j]   = 0.0;
            lamda[i][j] = V_p[i][j]*V_p[i][j]*DEN[i][j]-2.0*miu[i][j];
            V_s[i][j]   = sqrt(miu[i][j]*1.0/DEN[i][j]);
            V[i][j]     = sqrt(V_p[i][j]*V_p[i][j]+V_s[i][j]*V_s[i][j]);
            if(miu[i][j]!=0){param->acoustic_flag = 0;}
        }

        for(int j=param->PML+100;j<param->Nz;j++)
        {
            V_p[i][j]   = 3000.0;
            DEN[i][j]   = 2400.0;
            miu[i][j]   = 0.0;
            lamda[i][j] = V_p[i][j]*V_p[i][j]*DEN[i][j]-2.0*miu[i][j];
            V_s[i][j]   = sqrt(miu[i][j]*1.0/DEN[i][j]);
            V[i][j]     = sqrt(V_p[i][j]*V_p[i][j]+V_s[i][j]*V_s[i][j]);
            if(miu[i][j]!=0){param->acoustic_flag = 0;}
        }
        for(int j=param->PML+150;j<param->Nz;j++)
        {
            V_p[i][j]   = 2500.0;
            DEN[i][j]   = 2400.0;
            miu[i][j]   = 0.0;
            lamda[i][j] = V_p[i][j]*V_p[i][j]*DEN[i][j]-2.0*miu[i][j];
            V_s[i][j]   = sqrt(miu[i][j]*1.0/DEN[i][j]);
            V[i][j]     = sqrt(V_p[i][j]*V_p[i][j]+V_s[i][j]*V_s[i][j]);
            if(miu[i][j]!=0){param->acoustic_flag = 0;}
        }
        for(int j=param->PML+300;j<param->Nz;j++)
        {
            V_p[i][j]   = 3500.0;
            DEN[i][j]   = 2400.0;
            miu[i][j]   = 0.0;
            lamda[i][j] = V_p[i][j]*V_p[i][j]*DEN[i][j]-2.0*miu[i][j];
            V_s[i][j]   = sqrt(miu[i][j]*1.0/DEN[i][j]);
            V[i][j]     = sqrt(V_p[i][j]*V_p[i][j]+V_s[i][j]*V_s[i][j]);
            if(miu[i][j]!=0){param->acoustic_flag = 0;}
        }
    }

    /********************模型加载********************/
/*
    FILE *fpvp;
	int Nx_temp,Nz_temp;

	Nx_temp          = 1601;
	Nz_temp          = 601;

	float **V_p_temp = new float*[Nx_temp];
	for(int i = 0; i < Nx_temp; i++)
	V_p_temp[i] = new float[Nz_temp];

    fpvp=fopen("/home/ss/WINDOWS/Model/sigsbee_1601x601.dat","rb");
    for(int i=0; i<1601; i++)
    {
        fread(V_p_temp[i],sizeof(float),601,fpvp);
    }
    fclose(fpvp);


    for (int i=0; i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
			V_p[i][j] = V_p_temp[i][j];
            DEN[i][j]   = 2400.0;
            miu[i][j]   = 0.0;
            lamda[i][j] = V_p[i][j]*V_p[i][j]*DEN[i][j]-2.0*miu[i][j];
            V_s[i][j]   = sqrt(miu[i][j]*1.0/DEN[i][j]);
            V[i][j]     = sqrt(V_p[i][j]*V_p[i][j]+V_s[i][j]*V_s[i][j]);
        }
    }
*/
    /********************稳定性判断********************/
    stability_func(param, V_p, V_s);

    /********************模型保存********************/
    if(param->model_save_flag == 1)
    {
        FILE *fp;
        if((fp = fopen ("../file/Model", "wb"))!=NULL)
        {

            for (int i=param->PML; i<param->Nx-param->PML; i++)
            {
                for (int j=param->PML; j<param->Nz-param->PML; j++)
                {
                    fwrite (&V_p[i][j] , sizeof(float), 1, fp);
                }
            }
            fclose (fp);
        }
    }
    return 0.0;
}
