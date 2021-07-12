#include"../lib/RickerWavelet_func.h"

//double rickerwavelet_func(int Lw, float freq, float dt, float *signal)
double rickerwavelet_func(struct PARAMETER *param, float *signal)
{
    //雷克子波，freq为频率,f(t)为雷克子波
    int i;
    float t, t1;
    float temp[param->Lw];
    for(i=0;i<param->Lw;i++)
    {
        t=param->dt*i;
        t1=1.0/param->freq;//双边雷克子波
        //t1=0;//单边雷克子波
        temp[i]=50*(1-2*PI*PI*param->freq*param->freq*(t-t1)*(t-t1))
            *exp(-PI*PI*param->freq*param->freq*(t-t1)*(t-t1));
    }

    //导数形式
/*    signal[0]=0;*/
    //for (i=1;i<param->Lw;i++)
    //{
        //signal[i]=(temp[i]-temp[i-1])*1.0/param->dt;
/*    }*/

    //雷克子波形式
    for (i=0; i<param->Lw; i++)
    {
        signal[i]=temp[i];
    }



    if(param->wavelet_save_flag == 1)
    {
        //*雷克子波保存
        FILE *fp;
        if((fp = fopen ("../file/wavelet", "w"))!=NULL)
        {
            for (i=0;i<param->Lw;i++)
            {
                fprintf(fp,"%f\n",signal[i]);
            }
            fclose (fp);
        }
    }
    return 0.0;
}
