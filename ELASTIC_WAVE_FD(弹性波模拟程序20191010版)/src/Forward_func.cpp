#include "../lib/Forward_func.h"

double forward_func(struct PARAMETER *param)
{

    //---------------------------------------------//
    //                文件清除                     //
    //---------------------------------------------//

    remove("../file/record");
    remove("../file/Vx_time");
    remove("../file/Vz_time");
    remove("../file/p_time");
    remove("../file/Model");
    remove("../file/xPML");
    remove("../file/zPML");
    remove("../file/wavelet");
    //---------------------------------------------//
    //              变量申请及赋值                 //
    //---------------------------------------------//

    float **Vx_aft, **Vxx_aft, **Vxz_aft;
    float **Vx_now, **Vxx_now, **Vxz_now;
    float **Vz_aft, **Vzx_aft, **Vzz_aft;
    float **Vz_now, **Vzx_now, **Vzz_now;
    float **Txx_aft, **Txxx_aft, **Txxz_aft;
    float **Txx_now, **Txxx_now, **Txxz_now;
    float **Tzz_aft, **Tzzx_aft, **Tzzz_aft;
    float **Tzz_now, **Tzzx_now, **Tzzz_now;
    float **Txz_aft, **Txzx_aft, **Txzz_aft;
    float **Txz_now, **Txzx_now, **Txzz_now;
    float **DEN, **V_p,  **V_s, **V;
    float **lamda, **miu;
    float **absorbx, **absorbz, **record;
    float *signal;

    //*开辟动态数组*//
    Vx_now   = alloc2float(param->Nz,param->Nx);
    Vxx_now  = alloc2float(param->Nz,param->Nx);
    Vxz_now  = alloc2float(param->Nz,param->Nx);
    Vx_aft   = alloc2float(param->Nz,param->Nx);
    Vxx_aft  = alloc2float(param->Nz,param->Nx);
    Vxz_aft  = alloc2float(param->Nz,param->Nx);

    Vz_now   = alloc2float(param->Nz,param->Nx);
    Vzx_now  = alloc2float(param->Nz,param->Nx);
    Vzz_now  = alloc2float(param->Nz,param->Nx);
    Vz_aft   = alloc2float(param->Nz,param->Nx);
    Vzx_aft  = alloc2float(param->Nz,param->Nx);
    Vzz_aft  = alloc2float(param->Nz,param->Nx);

    Txx_aft  = alloc2float(param->Nz,param->Nx);
    Txxx_aft = alloc2float(param->Nz,param->Nx);
    Txxz_aft = alloc2float(param->Nz,param->Nx);
    Txx_now  = alloc2float(param->Nz,param->Nx);
    Txxx_now = alloc2float(param->Nz,param->Nx);
    Txxz_now = alloc2float(param->Nz,param->Nx);

    Txz_aft  = alloc2float(param->Nz,param->Nx);
    Txzx_aft = alloc2float(param->Nz,param->Nx);
    Txzz_aft = alloc2float(param->Nz,param->Nx);
    Txz_now  = alloc2float(param->Nz,param->Nx);
    Txzx_now = alloc2float(param->Nz,param->Nx);
    Txzz_now = alloc2float(param->Nz,param->Nx);

    Tzz_aft  = alloc2float(param->Nz,param->Nx);
    Tzzx_aft = alloc2float(param->Nz,param->Nx);
    Tzzz_aft = alloc2float(param->Nz,param->Nx);
    Tzz_now  = alloc2float(param->Nz,param->Nx);
    Tzzx_now = alloc2float(param->Nz,param->Nx);
    Tzzz_now = alloc2float(param->Nz,param->Nx);

    record   = alloc2float(param->Nt,param->Nr);
    DEN      = alloc2float(param->Nz,param->Nx);
    absorbx  = alloc2float(param->Nz,param->Nx);
    absorbz  = alloc2float(param->Nz,param->Nx);
    lamda    = alloc2float(param->Nz,param->Nx);
    miu      = alloc2float(param->Nz,param->Nx);
    V        = alloc2float(param->Nz,param->Nx);
    V_p      = alloc2float(param->Nz,param->Nx);
    V_s      = alloc2float(param->Nz,param->Nx);

    //------------------------------------//
    //            物理模型设计            //
    //------------------------------------//

    model_func(param,V, V_p, V_s, DEN ,lamda, miu);

    //------------------------------------//
    //            PML边界设计             //
    //------------------------------------//

    pml_func(param, absorbx, absorbz);

    //------------------------------------//
    //            子波设计                //
    //------------------------------------//

    signal = new float [param->Lw];
    rickerwavelet_func(param, signal);

    //-----------------------------------------//
    //          进度条显示设置                 //
    //-----------------------------------------//*/

    int i_timebar = 0;
    char timebar[25];
    const char *timebar_lable = "|/-\\";
    timebar[0] = 0;

    //-----------------------------------------//
    //          炮循环                         //
    //-----------------------------------------//
    for(int num_shot=0; num_shot<param->Ns; num_shot++)
    {
        //观测系统
        param->nx_location=param->Sx[num_shot];
        param->nz_location=param->Sz[num_shot];

        //初始化
        zeros_func(Vx_aft,Vxx_aft,Vxz_aft,
                Vx_now,Vxx_now,Vxz_now,
                Vz_aft,Vzx_aft,Vzz_aft,
                Vz_now,Vzx_now,Vzz_now,
                Txx_aft,Txxx_aft,Txxz_aft,
                Txx_now,Txxx_now,Txxz_now,
                Tzz_aft,Tzzx_aft,Tzzz_aft,
                Tzz_now,Tzzx_now,Tzzz_now,
                Txz_aft,Txzx_aft,Txzz_aft,
                Txz_now,Txzx_now,Txzz_now,
                record,param);

        //炮显示设置
        printf("No:%d/All:%d Shot\t",num_shot+1,param->Ns);
        printf("Position Of Source:[%fm,%fm]\n",(param->nx_location-param->PML)
                *param->dx,(param->nz_location-param->PML)*param->dz);
        //-----------------------------------------//
        //            波场时间递推计算             //
        //----------------------------------------//
        for(int k = 0; k < param->Nt; k++) //波场时间递推计算开始
        {

            //**进度条显示**//
            i_timebar=k*25/param->Nt;
            printf("[%-25s][%d%%][%c]\r", timebar, (i_timebar+1)*4, timebar_lable[k%4]);
            fflush(stdout);
            timebar[i_timebar] = '#';
            timebar[i_timebar+1] = 0;

            //*波场计算*//
            //计算速度
            core01_func(Vx_now, Vx_aft,
                    Vxx_now, Vxx_aft,
                    Vxz_now, Vxz_aft,
                    Vz_now, Vz_aft,
                    Vzx_now, Vzx_aft,
                    Vzz_now, Vzz_aft,
                    Txx_now,Txz_now,Tzz_now,
                    signal[k],DEN,absorbx,absorbz,
                    param);

            //计算应力
            core02_func(Txx_now, Txx_aft,
                    Txxx_now, Txxx_aft,
                    Txxz_now, Txxz_aft,
                    Tzz_now, Tzz_aft,
                    Tzzx_now, Tzzx_aft,
                    Tzzz_now, Tzzz_aft,
                    Txz_now, Txz_aft,
                    Txzx_now, Txzx_aft,
                    Txzz_now, Txzz_aft,
                    Vx_now, Vz_now,
                    lamda, miu,signal[k],
                    DEN,absorbx,absorbz,param);

            //-------------------------------------------//
            //               时间切片                    //
            //-------------------------------------------//

            if (param->time_vx_slice_save_flag == 1 && k%param->sdt == 0)
            {
                time_slice_IO_func(param, Vx_now, "../file/Vx_time");
            }
            if (param->time_vz_slice_save_flag == 1 && k%param->sdt == 0)
            {
                time_slice_IO_func(param, Vz_now, "../file/Vz_time");
            }
            if (param->time_p_slice_save_flag == 1 && k%param->sdt == 0)
            {
                time_slice_IO_func(param, Txx_now, "../file/p_time");
            }

            //-------------------------------------------//
            //               地震记录record              //
            //-------------------------------------------//
            for(int num_r = 0; num_r < param->Nr; num_r++)
            {
                int temp_num_rx=param->Rx[num_shot][num_r];
                int temp_num_rz=param->Rz[num_shot][num_r];

                if(param->acoustic_flag == 1){
                record[num_r][k] = Txx_now[temp_num_rx][temp_num_rz]
                                 + Tzz_now[temp_num_rx][temp_num_rz];
                }

                if(param->acoustic_flag == 0){
                record[num_r][k]=PS_IO_func(param, Vx_now, Vz_now, num_shot, num_r);}

            }
        }


        //-------------------------------------------//
        //       存储地表地震记录record              //
        //-------------------------------------------//

        if(param->record_save_flag == 1)
        {
            profile_IO_func(param, record);
        }
        printf("\n");
    }


    //---------------------------------------------//
    //                 释放内存                    //
    //---------------------------------------------//
    free2float(Vx_aft);
    free2float(Vxx_aft);
    free2float(Vxz_aft);
    free2float(Vx_now);
    free2float(Vxx_now);
    free2float(Vxz_now);
    free2float(Vz_aft);
    free2float(Vzx_aft);
    free2float(Vzz_aft);
    free2float(Vz_now);
    free2float(Vzx_now);
    free2float(Vzz_now);
    free2float(Txx_aft);
    free2float(Txxx_aft);
    free2float(Txxz_aft);
    free2float(Txx_now);
    free2float(Txxx_now);
    free2float(Txxz_now);
    free2float(Txz_aft);
    free2float(Txzx_aft);
    free2float(Txzz_aft);
    free2float(Txz_now);
    free2float(Txzx_now);
    free2float(Txzz_now);
    free2float(Tzz_aft);
    free2float(Tzzx_aft);
    free2float(Tzzz_aft);
    free2float(Tzz_now);
    free2float(Tzzx_now);
    free2float(Tzzz_now);
    free2float(record);
    free2float(DEN);
    free2float(absorbx);
    free2float(absorbz);
    free2float(lamda);
    free2float(miu);
    free2float(V);
    free2float(V_p);
    free2float(V_s);
    return 0.0;
}
