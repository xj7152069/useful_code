#include"../lib/Parameter_func.h"


double parameter_func(struct PARAMETER* param)
{
    //**模型参数
    param -> NX     = 200;    //x方向样点数
    param -> NZ     = 200;    //z方向样点数
    param -> Nt     = 2000;   //模拟时间样点数
    param -> dx     = 10.0;   //x方向取样
    param -> dz     = 10.0;   //z方向取样
    param -> dt     = 0.0005; //时间取样
    param -> PML    = 25+8;   //PML边界厚度
    param -> V_pmax = 0;      //最大纵波速递


    //**有限差分系数
    param -> C1      = 1.2340911;
    param -> C2      = -1.0664985e-01;
    param -> C3      = 2.3036367e-02;
    param -> C4      = -5.3423856e-03;
    param -> C5      = 1.0772712e-03;
    param -> C6      = -1.6641888e-04;
    param -> C7      = 1.7021711e-005;
    param -> C8      = -8.5234642e-007;//差分系数

    //计算空间大小
    param -> Nx    = param -> NX+2*param -> PML;//x方向采样点数(x-z方向各留了PML行用于边界处理)
    param -> Nz    = param -> NZ+2*param -> PML;

    //子波参数设置
    param -> Lw    = param->Nt;//子波长度
    param -> freq  = 20;       //主频

    //**观测系统
    param->nx_location = param->PML;
    param->nz_location = param->PML;
    //炮点&检波点坐标系设置
    param->Ns    = 1;                      //炮数
    param->dsx   = param->NX/param->Ns;    //x方向炮间距
    param->dsz   = 0;                      //z方向炮间距

    param->Nr    = param->NX;              //检波点数目
    param->drx   = 1;                      //x方向检波器间距
    param->drz   = 0;                      //z方向检波器间距

    param->Sx=alloc1int(param->Ns);
    param->Sz=alloc1int(param->Ns);
    param->Rx=alloc2int(param->Nr,param->Ns);
    param->Rz=alloc2int(param->Nr,param->Ns);

    for (int i = 0; i<param->Ns; i++)
    {
        param->Sx[i] = param->PML + i*param->dsx; //炮点坐标Sx
        param->Sz[i] = param->PML + i*param->dsz; //炮点坐标Sz
        for (int j = 0; j<param->Nr; j++)
        {
            param->Rx[i][j] = param->PML + j*param->drx; //每一炮对应的检波点坐标Rx
            param->Rz[i][j] = param->Sz[i] + j*param->drz;//每一炮对应的检波点坐标Rz
        }
    }


    //**数据保存设置
    param->record_save_flag        = 1;  //1表示保存地震剖面，其他数字表示不保存
    param->time_vx_slice_save_flag = 0;  //1表示保存时间切片，其他数字表示不保存
    param->time_vz_slice_save_flag = 0;  //1表示保存时间切片，其他数字表示不保存
    param->time_p_slice_save_flag  = 0;  //1表示保存时间切片，其他数字表示不保存
    param->sdt                     = 20; //时间切片保存间隔
    param->model_save_flag         = 0;  //1表示输出速度模型，其他数字表示不保存
    param->pml_save_flag           = 0;  //1表示保存PML模型,其他数字不保存
    param->wavelet_save_flag       = 1;  //1表示保存子波，其他数字不保存

    //**并行设计
    param->omp_ncores = 8; //单进程 线程数目


    return 0.0;
}
