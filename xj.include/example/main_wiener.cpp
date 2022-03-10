/*
    wiener filter: suit for 2d or 3d data
*/
//#define ARMA_DONT_USE_BLAS
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "./include/xjc.h"
using namespace std;
using namespace arma;

int main()
{
    int i,j,k,numi,narxmax;
    struct wiener3d par;           //结构体，包含数据和参数

    //int nz(1024),nx(200),ny(1);     //数据体采样点数
    //int nz(1024),nx(200),ny(200);   //数据体采样点数
    int nz(3501),nx(199),ny(14);    //数据体采样点数
    //int nz(3000),nx(974),ny(1);     //数据体采样点数
    
    wiener3d_parset(nx,ny,nz,par); //结构体初始化函数,并设置一组初始化参数
    par.nwx=21;par.nwy=9;          //数据窗求解误差逼近泛函
    par.narx=5;par.nary=3;         //wiener滤波器系数阶数
    par.dig_n=0.001;               //法方程求解对角稳定系数
    par.nmovex=par.nwx;            //滑动滤波的滑动步长，
    par.nmovey=par.nwy;            //这里设置为等于数据窗大小，避免多次覆盖滤波
    par.covermax=1;
    par.nf1=0,par.nf2=par.nf/4;    //滤波频率范围
    par.ncpu=6;
    numi=25;
    par.use_way=1;
    par.normal_bool=false;
    wiener3d_parupdate(par);       //调整参数后更新结构体

    ifstream inf,inf2;
    ofstream outf,outf2;
    //inf.open("./data/txz2.bin");
    //inf.open("./agcdata.3000.974.bin");
    inf.open("./data/data.3501.199.14.bin");
    //inf.open("./data3d.bin");

//注意处理实际数据的时候，可能一个测线包含的道数不相等;
//进行三维处理时需要重新调整;
    int nh(nx*ny);
    fmat basedata(nz,nh);
//读取数据并加入随机噪声，实际数据可以不用加
    dataread(basedata,nz,nh,inf);
    inf.close();
    //////////////
    for(j=0;j<nz;j++)
    {
        for(k=0;k<nh;k++)
        {
            basedata(j,k)+=0.001*(rand()%100)-0.05;
        }
    } 
    
    fmat dataline(nz,nx,fill::zeros);
    fmat filterdata(nz,nh,fill::zeros),dataline_t(nx,nz);
    k=0;
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            dataline.col(j)=basedata.col(k);
            k++;
        }
        par.data.col(i)=dataline.t();
        dataline.fill(0.0);
    }
    getagc_wiener3d(par);
    k=0;
    for(i=0;i<ny;i++)
    {
        dataline_t=par.data.col(i);
        dataline=dataline_t.t();
        for(j=0;j<nx;j++)
        {
            filterdata.col(k)=dataline.col(j);
            k++;
        }
        dataline_t.fill(0.0);
    }
    datawrite(filterdata,nz,nh,"./result/data.pretreatment.bin");

par.use_way=1;
{
    wiener3d_mid_thread(par,1); 
    wiener3d_cleardata(par);
}
/*
else
{
    wiener3d_mid_thread(par,numi);
    wiener3d_cleardata(par);
}*/

    outf.open("./result/data.wiener.bin");
    outf2.open("./result/err.wiener.bin");
    k=0;
    for(i=0;i<ny;i++)
    {
        dataline_t=par.realrebuildtx.col(i);
        dataline=dataline_t.t();
        for(j=0;j<nx;j++)
        {
            filterdata.col(k)=dataline.col(j);
            k++;
        }
        dataline_t.fill(0.0);
    }
    if(par.normal_bool)
    {
        filterdata=filterdata*par.normal_con;
    }
    datawrite(filterdata,nz,nh,outf);
    datawrite(filterdata=(filterdata-basedata),\
        nz,nh,outf2);
    outf.close();outf2.close();
    
    cout<<endl;
    cout<<"win_x = "<<par.nwx<<endl;
    cout<<"win_y = "<<par.nwy<<endl;
    cout<<"nar_x = "<<par.narx<<endl;
    cout<<"nar_y = "<<par.nary;
    
    return 0;
}

/*
struct wiener3d
{
    int nz,nx,ny,nf,nf1,nf2,nwx,nwy,narx,nary,\
        halfnarx,halfnary,nbackx,nbacky;
    float dz,dy,dx,df,dig_n;
    fcube data,realrebuildtx;
    cx_fcube rebuildfx,rebuildtx,datafx;
};
*/
