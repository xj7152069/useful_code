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
void dataread_fromsegy_to3dblock(segyhead & head, int nline, int keyw);

int main()
{
    int i,j,k,numi,narxmax,ncpu(6);
    struct wiener3d par;           //结构体，包含数据和参数
    ofstream outf_globel;
    ifstream inf,inf2;
    ofstream outf,outf2;
    outf_globel.open("./result/data.wiener.bin");
    int nz(3501),nx(199),ny(14);    //数据体采样点数
    char file[999];
    file[0]='\0';
    strcat(file,"/media/xj/TOSHIBA EXT/seismic_data/TPTD3D_CMPGATHER_LINE1692_1711.SEGY");

    segyhead head;
    head.filename[0]='\0';
    strcat(head.filename,file);
    segyhead_open(head);
    segyhead_endianget(head);
    fmat data;
    for(k=0;k<150000;k++)
    {
        segyhead_readonetrace_tofmat(head,data);
    }
    int nline(21),keyw(1),num(0);

while(!head.infile.eof() && num<30)
{
    cout<<"now is running: "<<num<<endl;
    num++;
    dataread_fromsegy_to3dblock(head, nline, keyw);
    nx=head.nx,ny=nline,nz=head.nz;
    cout<<nx<<"|"<<nz<<endl;
    wiener3d_parset(nx,ny,nz,par); //结构体初始化函数,并设置一组初始化参数
    //par.nwx=21;par.nwy=9;          //数据窗求解误差逼近泛函
    //par.narx=5;par.nary=3;         //wiener滤波器系数阶数
    par.dig_n=0.001;               //法方程求解对角稳定系数
    par.dig_n2=0.001;
    par.nmovex=par.nwx;            //滑动滤波的滑动步长，
    par.nmovey=par.nwy;            //这里设置为等于数据窗大小，避免多次覆盖滤波
    par.covermax=1;
    par.nf1=0,par.nf2=par.nf/8;    //滤波频率范围
    par.ncpu=ncpu;
    wiener3d_parupdate(par);       //调整参数后更新结构体
//////////////////////////////////////////
    inf.open("./data.oneblock.bin");
//注意处理实际数据的时候，可能一个测线包含的道数不相等;
//进行三维处理时需要重新调整;
    fmat basedata, dataline, dataline_t;
    int nh(nx*ny);
    basedata.zeros(nz,nh);
    dataread(basedata,nz,nh,inf);
    inf.close();
///////////////////
    dataline.zeros(nz,nx);
    k=0;
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            dataline.col(j)=basedata.col(k);
            for(int k2=0;k2<nz;k2++)
            {
                dataline(k2,j)+=0.001*(rand()%100)-0.05;
            }
            k++;
        }
        par.data.col(i)=dataline.t();
    }
////////////////////////////////
{
    par.use_way=1;
    wiener3d_mid_thread(par,ncpu); 
    wiener3d_cleardata(par);
}
///////////////////////////////
    
    dataline_t.zeros(nx,nz);
    k=0;
    for(i=0;i<ny;i++)
    {
        dataline_t=par.realrebuildtx.col(i);
        dataline=dataline_t.t();
        for(j=0;j<nx;j++)
        {
            basedata.col(k)=dataline.col(j);
            k++;
        }
    }
    datawrite(basedata,nz,nh,outf_globel);
}
cout<<head.nx<<"|"<<head.nz<<endl;
    outf_globel.close();
    return 0;
}

void dataread_fromsegy_to3dblock(segyhead & head, int nline, int keyw)
{
    //keyw==1, use cdp; keyw==2, use trace;
    int i,j,k,nz(3501),ny(15),cdp0,cdp1,trac0,trac1;
    ny=nline;
    fmat data,numnx(ny,1);
    ofstream outf;
    outf.open("./data.bin");
    int nx(1),maxnx(-1),nynum(0),begnx(0);
    for(k=0;k<203000;k++)
    {
        segyhead_readonetrace_tofmat(head,data);
        nz=data.n_rows;
        if(k==begnx){
            datawrite(data,nz,1,outf);
            if(head.endian=='l'){
                trac0=getendianchange(head.head2.tracl);
                cdp0=getendianchange(head.head2.cdp);
                }
            else{
                trac0=head.head2.tracl;cdp0=head.head2.cdp;
                }
            }

        if(k>begnx){
            datawrite(data,nz,1,outf);
            if(head.endian=='l'){
                trac1=getendianchange(head.head2.tracl);
                cdp1=getendianchange(head.head2.cdp);
                }
            else{
                trac1=head.head2.tracl;cdp1=head.head2.cdp;
                }

            if(keyw==2)
            {cdp1=trac1;cdp0=trac0;}
            if(cdp1==cdp0)
            {
                nx++;
            }
            else
            {
                if(maxnx<nx)
                {maxnx=nx;}
                //cout<<nx<<",";
                numnx(nynum,0)=nx;
                nx=1;
                cdp0=cdp1;
                nynum++;
                if(nynum>=(ny))
                {break;}
                if(head.infile.eof())
                {break;}
            }
        }
    }
    //head.infile.close();
    outf.close();
    //datawrite(numnx,ny,1,"numofnx.bin");
/////////////////////////////
    ifstream inf;
    if(maxnx<1)
    {maxnx=1;}
    data.resize(nz,maxnx);
    inf.open("./data.bin");
    outf.open("./data.oneblock.bin");
    for(i=0;i<ny;i++)
    {
        data.fill(0.0);
        dataread(data,nz,numnx(i,0),inf);
        datawrite(data,nz,maxnx,outf);
    }
    inf.close();
    outf.close();
    head.nz=nz,head.nx=maxnx;
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
