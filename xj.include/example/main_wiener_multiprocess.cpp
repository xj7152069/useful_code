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
fmat dataread_fromsegy_to3dblock_byCDP(segyhead & head, int nline,\
    segyhead2 *suHeadArray2d);

int main(int argc, char * argv[])
{
    int i,j,k,numi,narxmax,ncpu(6);
    float dt(0.002),df,fmax(150.0);
    struct wiener3d par;           //结构体，包含数据和参数
    ofstream outf_globel;
    ifstream inf,inf2;
    ofstream outf,outf2;
    int nz(2901),nx(999),ny(101);    //数据体采样点数
    char filein[399],fileout[399];
    filein[0]='\0';
    strcat(filein,argv[1]);
    fileout[0]='\0';
    strcat(fileout,argv[2]);
    outf_globel.open(fileout);
    strcat(fileout,".err");
    outf.open(fileout);

    segyhead head;
    bool sufile(true);
    head.filename[0]='\0';
    strcat(head.filename,filein);
    segyhead_open(head,sufile);
    //segyhead_endianget(head);
    int nline(91),ntracemax(299*nline),keyw(1),num(0);
    segyhead2 *suHeadArray2d;
    suHeadArray2d=new segyhead2[ntracemax];

while(!head.infile.eof())
{
    cout<<"now is running: "<<num<<endl;
    num++;
    fmat numnx;
    numnx=dataread_fromsegy_to3dblock_byCDP(head, nline,suHeadArray2d);
    nx=head.nx,ny=min(nline,head.begnx),nz=head.nz;
    cout<<nx<<"|"<<ny<<"|"<<nz<<endl;
    wiener3d_parset(nx,ny,nz,par); //结构体初始化函数,并设置一组初始化参数
    //par.nwx=21;par.nwy=21;          //数据窗求解误差逼近泛函
    //par.narx=5;par.nary=5;         //wiener滤波器系数阶数
    par.dig_n=0.001;               //法方程求解对角稳定系数
    par.dig_n2=0.001;
    par.nmovex=par.nwx;            //滑动滤波的滑动步长，
    par.nmovey=par.nwy;            //这里设置为等于数据窗大小，避免多次覆盖滤波
    par.covermax=1;
    df=1.0/nz/dt;
    par.nf1=1,par.nf2=fmax/df;    //滤波频率范围
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
    par.data=par.realrebuildtx-par.data;
    dataline_t.zeros(nx,nz);
    k=0;
    fmat dataonecol(nz,1);
    for(i=0;i<ny;i++)
    {
        dataline_t=par.realrebuildtx.col(i);
        dataline=dataline_t.t();
        for(j=0;j<numnx(i,0);j++){
            outf_globel.write((char*)\
                &suHeadArray2d[k],sizeof(suHeadArray2d[k]));
            dataonecol.col(0)=dataline.col(j);
            k++;
            datawrite(dataonecol,nz,1,outf_globel);
    }}
    k=0;
    for(i=0;i<ny;i++)
    {
        dataline_t=par.data.col(i);
        dataline=dataline_t.t();
        for(j=0;j<numnx(i,0);j++){
            outf.write((char*)\
                &suHeadArray2d[k],sizeof(suHeadArray2d[k]));
            dataonecol.col(0)=dataline.col(j);
            k++;
            datawrite(dataonecol,nz,1,outf);
    }}

}
cout<<head.nx<<"|"<<head.nz<<endl;
    outf_globel.close();
    outf.close();
    return 0;
}

fmat dataread_fromsegy_to3dblock_byCDP(segyhead & head, int nline,\
    segyhead2 *suHeadArray2d)
{
    //keyw==1, use cdp; keyw==2, use trace;
    int i,j,k,nz(3501),ny(15),cdp0,cdp1,trac0,trac1;
    ny=nline;
    fmat data,numnx(ny,1);
    ofstream outf;
    outf.open("./data.bin");
    int nx(0),maxnx(-1),nynum(0),begnx(0);
    for(k=0;k<2000000000;k++)
    {
        segyhead_readonetrace_tofmat(head,data);
        suHeadArray2d[k]=head.head2;
        nz=data.n_rows;
        if(k==begnx){
            datawrite(data,nz,1,outf);
            nx++;
            if(head.endian=='l'){
                cdp0=getendianchange(head.head2.cdp);
                }
            else{
                cdp0=head.head2.cdp;
                }
            }

        if(k>begnx){
            if(head.endian=='l'){
                cdp1=getendianchange(head.head2.cdp);
                }
            else{
                cdp1=head.head2.cdp;
                }

            if(abs(cdp1-cdp0)<0.5)
            {
                datawrite(data,nz,1,outf);
                nx++;
            }
            else
            {
                if(maxnx<nx)
                {maxnx=nx;}
                //cout<<nx<<",";
                numnx(nynum,0)=nx;
                nx=0;
                cdp0=cdp1;
                nynum++;
                if(nynum>=(ny))
                {
                    head.infile.seekg(-nz*1*sizeof(float),ios::cur);
                    head.infile.seekg(-240,ios::cur);
                    //nynum--;
                    break;
                }
                datawrite(data,nz,1,outf);
                nx++;
            }
        }
        float test;
        head.infile.read((char *)(&test), sizeof(test));
        if(head.infile.eof())
        {
            numnx(nynum,0)=nx;
            nynum++;
            if(maxnx<nx)
            {maxnx=nx;}
            break;
        }else{
            head.infile.seekg(-sizeof(test),ios::cur);
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
    for(i=0;i<nynum;i++)
    {
        data.fill(0.0);
        dataread(data,nz,numnx(i,0),inf);
        datawrite(data,nz,maxnx,outf);
    }
    inf.close();
    outf.close();
    head.nz=nz;
    head.nx=maxnx;
    head.begnx=nynum;
    return numnx;
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
