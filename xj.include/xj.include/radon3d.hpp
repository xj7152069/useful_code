/************update in 2021.01.07**************/
/*
    
    
***********************************************/

#ifndef RADON_BEAMFORMING_HPP
#define RADON_BEAMFORMING_HPP

#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;

#include "../xjc.h"

void cal_p_power_3d(struct linerradon3d & par);
void tx_to_fx3d(struct linerradon3d & par);
void tx_to_fx2d(cx_fmat & datafx, fmat & data, int ny, int nf);
void fx_to_tx3d(struct linerradon3d & par);
void fx_to_tx2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz);
void fp_to_tp3d(struct linerradon3d & par);
void fp_to_tp2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz);

/*beamforming/liner_radon变换 传递的参数：
nz（处理数据体Z方向采样点）、nx（处理数据体X方向采样点）、
nf（数据体频率域变换的频率采样）、np（Radon变换倾角采样点）、
par4（RLS_beamforming的迭代次数）、par6（我也忘了这参数干嘛的）、
dz（数据Z方向采样间隔）、dx（数据X方向采样间隔）、
df（数据变换后频率采样间隔）、dp（数据变换后倾角采样间隔）、
p0（数据变换后的中心倾角值）、data（原始数据）、
realdataTP（beamforming/liner_radon变换得到的数据实部）、
realrebuild（beamforming/liner_radon反变换重建的数据实部）、
realdatafft（beamforming/liner_radon变换得到的中间矩阵实部，用于检查算法）、
dataTP（beamforming/liner_radon变换得到的复数数据体）、
rebuild（beamforming/liner_radon反变换重建的复数数据体）、
dig_n（对角加权值）、
par1,par2,par3,par5（这几个par曾用于调试优化结果）、

allAreal[99],allAimag[99]
（是文件路径，用于存放一个转换矩阵A；在处理多批数据时需要重复计算A，
且A较大不好直接存在内存里，故将其存入文件以备读取使用）
*/
struct linerradon3d
{
    int nz,nx,ny,nf,npx,npy,nf1,nf2;
    float dz,dx,dy,df,dpx,dpy,p0x,p0y;
    fcube data,realdataTP,realrebuildtx,realdatafk;
    cx_fcube datafP,dataTP,rebuildfp,rebuildfx,\
        rebuildtx,datafx,datafk;
    fmat p_power;
    char allAreal[99],allAimag[99];

    //cx_fmat Acol1,Arow1;
};

//设置一组常用的初始化参数
void beamforming_parset(int nx, int ny, int nz, int nf,\
    struct linerradon3d & par)
{
    par.nx=nx,par.nz=nz,par.ny=ny,par.nf=nf;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.npx=par.nz,par.dpx=2*par.dz/par.nx/par.dx;
    par.npy=par.nz,par.dpy=2*par.dz/par.ny/par.dy;
    par.nf1=0;par.nf2=par.nf/2;
    par.allAreal[0]='\0';
    strcat(par.allAreal,"./allAreal.bin");
    par.allAimag[0]='\0';
    strcat(par.allAimag,"./allAimag.bin");
}

void beamforming_parupdate(struct linerradon3d & par)
{
    par.df=1.0/par.dz/par.nf;
    par.p0x=-par.dpx*int(par.npx/2);
    par.p0y=-par.dpy*int(par.npy/2);
    par.data.zeros(par.nx,par.ny,par.nz);
    par.datafk.zeros(par.nx,par.ny,par.nf);
    par.datafP.zeros(par.npx,par.npy,par.nf);
    par.datafx.zeros(par.nx,par.ny,par.nf);
    par.dataTP.zeros(par.npx,par.npy,par.nz);
    par.realdataTP.zeros(par.npx,par.npy,par.nz);
    par.realdatafk.zeros(par.nx,par.ny,par.nf);
    par.rebuildfp.zeros(par.npx,par.npy,par.nf);
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nz);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
    par.p_power.zeros(par.npx,par.npy);
}

void beamforming_cleardata(struct linerradon3d & par)
{
    par.datafk.fill(0.0);
    par.datafP.fill(0.0);
    par.datafx.fill(0.0);
    par.dataTP.fill(0.0);
    par.realdataTP.fill(0.0);
    par.realdatafk.fill(0.0);
    par.rebuildfx.fill(0.0);
    par.rebuildtx.fill(0.0);
    par.realrebuildtx.fill(0.0);
}
/*
//预备计算变换需要的A矩阵，便于后续使用
void linerandontransmat(struct linerradon3d & par)
{
    int i,j,k;
    float w,pi(3.1415926),f,dx(par.dx),\
        dp(par.dp),p1(par.p0);
    int nx(par.nx),np(par.np);
    fmat X(nx,1), P(1,np), forA(nx,np);
    cx_fmat forAw(nx,np), A(nx,np);
    ofstream outreal,outimag;
    outreal.open(par.allAreal);
    outimag.open(par.allAimag);
    for(k=0;k<par.nf;k++)  
    {
        w=2.0*pi*par.df*k;  
        for(i=0;i<nx;i++)
        {
            X(i,0)=(i)*dx-(nx*dx)/2.0;
        }
            
        for(i=0;i<np;i++)
        {
            P(0,i)=i*dp+p1;
        }
        forA=X*P;   

        for(i=0;i<nx;i++)
        {
            for(j=0;j<np;j++)
            {
                (forAw(i,j).real(0.0));
                (forAw(i,j).imag(0.0));
                (forAw(i,j).imag(forA(i,j)*w));
                //forAw(i,j)=forAw(i,j);
            }
        }
        A=exp(forAw);
        datawrite(forA=real(A),nx,np,outreal);
        datawrite(forA=imag(A),nx,np,outimag);
    }
    outreal.close(),outimag.close();
}
*/
//线性Radon变换
void linerradon(struct linerradon3d & par)
{
    int kf,kpx,kpy,kx,ky;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dy(par.dy),\
        dpx(par.dpx),dpy(par.dpy),p0x(par.p0x),\
        p0y(par.p0y),dz(par.dz);
    int nx(par.nx),npx(par.npx),nf(par.nf),\
        ny(par.ny),npy(par.npy);

    cx_fmat forA(nx,ny),A(nx,ny),a(1,1),B(nx,ny);
    forA.fill(0.0);
    ofstream outf;
    ifstream infreal,infimag;

    tx_to_fx3d(par);
    for(kf=par.nf1;kf<par.nf2;kf++)  
    {
        w=2.0*pi*par.df*kf; 
        cout<<"now is running kf = "<<kf<<endl;
        B=par.datafx.slice(kf);
        for(kpx=0;kpx<npx;kpx++)
        {
            fpx=kpx*dpx+p0x;
            for(kpy=0;kpy<npy;kpy++)
            {
                fpy=kpy*dpy+p0y;
                for(kx=0;kx<nx;kx++)
                {
                    fx=kx*dx-(nx*dx)/2.0;
                    for(ky=0;ky<ny;ky++)
                    {
                        fy=ky*dy-(ny*dy)/2.0;
                        forA(kx,ky).imag(w*(fx*fpx+fy*fpy));
                    }
                }
                A=exp(forA);
                a=cx_fmatmul(A,B);
                par.datafP(kpx,kpy,kf)=a(0,0);
            }
        }
        
    }
    fp_to_tp3d(par);
    par.realdataTP=real(par.dataTP);
}

//重建数据，重建之前可以对数据按倾角去噪等
void rebuildsignal(struct linerradon3d & par)
{
    int kf,kpx,kpy,kx,ky;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dy(par.dy),\
        dpx(par.dpx),dpy(par.dpy),p0x(par.p0x),\
        p0y(par.p0y),dz(par.dz);
    int nx(par.nx),npx(par.npx),nf(par.nf),\
        ny(par.ny),npy(par.npy);

    cx_fmat forA(npx,npy),A(npx,npy),a(1,1),B(npx,npy);
    forA.fill(0.0);
    ofstream outf;
    ifstream infreal,infimag;

    for(kf=par.nf1;kf<par.nf2;kf++)  
    {
        w=2.0*pi*par.df*kf; 
        cout<<"now is running kf = "<<kf<<endl;
        B=par.datafP.slice(kf);
        for(kx=0;kx<nx;kx++)
        {
            fx=kx*dx-(nx*dx)/2.0;
            for(ky=0;ky<ny;ky++)
            {
                fy=ky*dy-(ny*dy)/2.0;
                for(kpx=0;kpx<npx;kpx++)
                {
                    fpx=kpx*dpx+p0x;
                    for(kpy=0;kpy<npy;kpy++)
                    {
                        fpy=kpy*dpy+p0y;
                        forA(kpx,kpy).imag(w*(fx*fpx+fy*fpy));
                    }
                }
                A=exp(forA);
                a=cx_fmatmul(A,B);
                par.rebuildfx(kx,ky,kf)=a(0,0);
            }
        }
        
    }
    fx_to_tx3d(par);
    par.realrebuildtx=real(par.rebuildtx);
}

/////////////
void cal_p_power_3d(struct linerradon3d & par)
{
    int i,j,k;
    par.p_power.fill(0.0);

    for(i=0;i<par.npx;i++)
    {
        for(j=0;j<par.npy;j++)
        {
            for(k=par.nf1;k<par.nf2;k++)
            {
                par.p_power(i,j)=par.p_power(i,j)+\
                    abs(par.datafP(i,j,k));
            }
        }
    }
}

void tx_to_fx3d(struct linerradon3d & par)
{
    int i;
    fmat data(par.ny,par.nz);
    cx_fmat datafx(par.ny,par.nf);

    for(i=0;i<par.nx;i++)
    {
        tx_to_fx2d(datafx, data=par.data.row(i),\
            par.ny, par.nf);
        par.datafx.row(i)=datafx;
    } 
}

void tx_to_fx2d(cx_fmat & datafx, fmat & data, int ny, int nf)
{
    int i;
    for(i=0;i<ny;i++)
    {
        datafx.row(i)=fft(data.row(i),nf);
    } 
}

void fx_to_tx3d(struct linerradon3d & par)
{
    int i;
    cx_fmat data(par.ny,par.nz);
    cx_fmat datafx(par.ny,par.nf);

    for(i=0;i<par.nx;i++)
    {
        fx_to_tx2d(data,datafx=par.rebuildfx.row(i),\
            par.ny, par.nz);
        par.rebuildtx.row(i)=data;
    } 
}

void fx_to_tx2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz)
{
    int i;
    for(i=0;i<ny;i++)
    {
        data.row(i)=ifft(datafx.row(i),nz);
    } 
}

void fp_to_tp3d(struct linerradon3d & par)
{
    int i;
    cx_fmat data(par.npy,par.nz);
    cx_fmat datafp(par.npy,par.nf);

    for(i=0;i<par.npx;i++)
    {
        fp_to_tp2d(data,datafp=par.datafP.row(i),\
            par.npy, par.nz);
        par.dataTP.row(i)=data;
    } 
}

void fp_to_tp2d(cx_fmat & data, cx_fmat & datafx, int npy, int nz)
{
    int i;
    for(i=0;i<npy;i++)
    {
        data.row(i)=ifft(datafx.row(i),nz);
    } 
}

#endif
