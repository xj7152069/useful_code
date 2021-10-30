/************update in 2021.01.07**************/
/*
    
    
***********************************************/

#ifndef RADON3D_HPP
#define RADON3D_HPP
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <thread>
#include <future>
#include "../xjc.h"
using namespace std;
using namespace arma;

void cal_p_power_3d(struct linerradon3d & par);
void tx_to_fx3d_linerradon3d(struct linerradon3d & par);
void tx_to_fx2d(cx_fmat & datafx, fmat & data, int ny, int nf);
void tp_to_fp3d_linerradon3d(struct linerradon3d & par);
void fx_to_tx3d_linerradon3d(struct linerradon3d & par);
void fx_to_tx2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz);
void fp_to_tp3d_linerradon3d(struct linerradon3d & par);
void fp_to_tp2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz);
void linerradon(struct linerradon3d & par);

inline cx_fmat cx_fmatmul(cx_fmat & mat1, cx_fmat & mat2);
inline cx_fmat cx_hessianelement1d(float w,\
 int nx, float kx1, float dx, float dpx);
cx_fmat cx_hessianelement2d(float w,\
 int nx, float kx1, float dx, float dpx,\
 int ny, float ky1, float dy, float dpy);
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
    int nz,nx,ny,nf,npx,npy,nf1,nf2,numthread,rulef1,rulef2;
    float dz,dx,dy,df,dpx,dpy,p0x,p0y,px_center,py_center,dig_n;
    fcube data,realdataTP,realrebuildtx;
    cx_fcube datafx,rebuildfx,datafP,rebuildfp,dataTP,rebuildtx;
    fmat p_power,py_coord,px_coord,x_coord,y_coord;
    char allAreal[99],allAimag[99];
    //cx_fmat Acol1,Arow1;
};
//设置一组常用的初始化参数
void beamforming_parset(int nx, int ny, int nz, int nf,\
    struct linerradon3d & par)
{
    par.rulef1=5,par.rulef2=30;
    par.nx=nx,par.nz=nz,par.ny=ny,par.nf=nf;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.npx=par.nz,par.dpx=2*par.dz/par.nx/par.dx;
    par.npy=par.nz,par.dpy=2*par.dz/par.ny/par.dy;
    par.px_center=0;par.py_center=0;
    par.p0x=-par.dpx*int(par.npx/2)+par.px_center;
    par.p0y=-par.dpy*int(par.npy/2)+par.py_center;
    par.nf1=0;par.nf2=par.nf/2;
    par.numthread=1;
    par.dig_n=(0.001);
    par.allAreal[0]='\0';
    strcat(par.allAreal,"./allAreal.bin");
    par.allAimag[0]='\0';
    strcat(par.allAimag,"./allAimag.bin");
} 
void beamforming_parset(int nx, int ny, int nz,\
    struct linerradon3d & par)
{
    struct linerradon3d * ppar;
    ppar=&par;
    int nf(1);
    while (nf<nz)
    {
        nf*=2;
    }
    cout<<"nf = "<<nf<<endl;
    beamforming_parset(nx,ny,nz,nf,ppar[0]);
}

void beamforming_parupdate(struct linerradon3d & par)
{ 
    par.df=1.0/par.dz/par.nf;
    par.p0x=-par.dpx*int(par.npx/2)+par.px_center;
    par.p0y=-par.dpy*int(par.npy/2)+par.py_center;
    par.data.zeros(par.nx,par.ny,par.nz);
    par.datafP.zeros(par.npx,par.npy,par.nf);
    par.datafx.zeros(par.nx,par.ny,par.nf);
    par.dataTP.zeros(par.npx,par.npy,par.nz);
    par.realdataTP.zeros(par.npx,par.npy,par.nz);
    par.rebuildfp.zeros(par.npx,par.npy,par.nf);
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nz);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
    par.p_power.zeros(par.npx,par.npy);
    par.dig_n=par.dig_n*par.nx*par.ny;
    par.px_coord.zeros(par.npx,1);
    par.py_coord.zeros(par.npy,1);
    par.x_coord.zeros(par.nx,1);
    par.y_coord.zeros(par.ny,1);
    int k;
    for(k=0;k<par.npx;k++)
    {
        par.px_coord(k,0)=k*par.dpx+par.p0x;
    }
    for(k=0;k<par.npy;k++)
    {
        par.py_coord(k,0)=k*par.dpy+par.p0y;
    }
    for(k=0;k<par.nx;k++)
    {
        par.x_coord(k,0)=k*par.dx-(par.nx*par.dx)/2.0;
    }
    for(k=0;k<par.ny;k++)
    {
        par.y_coord(k,0)=k*par.dy-(par.ny*par.dy)/2.0;
    }
    par.rulef1=round(par.rulef1/par.df);
    par.rulef2=round(par.rulef2/par.df);
    par.nf1=round(par.nf1/par.df);
    par.nf2=round(par.nf2/par.df);
    if(par.nf2>par.nf/2)
    {
        par.nf2=par.nf/2;
    }
    if(par.rulef2>par.nf2/2)
    {
        par.rulef2=par.nf2/2;
    }
}

void beamforming_cleardata(struct linerradon3d & par)
{
    par.datafP.fill(0.0);
    par.datafx.fill(0.0);
    par.dataTP.fill(0.0);
    par.realdataTP.fill(0.0);
    par.rebuildfx.fill(0.0);
    par.rebuildtx.fill(0.0);
    par.realrebuildtx.fill(0.0);
}
/////////////////////////beamforminginv3d///////////////////////////
struct beamforminginv3d
{
    fmat digw_fmat_npxnpy, **base_fmatp2npxnpy_nxny;
    //cx_fmat *hessianinv_cxfmat_p1nf_npxynpxy;
    cx_fmat *hessianinv_cxfmat_p1ncpu_npxynpxy;
};

void beamforminginv3d_parset(struct linerradon3d & par,\
 struct beamforminginv3d & parinv)
{
    parinv.digw_fmat_npxnpy.zeros(par.npx,par.npy);
    int ncpu(par.numthread);
    parinv.hessianinv_cxfmat_p1ncpu_npxynpxy=new cx_fmat[ncpu];
    int i,j,k,n;
    n=par.npx*par.npy;
    for(k=0;k<ncpu;k++)
    {
        parinv.hessianinv_cxfmat_p1ncpu_npxynpxy[k].zeros(n,n);
    }
}
void beamforminginv3d_getdigfmat(struct linerradon3d & par,\
 struct beamforminginv3d & parinv)
{
    parinv.digw_fmat_npxnpy.fill(0.0);
    int k;
    float maxpower;
    for(k=par.rulef1;k<par.rulef2;k++)
    {
        parinv.digw_fmat_npxnpy+=abs(par.datafP.slice(k));
    }
    maxpower=max(max(parinv.digw_fmat_npxnpy));
    parinv.digw_fmat_npxnpy=(0.001/par.nx/par.ny\
        +parinv.digw_fmat_npxnpy/maxpower);
    parinv.digw_fmat_npxnpy=1.0/parinv.digw_fmat_npxnpy;
    maxpower=max(max(parinv.digw_fmat_npxnpy));
    parinv.digw_fmat_npxnpy=parinv.digw_fmat_npxnpy/maxpower;
    parinv.digw_fmat_npxnpy=parinv.digw_fmat_npxnpy*par.dig_n;
    par.p_power=parinv.digw_fmat_npxnpy;
}

void beamforminginv3d_hessianget_thread(struct linerradon3d * par,\
 struct beamforminginv3d * par2, int pncpu, int pnf)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926535897932384626);
    float df(par[0].df),dx(par[0].dx),dy(par[0].dy),\
        dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x),\
        p0y(par[0].p0y),dz(par[0].dz),dig_n(par[0].dig_n);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy);

    float ky1(-(ny*dy)/2.0);
    float kx1(-(nx*dx)/2.0);
    cx_fmat a(1,1);

    kf=pnf;
    if(kf<par[0].nf2) 
    {
        w=2.0*pi*df*(kf);  //!!!take care of when (par.nf1 != 0)
        for(kpx=0;kpx<npx;kpx++)
        {
        for(kpy=0;kpy<npy;kpy++)
        {
            i=kpx*npy+kpy;
            float fpx1(kpx*dpx+p0x);
            float fpy1(kpy*dpy+p0y);
            for(kpx2=0;kpx2<npx;kpx2++)
            {
            for(kpy2=0;kpy2<npy;kpy2++)
            {
                j=kpx2*npy+kpy2;
                float fpx2(kpx2*dpx+p0x);
                float fpy2(kpy2*dpy+p0y);
                a=cx_hessianelement2d(w,nx,kx1,dx,fpx1-fpx2,\
                    ny,ky1,dy,fpy1-fpy2);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,j)=a(0,0);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](j,i)=a(0,0); 
            }
            }        
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,i)+=dig_n; //hessian dig
            //par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,i)+=par->p_power(kpx,kpy);
            kpx2=npx;kpy2=npy;
        }
        } 
    }
}

void beamforminginv3d_beamget_thread(struct linerradon3d * par,\
 struct beamforminginv3d * par2, int pncpu, int pnf)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float df(par[0].df),dx(par[0].dx),dy(par[0].dy),\
        dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x),\
        p0y(par[0].p0y),dz(par[0].dz);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy);

    float ky1(-(ny*dy)/2.0);
    float kx1(-(nx*dx)/2.0);

    kf=pnf;
    if(kf<par[0].nf2) 
    {
        cx_fmat datafp(npx,npy),a(1,1);
        for(kpx=0;kpx<npx;kpx++)
        {
        for(kpy=0;kpy<npy;kpy++)
        {
            i=kpx*npy+kpy;
            a(0,0).real(0.0);
            a(0,0).imag(0.0);
            for(kpx2=0;kpx2<npx;kpx2++)
            {
            for(kpy2=0;kpy2<npy;kpy2++)
            {
                j=kpx2*npy+kpy2;
                a(0,0)+=par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,j)\
                    *par[0].datafP(kpx2,kpy2,kf);
            }
            }
            datafp(kpx,kpy)=a(0,0);
        }
        }
        par[0].datafP.slice(kf)=datafp;
    }
}

void beamforminginv3d(struct linerradon3d & par)
{
    int numf(par.nf2-par.nf1),ncpu(par.numthread);
    int kf,kcpu,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dy(par.dy),\
        dpx(par.dpx),dpy(par.dpy),p0x(par.p0x),\
        p0y(par.p0y),dz(par.dz);
    int nx(par.nx),npx(par.npx),nf(par.nf),\
        ny(par.ny),npy(par.npy);

    struct beamforminginv3d par2;
    struct linerradon3d * ppar;
    ppar=&par;
    linerradon(ppar[0]);
    beamforminginv3d_parset(ppar[0], par2);
    beamforminginv3d_getdigfmat(ppar[0], par2);

    cout<<"now is running hessian_inv : "<<endl;
    //beamforminginv3d_hessianinv(ppar[0], par2);

thread *pcal;
pcal=new thread[ncpu];
for(kf=0;kf<numf;kf+=ncpu)  
{
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu]=thread(beamforminginv3d_hessianget_thread,\
            &par,&par2,kcpu,kf+kcpu+par.nf1);
    }
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu].join();
    }
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
    if((kf+kcpu+par.nf1)<par.nf2) 
    {
        cout<<"hessian_inv kf = "<<kf+kcpu+par.nf1<<endl;
        par2.hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu]=\
            inv(par2.hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu]);
    }
    }
//////////////////////////////////
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu]=thread(beamforminginv3d_beamget_thread,\
            &par,&par2,kcpu,kf+kcpu+par.nf1);
    }
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu].join();
    }
}
    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
}

void beamforminginv3d_thread(struct linerradon3d & par)
{
    int numf(par.nf2-par.nf1),ncpu(par.numthread);
    int kf,kcpu,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dy(par.dy),\
        dpx(par.dpx),dpy(par.dpy),p0x(par.p0x),\
        p0y(par.p0y),dz(par.dz);
    int nx(par.nx),npx(par.npx),nf(par.nf),\
        ny(par.ny),npy(par.npy);

    struct beamforminginv3d par2;
    struct linerradon3d * ppar;
    ppar=&par;
    linerradon(ppar[0]);
    beamforminginv3d_parset(ppar[0], par2);
    beamforminginv3d_getdigfmat(ppar[0], par2);

    cout<<"now is running hessian_inv : "<<endl;
    //beamforminginv3d_hessianinv(ppar[0], par2);

thread *pcal;
pcal=new thread[ncpu];
for(kf=0;kf<numf;kf+=ncpu)  
{
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu]=thread(beamforminginv3d_hessianget_thread,\
            &par,&par2,kcpu,kf+kcpu+par.nf1);
    }
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu].join();
    }
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
    if((kf+kcpu+par.nf1)<par.nf2) 
    {
        cout<<"hessian_inv kf = "<<kf+kcpu+par.nf1<<endl;
        par2.hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu]=\
            inv(par2.hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu]);
    }
    }
//////////////////////////////////
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu]=thread(beamforminginv3d_beamget_thread,\
            &par,&par2,kcpu,kf+kcpu+par.nf1);
    }
    for(kcpu=0;kcpu<ncpu;kcpu++)  
    {
        pcal[kcpu].join();
    }
}
    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
}

////////////////////linerradon 3D transform/////////////////////////
void linerradon_getA(cx_fmat & A, cx_fmat & basex,\
    cx_fmat & basey, int kpx, int kpy, int nx, int ny)
{
    int kx,ky;
    for(kx=0;kx<nx;kx++)
    {for(ky=0;ky<ny;ky++)
    {
        A(kx,ky)=basex(kx,kpx)*basey(ky,kpy);
    }}
}
//线性Radon变换，无反演
void linerradon_fthread(struct linerradon3d * par,int pnf1, int pnf2)
{
    int kf,kpx,kpy,kx,ky;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par->df),dx(par->dx),dy(par->dy),\
        dpx(par->dpx),dpy(par->dpy),p0x(par->p0x),\
        p0y(par->p0y),dz(par->dz);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy);

    cx_fmat forA(nx,ny),A(nx,ny),a(1,1),B(nx,ny);
    cx_fmat basey(ny,npy,fill::zeros),basex(nx,npx,fill::zeros);

    forA.fill(0.0);
    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        //cout<<"now is running kf = "<<kf<<endl;
        basey.zeros(ny,npy);basex.zeros(nx,npx);
        basey.set_imag(w*(par->y_coord*par->py_coord.t()));  
        basex.set_imag(w*(par->x_coord*par->px_coord.t()));  
        basey=exp(basey);basex=exp(basex);
        B=par->datafx.slice(kf);
        for(kpx=0;kpx<npx;kpx++)
        {
            for(kpy=0;kpy<npy;kpy++)
            {
                A=basex.col(kpx)*basey.col(kpy).st();
                //linerradon_getA(A,basex,basey,kpx,kpy,nx,ny);
                a=cx_fmatmul(A,B);
                par->datafP(kpx,kpy,kf)=a(0,0);
            }
        }
    }
}

void linerradon(struct linerradon3d & par)
{
    tx_to_fx3d_linerradon3d(par);
    int pn(par.numthread),pnf1,pnf2,k;
    float dnf;
    thread *pcal;
    pcal=new thread[pn];
    dnf=float(par.nf2-par.nf1)/pn;

    for(k=0;k<pn-1;k++)
    {
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(linerradon_fthread,&par,pnf1,pnf2);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(linerradon_fthread,&par,pnf1,pnf2);

    for(k=0;k<pn;k++)
    {
        pcal[k].join();
    }

    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP)/(par.npx*par.npy);
    delete [] pcal;
}

////////////////////rebuild 3D data/////////////////////////
//重建数据，重建之前可以对数据按倾角去噪等
void rebuildsignal_getA(cx_fmat & A, cx_fmat & basex,\
    cx_fmat & basey, int kx, int ky, int npx, int npy)
{
    int kpx,kpy;
    for(kpx=0;kpx<npx;kpx++)
    {for(kpy=0;kpy<npy;kpy++)
    {
        A(kpx,kpy)=basex(kx,kpx)*basey(ky,kpy);
    }}
}
void rebuildsignal_fthread(struct linerradon3d * par,int pnf1, int pnf2)
{
    int kf,kpx,kpy,kx,ky;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par->df),dx(par->dx),dy(par->dy),\
        dpx(par->dpx),dpy(par->dpy),p0x(par->p0x),\
        p0y(par->p0y),dz(par->dz);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy);

    cx_fmat forA(npx,npy),A(npx,npy),a(1,1),B(npx,npy);
    cx_fmat basey(ny,npy,fill::zeros),basex(nx,npx,fill::zeros);

    forA.fill(0.0);
    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        basey.zeros(ny,npy);basex.zeros(nx,npx);
        basey.set_imag(w*(par->y_coord*par->py_coord.t()));  
        basex.set_imag(w*(par->x_coord*par->px_coord.t()));  
        basey=exp(basey);basex=exp(basex);
        //cout<<"now is running kf = "<<kf<<endl;
        B=par->datafP.slice(kf);
        for(kx=0;kx<nx;kx++)
        {
            for(ky=0;ky<ny;ky++)
            {
                A=basex.row(kx).st()*basey.row(ky);
                //rebuildsignal_getA(A,basex,basey,kx,ky,npx,npy);
                A.set_imag(-imag(A));
                a=cx_fmatmul(A,B);
                par->rebuildfx(kx,ky,kf).real(4*real(a(0,0)));
                par->rebuildfx(kx,ky,kf).imag(4*imag(a(0,0)));
            }
        }
    }

}

void rebuildsignal(struct linerradon3d & par)
{
    tp_to_fp3d_linerradon3d(par);

    int pn(par.numthread),pnf1,pnf2,k;
    float dnf;
    thread *pcal;
    pcal=new thread[pn];
    dnf=float(par.nf2-par.nf1)/pn;

    for(k=0;k<pn-1;k++)
    {
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(rebuildsignal_fthread,&par,pnf1,pnf2);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(rebuildsignal_fthread,&par,pnf1,pnf2);

    for(k=0;k<pn;k++)
    {
        pcal[k].join();
    }
   
    fx_to_tx3d_linerradon3d(par);
    par.realrebuildtx=real(par.rebuildtx);
    delete [] pcal;
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

void tp_to_fp3d_linerradon3d(struct linerradon3d & par)
{
    int i;
    fmat data(par.npy,par.nz);
    cx_fmat datafp(par.npy,par.nf);

    for(i=0;i<par.npx;i++)
    {
        tx_to_fx2d(datafp, data=par.realdataTP.row(i),\
            par.npy, par.nf);
        par.datafP.row(i)=datafp;
    } 
}

void tx_to_fx3d_linerradon3d(struct linerradon3d & par)
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

void fx_to_tx3d_linerradon3d(struct linerradon3d & par)
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

void fp_to_tp3d_linerradon3d(struct linerradon3d & par)
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

inline cx_fmat cx_fmatmul(cx_fmat & mat1, cx_fmat & mat2)
{
    int nz,nx;
    nz=mat1.n_rows;
    nx=mat1.n_cols;

    cx_fmat a(1,1,fill::zeros);
    int i,j;
    for(i=0;i<nz;i++)
    {
        a+=mat1.row(i)*mat2.row(i).st();
    }
    /*
    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
        {
            a(0,0)+=mat1(i,j)*mat2(i,j);
        }
    }
    */
    return a;
}

cx_fmat cx_hessianelement2d(float w,\
 int nx, float kx1, float dx, float dpx,\
 int ny, float ky1, float dy, float dpy)
{
    cx_fmat a(1,1,fill::zeros),b(1,1,fill::zeros);
    //float n(1/((nx+1)*(ny+1)));
    float n2(0.0000000000001);
    float pdx(abs(cos(w*dpx*dx)-1.0)),pdy(abs(cos(w*dpy*dy)-1.0));

    if(pdx>n2 && pdy>n2)
    {
        a=cx_hessianelement1d(w, nx, kx1, dx, dpx);
        b=cx_hessianelement1d(w, ny, ky1, dy, dpy);
    }
    else if(pdx<=n2 && pdy>n2)
    {
        a=nx;
        b=cx_hessianelement1d(w, ny, ky1, dy, dpy);
    }
    else if(pdx>n2 && pdy<=n2)
    {
        a=cx_hessianelement1d(w, nx, kx1, dx, dpx);
        b=ny;
    }
    else if(pdx<=n2 && pdy<=n2)
    {
        a=nx;
        b=ny;
    }
    
    a=a*b;
    return a;
}

inline cx_fmat cx_hessianelement1d(float w,\
 int nx, float kx1, float dx, float dpx)
{
    cx_fmat a(1,1,fill::zeros);
    float a1,a2,b1,b2,c1,c2,d1,d2;

    a1=cos(w*dpx*kx1);
    a2=sin(w*dpx*kx1);
    
    b1=cos(w*dpx*nx*dx)-1.0;
    b2=sin(w*dpx*nx*dx);
    
    c1=cos(w*dpx*dx)-1.0;
    c2=sin(w*dpx*dx);

    d1=a1*b1-a2*b2;
    d2=a2*b1+a1*b2;
    a1=(d1*c1+d2*c2)/(c1*c1+c2*c2);
    a2=(d2*c1-d1*c2)/(c1*c1+c2*c2);
    
    a(0,0).real(a1);
    a(0,0).imag(a2);  //take care of the error!!!
    return a;
}

/////////////////////////beamformingCG3d///////////////////////////

struct beamformingCG3d
{
    int nx,nz,npx,npz,numi;
    float *a;
    fmat *fmatdigwnpxy;
    cx_fmat *b_cxfmatp1ncpunpxy;
    cx_fmat **r_k1_cxfmatp1ncpunpxy;
    cx_fmat **r_k2_cxfmatp1ncpunpxy;
    cx_fmat **p_k1_cxfmatp1ncpunpxy;
    cx_fmat **p_k2_cxfmatp1ncpunpxy;
    cx_fmat **x_k1_cxfmatp1ncpunpxy;
    cx_fmat **x_k2_cxfmatp1ncpunpxy;
    struct linerradon3d * radon3d_par;
    struct beamforminginv3d * inv3d_par;
};

void beamformingCG3d_parset(struct beamformingCG3d & cg,
    struct linerradon3d & par, struct beamforminginv3d & par2)
{
    cg.radon3d_par=&par;
    cg.inv3d_par=&par2;
    int i,j,k,ncpu(par.numthread);
    int kf,kpx,kpy,kx,ky;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dy(par.dy),\
        dpx(par.dpx),dpy(par.dpy),p0x(par.p0x),\
        p0y(par.p0y),dz(par.dz);
    int nx(par.nx),npx(par.npx),nf(par.nf),\
        ny(par.ny),npy(par.npy);
    cg.nz=(par.ny),cg.nx=(par.nx),cg.npx=(par.npx),cg.npz=(par.npy);
    cg.numi=5;cg.a=new float[ncpu];

    cg.fmatdigwnpxy=new fmat[1];
    cg.fmatdigwnpxy[0].zeros(npx,npy);

    cg.b_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.r_k1_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.r_k2_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.p_k1_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.p_k2_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.x_k1_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.x_k2_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    for(i=0;i<ncpu;i++)
    {
        cg.r_k1_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.r_k2_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.p_k1_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.p_k2_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.x_k1_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.x_k2_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.b_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.r_k1_cxfmatp1ncpunpxy[i][0].zeros(npx,npy);
        cg.r_k2_cxfmatp1ncpunpxy[i][0].zeros(npx,npy);
        cg.p_k1_cxfmatp1ncpunpxy[i][0].zeros(npx,npy);
        cg.p_k2_cxfmatp1ncpunpxy[i][0].zeros(npx,npy);
        cg.x_k1_cxfmatp1ncpunpxy[i][0].zeros(npx,npy);
        cg.x_k2_cxfmatp1ncpunpxy[i][0].zeros(npx,npy);
    }
}

void cxfmatget_Ap(cx_fmat & Ap,cx_fmat & A,cx_fmat & p)
{
    int npx(Ap.n_rows),npy(Ap.n_cols);
    int i,j,kpx,kpx2,kpy,kpy2;
    cx_fmat a(1,1),b(1,1);

    for(kpx=0;kpx<npx;kpx++)
    {
    for(kpy=0;kpy<npy;kpy++)
    {
        i=kpx*npy+kpy;
        a(0,0).real(0.0);
        a(0,0).imag(0.0);
        for(kpx2=0;kpx2<npx;kpx2++)
        {
        for(kpy2=0;kpy2<npy;kpy2++)
        {
            j=kpx2*npy+kpy2;
            a(0,0)+=A(i,j)*p(kpx2,kpy2);
        }
        }
        Ap(kpx,kpy)=a(0,0);
    }
    }
}

int beamformingCG3d_once(struct beamformingCG3d & cg, int kcpu, float aaa)
{
    //struct Axdata Ax;
    float bc; 
    //Ax.nz=cg.nz,Ax.nx=cg.nx,Ax.npx=cg.npx,Ax.npz=cg.npz;
    cx_fmat me(1,1),ak(1,1),bk(1,1),me_cxfmatnpxy;
    me_cxfmatnpxy.copy_size(cg.b_cxfmatp1ncpunpxy[kcpu]); 
    //Ax.outcxfmatp1npxy=&me_cxfmatnpxy; //!!!!!
    //Ax.fmatp1digwnpxy=cg.fmatdigwnpxy;
    //Ax.cxfmatp2npxynxy=cg.RH_cxfmatp3ncpunpxynxy[kcpu];
    cx_fmat Ap(cg.npx,cg.npz);
    cx_float a,an;
    an.real(-1.0),an.imag(-1.0);
    //float a;

//cal_ak || bc
    me_cxfmatnpxy.set_real(real(cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]));
    me_cxfmatnpxy.set_imag(-imag(cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]));
    ak=cx_fmatmul(cg.r_k1_cxfmatp1ncpunpxy[kcpu][0],me_cxfmatnpxy);
    cxfmatget_Ap(Ap,cg.inv3d_par->hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu],\
            cg.p_k1_cxfmatp1ncpunpxy[kcpu][0]);
    me_cxfmatnpxy.set_real(real(Ap));
    me_cxfmatnpxy.set_imag(-imag(Ap));
    me=cx_fmatmul(cg.p_k1_cxfmatp1ncpunpxy[kcpu][0],me_cxfmatnpxy);
    bc=abs(me(0,0));//cout<<bc<<"|";
    if(bc==0)
    {return 1;}
    cg.a[kcpu]=bc;
    a=((ak(0,0))/(me(0,0)));
    //a=real((ak(0,0))/(me(0,0)));
    //if(bc>aaa)
    if(false)
    {
    cg.x_k2_cxfmatp1ncpunpxy[kcpu][0]=cg.x_k1_cxfmatp1ncpunpxy[kcpu][0]\
        +a*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
    cxfmatget_Ap(Ap,cg.inv3d_par->hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu],\
        cg.x_k2_cxfmatp1ncpunpxy[kcpu][0]);
    cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]=cg.b_cxfmatp1ncpunpxy[kcpu]-Ap;
    cg.p_k2_cxfmatp1ncpunpxy[kcpu][0]=cg.r_k2_cxfmatp1ncpunpxy[kcpu][0];
    }
    else
    { 
//cal_xk2
    cg.x_k2_cxfmatp1ncpunpxy[kcpu][0]=cg.x_k1_cxfmatp1ncpunpxy[kcpu][0]\
        +a*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
//cal_rk2
    cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]=cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]\
        -a*Ap;
    //if(sum(sum(abs(cg.r_k2_cxfmatp1ncpunpxy[kcpu])))==0)
    //{return 1;}
//cal_bk || bc
    me_cxfmatnpxy.set_real(real(cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]));
    me_cxfmatnpxy.set_imag(-imag(cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]));
    bk=cx_fmatmul(cg.r_k2_cxfmatp1ncpunpxy[kcpu][0],me_cxfmatnpxy);
    bc=real(bk(0,0))/real(ak(0,0));
//cal_pk2
    cg.p_k2_cxfmatp1ncpunpxy[kcpu][0]=cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]\
        +bc*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
    }

//swap:
{
    cx_fmat *swap;
    swap=cg.p_k1_cxfmatp1ncpunpxy[kcpu];
    cg.p_k1_cxfmatp1ncpunpxy[kcpu]=cg.p_k2_cxfmatp1ncpunpxy[kcpu];
    cg.p_k2_cxfmatp1ncpunpxy[kcpu]=swap;
    swap=cg.r_k1_cxfmatp1ncpunpxy[kcpu];
    cg.r_k1_cxfmatp1ncpunpxy[kcpu]=cg.r_k2_cxfmatp1ncpunpxy[kcpu];
    cg.r_k2_cxfmatp1ncpunpxy[kcpu]=swap;
    swap=cg.x_k1_cxfmatp1ncpunpxy[kcpu];
    cg.x_k1_cxfmatp1ncpunpxy[kcpu]=cg.x_k2_cxfmatp1ncpunpxy[kcpu];
    cg.x_k2_cxfmatp1ncpunpxy[kcpu]=swap;
}
    return 0;
}

void beamformingCG3d_fthread(struct beamformingCG3d * cg, \
    struct linerradon3d * par, int pnf1, int pnf2, int kcpu)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j,k;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy,dig_n(par->dig_n);
    float w,pi(3.1415926),aaa(-1);
    float df(par->df),dx(par->dx),dy(par->dy),\
        dpx(par->dpx),dpy(par->dpy),p0x(par->p0x),\
        p0y(par->p0y),dz(par->dz);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy),cgstop;
    float ky1(-(ny*dy)/2.0);
    float kx1(-(nx*dx)/2.0);

    cx_fmat a(1,1),Ap(npx,npy);
    struct beamforminginv3d * par2;
    par2=cg->inv3d_par;

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        cout<<"now is running kf = "<<kf<<endl;
        cg->b_cxfmatp1ncpunpxy[kcpu]=par->datafP.slice(kf); //b_cxfmat 
        for(kpx=0;kpx<npx;kpx++)
        {
        for(kpy=0;kpy<npy;kpy++)
        {
            i=kpx*npy+kpy;
            float fpx1(kpx*dpx+p0x);
            float fpy1(kpy*dpy+p0y);
            for(kpx2=0;kpx2<npx;kpx2++)
            {
            for(kpy2=0;kpy2<npy;kpy2++)
            {
                j=kpx2*npy+kpy2;
                float fpx2(kpx2*dpx+p0x);
                float fpy2(kpy2*dpy+p0y);
                a=cx_hessianelement2d(w,nx,kx1,dx,fpx1-fpx2,\
                    ny,ky1,dy,fpy1-fpy2);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu](i,j)=a(0,0);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu](j,i)=a(0,0); 
            }
            }        
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu](i,i)+=dig_n; //hessian dig
            //par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu](i,i)+=par->p_power(kpx,kpy); //hessian dig
        }
        }
        //cg->x_k1_cxfmatp1ncpunpxy[kcpu][0]=cg->b_cxfmatp1ncpunpxy[kcpu];
        cg->x_k1_cxfmatp1ncpunpxy[kcpu][0].zeros(npx,npy);;
        cxfmatget_Ap(Ap,par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu],\
            cg->x_k1_cxfmatp1ncpunpxy[kcpu][0]);
        cg->r_k1_cxfmatp1ncpunpxy[kcpu][0]=cg->b_cxfmatp1ncpunpxy[kcpu]-Ap;
        cg->p_k1_cxfmatp1ncpunpxy[kcpu][0]=cg->r_k1_cxfmatp1ncpunpxy[kcpu][0];
        cgstop=beamformingCG3d_once(cg[0],kcpu,aaa);aaa=cg->a[kcpu];
        for(i=0;i<cg->numi;i++)
        {
            cgstop=beamformingCG3d_once(cg[0],kcpu,aaa);aaa=cg->a[kcpu];
            if(cgstop==1)
            {break;} 
        }
        aaa=-1;
        cout<<"CG times:"<<i<<endl;

        par->datafP.slice(kf)=cg->x_k1_cxfmatp1ncpunpxy[kcpu][0];
    }
}

void beamformingCG3d(struct linerradon3d & par, int numi=9)
{
    tx_to_fx3d_linerradon3d(par);
    fmat wmat(par.npx,par.npy); 
    struct beamformingCG3d cg;
    struct beamforminginv3d par2;
    {
        struct linerradon3d * ppar;
        ppar=&par;
        linerradon(ppar[0]);
        beamforminginv3d_parset(ppar[0], par2);
        beamforminginv3d_getdigfmat(ppar[0], par2);
        beamformingCG3d_parset(cg, ppar[0], par2);
    }
    
    wmat=par.p_power/max(max(par.p_power));
    wmat+=par.dig_n;
    wmat=1.0/wmat;
    wmat=wmat/max(max(wmat));
    wmat*=par.dig_n;
    cg.fmatdigwnpxy[0]=wmat;
    cg.numi=numi;

    int pn(par.numthread),pnf1,pnf2,k;
    float dnf;
    thread *pcal;
    pcal=new thread[pn];
    dnf=float(par.nf2-par.nf1)/pn;

    for(k=0;k<pn-1;k++)
    {
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(beamformingCG3d_fthread,&cg,&par,pnf1,pnf2,k);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(beamformingCG3d_fthread,&cg,&par,pnf1,pnf2,k);

    for(k=0;k<pn;k++)
    {
        pcal[k].join();
    }

    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
    delete [] pcal;
}

fmat get_blackman_leftwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(j<=w)
            {
                n=j;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}
fmat get_blackman_rightwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(j>=n2-w-1)
            {
                n=n2-1-j;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

fmat get_blackman_downwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i>=n1-w-1)
            {
                n=n1-1-i;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

fmat get_blackman_upwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i<=w)
            {
                n=i;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

#endif
