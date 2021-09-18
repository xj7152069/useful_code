/************update in 2021.01.07**************/
/*
    
    
***********************************************/

#ifndef RADON3D_HPP
#define RADON3D_HPP

#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

#include <thread>
#include <future>
using namespace std;
using namespace arma;

#include "../xjc.h"

void cal_p_power_3d(struct linerradon3d & par);
void tx_to_fx3d(struct linerradon3d & par);
void tp_to_fp3d(struct linerradon3d & par);
void tx_to_fx2d(cx_fmat & datafx, fmat & data, int ny, int nf);
void fx_to_tx3d(struct linerradon3d & par);
void fx_to_tx2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz);
void fp_to_tp3d(struct linerradon3d & par);
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
    int nz,nx,ny,nf,npx,npy,nf1,nf2,numthread;
    float dz,dx,dy,df,dpx,dpy,p0x,p0y,px_center,py_center,dig_n;
    fcube data,realdataTP,realrebuildtx;
    cx_fcube datafx,rebuildfx,datafP,rebuildfp,dataTP,rebuildtx\
    ;
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
    par.dig_n*=par.nx*par.ny;
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
    fp_to_tp3d(par);
    par.realdataTP=real(par.dataTP);
}

////////////////////linerradon 3D transform/////////////////////////
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
    forA.fill(0.0);
    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        //cout<<"now is running kf = "<<kf<<endl;
        B=par->datafx.slice(kf);
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
                par->datafP(kpx,kpy,kf)=a(0,0);
            }
        }
    }
}

void linerradon(struct linerradon3d & par)
{
    tx_to_fx3d(par);

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

    fp_to_tp3d(par);
    par.realdataTP=real(par.dataTP);
    delete [] pcal;
}

////////////////////rebuild 3D data/////////////////////////
//重建数据，重建之前可以对数据按倾角去噪等
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
    forA.fill(0.0);

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        //cout<<"now is running kf = "<<kf<<endl;
        B=par->datafP.slice(kf);
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
    tp_to_fp3d(par);

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
   
    fx_to_tx3d(par);
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

void tp_to_fp3d(struct linerradon3d & par)
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

inline cx_fmat cx_fmatmul(cx_fmat & mat1, cx_fmat & mat2)
{
    int nz,nx;
    nz=mat1.n_rows;
    nx=mat1.n_cols;

    cx_fmat a(1,1,fill::zeros);
    int i,j;
    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
        {
            a(0,0)+=mat1(i,j)*mat2(i,j);
        }
    }

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

/////////////////////cal hessian old////////////////////////
/*
void beamforminginv3d_hessianget_thread_old(struct linerradon3d * par,\
 struct beamforminginv3d * par2, int pnf1, int pnf2)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    int numf(par[0].nf2-par[0].nf1);
    float df(par[0].df),dx(par[0].dx),dy(par[0].dy),\
        dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x),\
        p0y(par[0].p0y),dz(par[0].dz);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy);

    float ky1(-(ny*dy)/2.0);
    float kx1(-(nx*dx)/2.0);
    cx_fmat a(1,1);
    cx_fmat A(nx,ny),B(nx,ny);
    cx_fcube C(nx,ny,npx*npy);

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*df*(kf+par[0].nf1);  //!!!take care of when (par.nf1 != 0)
        cout<<"now is running kf = "<<kf<<endl;
        for(kpx=0;kpx<npx;kpx++)
        {
        for(kpy=0;kpy<npy;kpy++)
        {
            i=kpx*npy+kpy;
            C.slice(i).fill(0.0);
            C.slice(i).set_imag(w*par2[0].base_fmatp2npxnpy_nxny[kpx][kpy]);
            C.slice(i)=exp(C.slice(i));
        }
        }
        for(i=0;i<npx*npy;i++)
        {
            A=C.slice(i);
            for(j=0;j<=i;j++)
            {
                B=C.slice(j);
                B.set_imag(-imag(B));
                a=cx_fmatmul(A,B);
                par2[0].hessianinv_cxfmat_p1nf_npxynpxy[kf](i,j)=a(0,0);
                par2[0].hessianinv_cxfmat_p1nf_npxynpxy[kf](j,i)=a(0,0);
            }
            j--;  //take care of the for(;;)
            par2[0].hessianinv_cxfmat_p1nf_npxynpxy[kf](i,j)*=2; //hessian dig
        }
    }
}

void beamforminginv3d_hessianinv_old(struct linerradon3d & par,\
 struct beamforminginv3d & par2)
{
    int numf(par.nf2-par.nf1);
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dy(par.dy),\
        dpx(par.dpx),dpy(par.dpy),p0x(par.p0x),\
        p0y(par.p0y),dz(par.dz);
    int nx(par.nx),npx(par.npx),nf(par.nf),\
        ny(par.ny),npy(par.npy);
    int pn(par.numthread),pnf1,pnf2,k,kf;
    float dnf;
    par2.base_fmatp2npxnpy_nxny=new fmat*[par.npx];
    for(j=0;j<par.npx;j++)
    {
        par2.base_fmatp2npxnpy_nxny[j]=new fmat[par.npy];
        for(i=0;i<par.npy;i++)
        {
            par2.base_fmatp2npxnpy_nxny[j][i].zeros(par.nx,par.ny);
        }
    }
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
                    par2.base_fmatp2npxnpy_nxny[kpx][kpy]\
                        (kx,ky)=(fx*fpx+fy*fpy);
                    //A[kpx][kpy](kx,ky).imag(w*(fx*fpx+fy*fpy));
                }
            }
        }
    }
    thread *pcal;
    pcal=new thread[pn];
    dnf=float(par.nf2-par.nf1)/pn;

    for(k=0;k<pn-1;k++)
    {
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(beamforminginv3d_hessianget_thread_old,\
            &par,&par2,pnf1,pnf2);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(beamforminginv3d_hessianget_thread_old,\
            &par,&par2,pnf1,pnf2);

    for(k=0;k<pn;k++)
    {
        pcal[k].join();
    }
    delete [] pcal;
    for(kf=0;kf<numf;kf++)  
    {
        par2.hessianinv_cxfmat_p1nf_npxynpxy[kf]=\
            inv(par2.hessianinv_cxfmat_p1nf_npxynpxy[kf]);
    }
}
*/
/////////////////////////beamformingCG3d///////////////////////////
/*
struct beamformingCG3d
{
    int nx,nz,npx,npz,numi;
    fmat **base_fmatp2npxynxy, *fmatdigwnpxy;
    cx_fmat ***RH_cxfmatp3ncpunpxynxy;
    cx_fmat *b_cxfmatp1ncpunpxy;
    cx_fmat *r_k1_cxfmatp1ncpunpxy;
    cx_fmat *r_k2_cxfmatp1ncpunpxy;
    cx_fmat *p_k1_cxfmatp1ncpunpxy;
    cx_fmat *p_k2_cxfmatp1ncpunpxy;
    cx_fmat *x_k1_cxfmatp1ncpunpxy;
    cx_fmat *x_k2_cxfmatp1ncpunpxy;
};
struct Axdata
{
    int nx,nz,npx,npz;
    fmat *fmatp1digwnpxy;
    cx_fmat **cxfmatp2npxynxy,*cxfmatp1npxy;
    cx_fmat *outcxfmatp1npxy;
};
void beamformingCG3d_parset(struct beamformingCG3d & cg,
    struct linerradon3d & par)
{
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
    cg.numi=9;

    cg.RH_cxfmatp3ncpunpxynxy=new cx_fmat**[ncpu]; 
    for(i=0;i<ncpu;i++)
    {
        cg.RH_cxfmatp3ncpunpxynxy[i]=new cx_fmat*[npx];
        for(j=0;j<npx;j++)
        {
            cg.RH_cxfmatp3ncpunpxynxy[i][j]=new cx_fmat[npy];
            for(k=0;k<npy;k++)
            {
                cg.RH_cxfmatp3ncpunpxynxy[i][j][k].zeros(nx,ny);
            }
        }
    }
    cg.base_fmatp2npxynxy=new fmat*[npx];
    for(j=0;j<npx;j++)
    {
        cg.base_fmatp2npxynxy[j]=new fmat[npy];
        for(k=0;k<npy;k++)
        {
            cg.base_fmatp2npxynxy[j][k].zeros(nx,ny);
        }
    }
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
                    cg.base_fmatp2npxynxy[kpx][kpy](kx,ky)=(fx*fpx+fy*fpy);
                }
            }
        }
    }

    cg.fmatdigwnpxy=new fmat[1];
    cg.fmatdigwnpxy[0].zeros(npx,npy);

    cg.b_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.r_k1_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.r_k2_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.p_k1_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.p_k2_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.x_k1_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    cg.x_k2_cxfmatp1ncpunpxy=new cx_fmat[ncpu];
    for(i=0;i<ncpu;i++)
    {
        cg.b_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.r_k1_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.r_k2_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.p_k1_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.p_k2_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.x_k1_cxfmatp1ncpunpxy[i].zeros(npx,npy);
        cg.x_k2_cxfmatp1ncpunpxy[i].zeros(npx,npy);
    }
}

void cxfmatget_Ap(struct Axdata & p)
{
    int nz(p.nz),nx(p.nx),npx(p.npx),npz(p.npz);
    int knpx2,knpz2,knpx,knpz;
    cx_fmat a(1,1),b(1,1),cxfmatnxy(nx,nz);

    for(knpx=0;knpx<npx;knpx++)
    {
    for(knpz=0;knpz<npz;knpz++)
    {
        b.fill(0.0);
        for(knpx2=0;knpx2<npx;knpx2++)
        {
        for(knpz2=0;knpz2<npz;knpz2++)
        {
            cxfmatnxy.set_real(real(p.cxfmatp2npxynxy[knpx2][knpz2]));
            cxfmatnxy.set_imag(-imag(p.cxfmatp2npxynxy[knpx2][knpz2]));
            a=cx_fmatmul(p.cxfmatp2npxynxy[knpx][knpz],cxfmatnxy); 
            if(knpx==knpx2 && knpz==knpz2)
            {
                a+=p.fmatp1digwnpxy[0](knpx,knpz);
            }
            b(0,0)=b(0,0)+a(0,0)*p.cxfmatp1npxy[0](knpx2,knpz2);
        }
        }
        p.outcxfmatp1npxy[0](knpx,knpz)=b(0,0);
    }
    }
}

int beamformingCG3d_once(struct beamformingCG3d & cg, int kcpu)
{
    struct Axdata Ax;
    float bc;
    Ax.nz=cg.nz,Ax.nx=cg.nx,Ax.npx=cg.npx,Ax.npz=cg.npz;
    cx_fmat me(1,1),ak(1,1),bk(1,1),an(1,1),me_cxfmatnpxy;
    me_cxfmatnpxy.copy_size(cg.b_cxfmatp1ncpunpxy[kcpu]); 
    Ax.outcxfmatp1npxy=&me_cxfmatnpxy; //!!!!!
    Ax.fmatp1digwnpxy=cg.fmatdigwnpxy;
    Ax.cxfmatp2npxynxy=cg.RH_cxfmatp3ncpunpxynxy[kcpu];

//cal_ak || bc
    me_cxfmatnpxy.set_real(real(cg.r_k1_cxfmatp1ncpunpxy[kcpu]));
    me_cxfmatnpxy.set_imag(-imag(cg.r_k1_cxfmatp1ncpunpxy[kcpu]));
    ak=cx_fmatmul(cg.r_k1_cxfmatp1ncpunpxy[kcpu],me_cxfmatnpxy);
    Ax.cxfmatp1npxy=cg.p_k1_cxfmatp1ncpunpxy; 
    cxfmatget_Ap(Ax);
    me_cxfmatnpxy.set_imag(-imag(me_cxfmatnpxy)); 
    me=cx_fmatmul(cg.p_k1_cxfmatp1ncpunpxy[kcpu],me_cxfmatnpxy);
    if(real(me(0,0))==0)
    {return 1;}
    bc=real(ak(0,0))/real(me(0,0));
//cal_xk2
    cg.x_k2_cxfmatp1ncpunpxy[kcpu]=cg.x_k1_cxfmatp1ncpunpxy[kcpu]\
        +bc*cg.p_k1_cxfmatp1ncpunpxy[kcpu];
//cal_rk2
    me_cxfmatnpxy.set_imag(-imag(me_cxfmatnpxy)); 
    cg.r_k2_cxfmatp1ncpunpxy[kcpu]=cg.r_k1_cxfmatp1ncpunpxy[kcpu]\
        +bc*me_cxfmatnpxy[kcpu];
    if(sum(sum(abs(cg.r_k2_cxfmatp1ncpunpxy[kcpu])))==0)
    {return 1;}
//cal_bk || bc
    me_cxfmatnpxy.set_real(real(cg.r_k2_cxfmatp1ncpunpxy[kcpu]));
    me_cxfmatnpxy.set_imag(-imag(cg.r_k2_cxfmatp1ncpunpxy[kcpu]));
    bk=cx_fmatmul(cg.r_k2_cxfmatp1ncpunpxy[kcpu],me_cxfmatnpxy);
    me_cxfmatnpxy.set_real(real(cg.r_k1_cxfmatp1ncpunpxy[kcpu]));
    me_cxfmatnpxy.set_imag(-imag(cg.r_k1_cxfmatp1ncpunpxy[kcpu]));
    me=cx_fmatmul(cg.r_k1_cxfmatp1ncpunpxy[kcpu],me_cxfmatnpxy);
    bc=real(bk(0,0))/real(me(0,0));
//cal_pk2
    cg.p_k2_cxfmatp1ncpunpxy[kcpu]=cg.r_k2_cxfmatp1ncpunpxy[kcpu]\
        +bc*cg.p_k1_cxfmatp1ncpunpxy[kcpu];
//swap:
    cg.p_k1_cxfmatp1ncpunpxy[kcpu]=cg.p_k2_cxfmatp1ncpunpxy[kcpu];
    cg.r_k1_cxfmatp1ncpunpxy[kcpu]=cg.r_k2_cxfmatp1ncpunpxy[kcpu];
    cg.x_k1_cxfmatp1ncpunpxy[kcpu]=cg.x_k2_cxfmatp1ncpunpxy[kcpu];

    return 0;
}

void beamformingCG3d_fthread(struct beamformingCG3d * cg, \
    struct linerradon3d * par, int pnf1, int pnf2, int kcpu)
{
    int kf,kpx,kpy,kx,ky,i;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par->df),dx(par->dx),dy(par->dy),\
        dpx(par->dpx),dpy(par->dpy),p0x(par->p0x),\
        p0y(par->p0y),dz(par->dz);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy),cgstop;

    cx_fmat a(1,1),B(nx,ny),A(npx,npy);
    struct Axdata Ax;
    Ax.nz=cg->nz,Ax.nx=cg->nx,Ax.npx=cg->npx,Ax.npz=cg->npz;
    Ax.outcxfmatp1npxy=&A; //!!!!!

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        cout<<"now is running kf = "<<kf<<endl;
        B=par->datafx.slice(kf);
        for(kpx=0;kpx<npx;kpx++)
        {
        for(kpy=0;kpy<npy;kpy++)
        {
            cg->RH_cxfmatp3ncpunpxynxy[kcpu][kpx][kpy]\
                .set_imag(w*cg->base_fmatp2npxynxy[kpx][kpy]);
            cg->RH_cxfmatp3ncpunpxynxy[kcpu][kpx][kpy]\
                =exp(cg->RH_cxfmatp3ncpunpxynxy[kcpu][kpx][kpy]);
            a=cx_fmatmul(cg->RH_cxfmatp3ncpunpxynxy[kcpu][kpx][kpy],B);
            cg->b_cxfmatp1ncpunpxy[kcpu](kpx,kpy)=a(0,0);
        }
        }
        cg->x_k1_cxfmatp1ncpunpxy[kcpu]=cg->b_cxfmatp1ncpunpxy[kcpu];
        Ax.fmatp1digwnpxy=cg->fmatdigwnpxy;
        Ax.cxfmatp2npxynxy=cg->RH_cxfmatp3ncpunpxynxy[kcpu];
        Ax.cxfmatp1npxy=cg->x_k1_cxfmatp1ncpunpxy;
        cxfmatget_Ap(Ax); 
        cg->r_k1_cxfmatp1ncpunpxy[kcpu]=cg->b_cxfmatp1ncpunpxy[kcpu]-A;
        for(i=0;i<cg->numi;i++)
        {
            cout<<"CG times:"<<i<<endl;
            cgstop=beamformingCG3d_once(cg[0],kcpu);
            if(cgstop==1)
            {break;} 
        }
        par->datafP.slice(kf)=cg->x_k2_cxfmatp1ncpunpxy[kcpu];
    }
}

void beamformingCG3d(struct linerradon3d & par, int numi=9)
{
    tx_to_fx3d(par);
    fmat wmat(par.npx,par.npy); 
    struct beamformingCG3d cg;
    struct linerradon3d * ppar;
    ppar=&par;
    beamformingCG3d_parset(cg, ppar[0]);
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

    fp_to_tp3d(par);
    par.realdataTP=real(par.dataTP);
    delete [] pcal;
}
*/

#endif
