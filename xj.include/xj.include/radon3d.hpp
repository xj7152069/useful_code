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

fmat smoothdig(fmat dig,int l,int n);
void cal_p_power_3d(struct linerradon3d & par);
void tx_to_fx3d_linerradon3d(struct linerradon3d & par);
void tp_to_fp3d_linerradon3d(struct linerradon3d & par);
void fx_to_tx3d_linerradon3d(struct linerradon3d & par);
void fp_to_tp3d_linerradon3d(struct linerradon3d & par);
void linerradon(struct linerradon3d & par);
inline cx_fmat cx_fmatmul(cx_fmat & mat1, cx_fmat & mat2);
inline cx_dmat cx_hessianelement1d(float w,\
 int nx, float kx1, float dx, float dpx);
cx_fmat cx_hessianelement2d(float w,\
 int nx, float kx1, float dx, float dpx,\
 int ny, float ky1, float dy, float dpy);

/*beamforming/liner_radon变换 传递的参数：
nz（处理数据体Z/T方向采样点）、nx（处理数据体X方向采样点）、
ny（3D数据体测线数）、
nf（数据频率域变换的频率采样点）、np（Radon变换倾角采样点）、
dz（数据Z方向采样间隔）、dx（数据X方向采样间隔）、
df（数据变换后频率采样间隔）、dp（数据变换后倾角采样间隔）、
p0（数据变换后的中心倾角值）、data（原始数据）、
realdataTP（beamforming/liner_radon变换得到的数据实部）、
realrebuild（beamforming/liner_radon反变换重建的数据实部）、
realdatafft（beamforming/liner_radon变换得到的中间矩阵实部，用于检查算法）、
dataTP（beamforming/liner_radon变换得到的复数数据体）、
rebuild（beamforming/liner_radon反变换重建的复数数据体）、
dig_n（对角加权值）、

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
    float dig_n1,dig_n2,kerpar1;
    float x_center,y_center;
    //cx_fmat Acol1,Arow1;
};
//设置一组常用的初始化参数
void beamforming_parset(int nx, int ny, int nz, int nf,\
    struct linerradon3d & par)
{
    par.rulef1=5,par.rulef2=30;
    par.kerpar1=0;
    par.nx=nx,par.nz=nz,par.ny=ny,par.nf=nf;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.npx=par.nz,par.dpx=2*par.dz/par.nx/par.dx;
    par.npy=par.nz,par.dpy=2*par.dz/par.ny/par.dy;
    par.x_center=0;par.y_center=0;
    par.px_center=0;par.py_center=0;
    par.p0x=-par.dpx*int(par.npx/2)+par.px_center;
    par.p0y=-par.dpy*int(par.npy/2)+par.py_center;
    par.nf1=0;par.nf2=par.nf/2;
    par.numthread=1;
    par.dig_n=(1);
    par.dig_n1=(1);
    par.dig_n2=(1);
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
    //处理频率采样数，使其为2的幂次方
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
    par.dig_n1=par.dig_n1*par.nx*par.ny;
    par.dig_n2=par.dig_n2*par.nx*par.ny;
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
        par.x_coord(k,0)=k*par.dx-(int(par.nx/2)*par.dx)+par.x_center;
    }
    for(k=0;k<par.ny;k++)
    {
        par.y_coord(k,0)=k*par.dy-(int(par.ny/2)*par.dy)+par.y_center;
    }

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
    par.dataTP.fill(0.0);
    par.realdataTP.fill(0.0);
    par.rebuildfx.fill(0.0);
    par.rebuildtx.fill(0.0);
    par.realrebuildtx.fill(0.0);
}

////////////function: tx2fx or fx2tx/////////////
void tx_to_fx3d_radon3d_pthread(int i, struct linerradon3d * par)
{
    int j;
    cx_fmat datafx(par[0].ny,par[0].nf);
    fmat data(par[0].ny,par[0].nz);
    data=par[0].data.row(i);
    for(j=0;j<par[0].ny;j++)
    {
        datafx.row(j)=fft(data.row(j),par[0].nf);
    }
    par[0].datafx.row(i)=datafx;
}
void tx_to_fx3d_radon3d_thread(struct linerradon3d & par)
{
    int i,k;
    thread *pcal;
    pcal=new thread[par.numthread];
    
    for(i=0;i<(par.nx-par.numthread);i+=par.numthread)
    {
        for(k=0;k<par.numthread;k++)
        {
            pcal[k]=thread(tx_to_fx3d_radon3d_pthread,i+k,&par);
        }
        for(k=0;k<par.numthread;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(par.nx-par.numthread);i<(par.nx);i++)
    {
        k=i-(par.nx-par.numthread);
        pcal[k]=thread(tx_to_fx3d_radon3d_pthread,i,&par);
    }
    for(i=(par.nx-par.numthread);i<(par.nx);i++)
    {
        k=i-(par.nx-par.numthread);
        pcal[k].join();
    }

    for(i=0;i<par.ny;i++)
    {/*
        par.datafx.col(i)=get_blackman_leftwin2d\
            (par.datafx.col(i),par.rulef1);
        par.datafx.col(i)=get_blackman_rightwin2d\
            (par.datafx.col(i),30,par.rulef2);*/
    }
    
    fmat datapower(par.nf,1);
    for(i=0;i<par.nf;i++)
    {
        datapower(i,0)=sum(sum(abs(par.datafx.slice(i))));
    }
    datawrite(datapower,par.nf,1,"./fpower.bin");
}

void fx_to_tx3d_radon3d_pthread(int i, struct linerradon3d * par)
{
    cx_fmat data(par[0].ny,par[0].nf),datafx(par[0].ny,par[0].nf),\
        data2(par[0].ny,par[0].nz);
    int j;
    datafx=par[0].rebuildfx.row(i);
    for(j=0;j<par[0].ny;j++)
    {
        data.row(j)=ifft(datafx.row(j),par[0].nf);
    }
    for(j=0;j<par[0].nz;j++)
    {
        data2.col(j)=data.col(j);
    }
    par[0].rebuildtx.row(i)=data2;
    par[0].realrebuildtx.row(i)=real(data2);
}
void fx_to_tx3d_radon3d_thread(struct linerradon3d & par)
{
    int i,k;
    thread *pcal;
    pcal=new thread[par.numthread];
    
    for(i=0;i<(par.nx-par.numthread);i+=par.numthread)
    {
        for(k=0;k<par.numthread;k++)
        {
            pcal[k]=thread(fx_to_tx3d_radon3d_pthread,i+k,&par);
        }
        for(k=0;k<par.numthread;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(par.nx-par.numthread);i<(par.nx);i++)
    {
        k=i-(par.nx-par.numthread);
        pcal[k]=thread(fx_to_tx3d_radon3d_pthread,i,&par);
    }
    for(i=(par.nx-par.numthread);i<(par.nx);i++)
    {
        k=i-(par.nx-par.numthread);
        pcal[k].join();
    }
}

////////////beamforminginv3d use c++ thread/////////////////
struct beamforminginv3d
{
    fmat digw_fmat_npxnpy, **base_fmatp2npxnpy_nxny;
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

//kernel function of diagonally weighted matrix W
void beamforminginv3d_getdigfmat(struct linerradon3d & par,\
 struct beamforminginv3d & parinv)
{
    parinv.digw_fmat_npxnpy.fill(0.0);
    int k,i;
    float maxpower,minpower(0);
    cx_fmat cxmat(par.npx,par.npy);
    for(k=par.rulef1;k<par.rulef2;k++)
    {
        parinv.digw_fmat_npxnpy+=abs(cxmat=par.datafP.slice(k));
    } 

    maxpower=((parinv.digw_fmat_npxnpy.max()));
    minpower=((parinv.digw_fmat_npxnpy.min()));
    parinv.digw_fmat_npxnpy=(parinv.digw_fmat_npxnpy-minpower)/\
        (maxpower-minpower);
    parinv.digw_fmat_npxnpy+=sum(sum(parinv.digw_fmat_npxnpy))/\
        parinv.digw_fmat_npxnpy.n_elem;
    parinv.digw_fmat_npxnpy=1.0/parinv.digw_fmat_npxnpy;
    maxpower=((parinv.digw_fmat_npxnpy.max()));
    minpower=((parinv.digw_fmat_npxnpy.min()));
    parinv.digw_fmat_npxnpy=(parinv.digw_fmat_npxnpy-minpower)/\
        (maxpower-minpower);

    par.p_power=parinv.digw_fmat_npxnpy;
    if(par.npy>=5){
        par.p_power=fmatsmooth(par.p_power,par.npx,par.npy,1);
        maxpower=((par.p_power.max()));
        minpower=((par.p_power.min()));
        par.p_power=(par.p_power-minpower)/\
            (maxpower-minpower);
        for(k=0;k<par.npy;k++){
            for(i=0;i<par.npx;i++){
                par.p_power(i,k)=1.0/(1.0+exp(100*(0.5-par.p_power(i,k))));
            }
            datawrite(par.p_power,par.npx,par.npy,"dig.bin");
        }
    }
    else{
        fmat digline(par.npx,1);
        for(k=0;k<par.npy;k++){
            par.p_power.col(k)=smoothdig(par.p_power.col(k),par.npx,99);
        }
        maxpower=((par.p_power.max()));
        minpower=((par.p_power.min()));
        par.p_power=(par.p_power-minpower)/\
            (maxpower-minpower);
        par.p_power+=par.kerpar1;
        for(k=0;k<par.npy;k++){
            for(i=0;i<par.npx;i++){
                par.p_power(i,k)=1.0/(1.0+exp(10*(0.5-par.p_power(i,k))));
            }
            datawrite(digline=par.p_power.col(k),par.npx,1,"dig.bin");
        }
    }

    //get diagonally weighted matrix W
    parinv.digw_fmat_npxnpy=par.dig_n2+par.dig_n1*par.p_power;
    par.p_power=parinv.digw_fmat_npxnpy;
}

//get mat LTL
void beamforminginv3d_hessianget_thread(struct linerradon3d * par,\
 struct beamforminginv3d * par2, int pncpu, int pnf)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par[0].df),dx(par[0].dx),dy(par[0].dy),\
        dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x),\
        p0y(par[0].p0y),dz(par[0].dz),dig_n(par[0].dig_n);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy);

    float ky1(-(int(ny/2)*dy));
    float kx1(-(int(nx/2)*dx));
    cx_fmat a(1,1),A(nx,ny),B(nx,ny);

    kf=pnf;
    if(kf<par[0].nf2) {
        w=2.0*pi*df*(kf); 
        for(kpx=0;kpx<npx;kpx++){
        for(kpy=0;kpy<npy;kpy++){
            i=kpx*npy+kpy;
            //A=basex.col(kpx)*basey.col(kpy).st(); //
            float fpx1(par[0].px_coord(kpx,0));
            float fpy1(par[0].py_coord(kpy,0));
            for(kpx2=0;kpx2<npx;kpx2++){
            for(kpy2=0;kpy2<npy;kpy2++){
                j=kpx2*npy+kpy2;
                //B=basex.col(kpx2)*basey.col(kpy2).st();
                //B.set_imag(-imag(B));
                float fpx2(par[0].px_coord(kpx2,0));
                float fpy2(par[0].py_coord(kpy2,0));
                a=cx_hessianelement2d(w,nx,kx1,dx,fpx1-fpx2,\
                    ny,ky1,dy,fpy1-fpy2);
                //a=cx_fmatmul(A,B);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,j)=a(0,0);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](j,i)=a(0,0); 
            }
            }        
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,i)+=par->p_power(kpx,kpy);
            kpx2=npx;kpy2=npy;
        }
        } 
        
        par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]=\
            inv(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]);
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

//inv function of radon
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
    beamforminginv3d_parset(ppar[0], par2);
    beamforminginv3d_getdigfmat(ppar[0], par2);

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

////////////////////////////////////////
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

/*
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
*/
///////////////linerradon 3D transform: slant stack///////////////
void linerradon_fthread(struct linerradon3d * par,int pnf1, int pnf2)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    float df(par->df),dx(par->dx),dy(par->dy),\
        dpx(par->dpx),dpy(par->dpy),p0x(par->p0x),\
        p0y(par->p0y),dz(par->dz);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy);

    cx_fmat forA(nx,ny),A(nx,ny),a(1,1),B(nx,ny),hesscxmat(npx*npy,npx*npy);
    cx_fmat basey(ny,npy,fill::zeros),basex(nx,npx,fill::zeros);
    cx_fmat ap(nx,ny,fill::zeros),bp(nx,ny,fill::zeros);
    fmat hessfmat(npx*npy,npx*npy);

    forA.fill(0.0);
    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
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
                a=cx_fmatmul(A,B);
                par->datafP(kpx,kpy,kf)=a(0,0);
            }
        }
    }
}

void linerradon(struct linerradon3d & par)
{
    tx_to_fx3d_radon3d_thread(par);
    
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
    par.realdataTP=real(par.dataTP);
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
   
    //fx_to_tx3d_linerradon3d(par);
    fx_to_tx3d_radon3d_thread(par);
    par.realrebuildtx=real(par.rebuildtx);
    delete [] pcal;
}

///////////////////////////////////////////////////////////
inline void tx_to_fx2d_pthread(cx_fmat * datafx, fmat * data, int n, int nf)
{
    datafx[0].row(n)=fft(data[0].row(n),nf);
}
void tx_to_fx2d(cx_fmat & datafx, fmat & data, int ny, int nf, int numthread=1)
{
    int i,k;
    thread *pcal;
    pcal=new thread[numthread];
    cx_fmat * pdatafx;fmat * pdata;
    pdatafx=&datafx;pdata=&data;
    
    for(i=0;i<(ny-numthread);i+=numthread)
    {
        for(k=0;k<numthread;k++)
        {
            pcal[k]=thread(tx_to_fx2d_pthread,pdatafx,pdata,i+k,nf);
        }
        for(k=0;k<numthread;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(ny-numthread);i<(ny);i++)
    {
        k=i-(ny-numthread);
        pcal[k]=thread(tx_to_fx2d_pthread,pdatafx,pdata,i,nf);
    }
    for(i=(ny-numthread);i<(ny);i++)
    {
        k=i-(ny-numthread);
        pcal[k].join();
    }
}
inline void fx_to_tx2d_pthread(cx_fmat * data2, cx_fmat * datafx, int n, int nf)
{
    data2[0].row(n)=ifft(datafx[0].row(n),nf);
}
void fx_to_tx2d(cx_fmat & data, cx_fmat & datafx, int ny, int nz,int numthread=1)
{
    int nf(datafx.n_cols);
    int i,k;
    thread *pcal;
    pcal=new thread[numthread];
    cx_fmat data2(ny,nf);
    cx_fmat * pdatafx;
    pdatafx=&datafx;
    
    for(i=0;i<(ny-numthread);i+=numthread)
    {
        for(k=0;k<numthread;k++)
        {
            pcal[k]=thread(fx_to_tx2d_pthread,&data2,pdatafx,i+k,nf);
        }
        for(k=0;k<numthread;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(ny-numthread);i<(ny);i++)
    {
        k=i-(ny-numthread);
        pcal[k]=thread(fx_to_tx2d_pthread,&data2,pdatafx,i,nf);
    }
    for(i=(ny-numthread);i<(ny);i++)
    {
        k=i-(ny-numthread);
        pcal[k].join();
    }
    
    for(i=0;i<nz;i++)
    {
        data.col(i)=data2.col(i);
    } 
}
void tp_to_fp3d_linerradon3d(struct linerradon3d & par)
{
    int i;
    fmat data(par.npx,par.nz);
    cx_fmat datafp(par.npx,par.nf);

    for(i=0;i<par.npy;i++)
    {
        tx_to_fx2d(datafp, data=par.realdataTP.col(i),\
            par.npx, par.nf,par.numthread);
        par.datafP.col(i)=datafp;
    }    
}

void tx_to_fx3d_linerradon3d(struct linerradon3d & par)
{
    int i;
    fmat data(par.nx,par.nz),datapower(par.nf,1,fill::zeros);
    cx_fmat datafx(par.nx,par.nf);

    for(i=0;i<par.ny;i++)
    {
        tx_to_fx2d(datafx, data=par.data.col(i),\
            par.nx, par.nf,par.numthread);
        par.datafx.col(i)=datafx;
    }
    for(i=0;i<par.nf;i++)
    {
        datapower(i,0)=sum(sum(abs(par.datafx.slice(i))));
    }
    datawrite(datapower,par.nf,1,"./fpower.bin");
}
void fx_to_tx3d_linerradon3d(struct linerradon3d & par)
{
    int i;
    cx_fmat data(par.nx,par.nz);
    cx_fmat datafx(par.nx,par.nf);

    for(i=0;i<par.ny;i++)
    {
        fx_to_tx2d(data,datafx=par.rebuildfx.col(i),\
            par.nx, par.nz,par.numthread);
        par.rebuildtx.col(i)=data;
    } 
}
void fp_to_tp3d_linerradon3d(struct linerradon3d & par)
{
    int i;
    cx_fmat data(par.npx,par.nz);
    cx_fmat datafp(par.npx,par.nf);

    for(i=0;i<par.npy;i++)
    {
        fx_to_tx2d(data,datafp=par.datafP.col(i),\
            par.npx, par.nz,par.numthread);
        par.dataTP.col(i)=data;
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
    return a;
}

cx_fmat cx_hessianelement2d(float w,\
 int nx, float kx1, float dx, float dpx,\
 int ny, float ky1, float dy, float dpy)
{
    cx_dmat a(1,1,fill::zeros),b(1,1,fill::zeros);
    float n2(0.000000000000);
    float pdx(abs(dpx)),pdy(abs(dpy));

    if(pdx>0 && pdy>0){
        a=cx_hessianelement1d(w, nx, kx1, dx, dpx);
        b=cx_hessianelement1d(w, ny, ky1, dy, dpy);
    }
    else if(pdx==0 && pdy>n2){
        a(0,0).real(nx);
        b=cx_hessianelement1d(w, ny, ky1, dy, dpy);
    }
    else if(pdx>0 && pdy==0){
        a=cx_hessianelement1d(w, nx, kx1, dx, dpx);
        b(0,0).real(ny);
    }
    else if(pdx==0 && pdy==0){
        a(0,0).real(nx);
        b(0,0).real(ny);
    }
    
    a=a*b;
    cx_fmat af(1,1,fill::zeros);
    af(0,0).real(a(0,0).real());
    af(0,0).imag(a(0,0).imag());
    return af;
}

inline cx_dmat cx_hessianelement1d(float w,\
 int nx, float kx1, float dx, float dpx)
{
    cx_dmat a(1,1,fill::zeros);
    cx_dmat b(1,1,fill::zeros);
    cx_dmat c(1,1,fill::zeros);
    cx_dmat d(1,1,fill::zeros);
    cx_dmat o(1,1,fill::zeros);
    o(0,0).real(1.0);

    a(0,0).imag(w*dpx*kx1);
    b=exp(a);
    a(0,0).imag(w*dpx*nx*dx);
    c=exp(a)-o;
    a(0,0).imag(w*dpx*dx);
    d=exp(a)-o;

    a=b*c/d;
    return a;
}

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
/////////////beamformingCG3d: not stable !/////////////

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

fmat smoothdig(fmat dig,int l,int n)
{
    fmat d1(l,1),d2(l,1);
    d1=dig,d2=dig;
    int j,k;
    for(j=0;j<n;j++)
    {
        for(k=1;k<l-1;k++)
        {
            d2(k,0)=(0.5*d1(k-1,0)+0.5*d1(k+1,0)+d1(k,0))/2;
        }
        d1=d2;
    }
    dig=d2;
    return dig;
}
#endif
