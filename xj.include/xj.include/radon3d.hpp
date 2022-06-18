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
// 2-D Linear Radon Transform in the frequence - space domain.
int Beamforming_LS_2D(float **trace, int ntrace, int ns, float dt,\
 float *coor, float x0, float **tauppanel, int npsr,float psrmin, float dpsr,\
 float fmax,float frule,int ncpu, bool regularization,\
 float factor_L2, int iterations_num, float residual_ratio);
/* discription.
This function will do the local slant-stack over local traces .
float
**trace;
the INPUT local seismic-gather, i.e. trace[ntrace][ns].
int ntrace;
the trace-number of the INPUT gather.
int ns;
the trace length.
float dt;
time- sampling interval, (unit of dt should be second) ,
float *coor;
coor[ntrace] stores the offset of each trace .
float x0;
the beam-center coordinate .
float **tauppanel; the OUTPUT tau-p spectrum, i.e. tauppanel[npsr][ns].
int npsr;
the ray- parameter number of tau-p panel.
float psrmin;
the minimum ray-parameter, (unit is s/m) .
float dpsr;
the interval of ray-parameter， (unit is s/m) .*/
///////////////////////////////////////////////////////////////////
fmat smoothdig(fmat dig,int l,int n);
void cal_p_power_3d(struct linerradon3d & par);
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
    float dz,dx,dy,df,dpx,dpy,p0x,p0y,dig_n;
/*****************************************************/
    fcube data,realdataTP,realrebuildtx;
    cx_fcube datafx,rebuildfx,datafP;
/*****************************************************/
    fmat p_power,py_coord,px_coord,x_coord,y_coord,coordx3d,coordy3d;
    char allAreal[99],allAimag[99];
    float dig_n1,dig_n2,kerpar1;
    float coordx0,coordy0;
    bool regularization;
    //cx_fmat Acol1,Arow1;
};
//设置一组常用的初始化参数
void beamforming_parset(int nx, int ny, int nz,\
    struct linerradon3d & par)
{
    par.rulef1=5,par.rulef2=30;
    par.kerpar1=0;
    par.nx=nx,par.nz=nz,par.ny=ny,par.nf=nz;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.nf1=0;par.nf2=par.nf/2;
    par.rulef1=0;par.rulef2=par.nf/2;
    par.numthread=1;
    par.dig_n=(1);
    par.dig_n1=(1);
    par.dig_n2=(1);

    par.dig_n=par.dig_n*par.nx*par.ny;
    par.dig_n1=par.dig_n1*par.nx*par.ny;
    par.dig_n2=par.dig_n2*par.nx*par.ny;
    par.coordx3d.zeros(par.nx,par.ny);
    par.coordy3d.zeros(par.nx,par.ny);
    par.x_coord.zeros(par.nx,1);
    par.y_coord.zeros(par.ny,1);
    int k;
    
    par.dig_n1=0.0;  //L1
    par.kerpar1=0.0;     //kernel function par
    par.regularization=true;
} 

void beamforming_parupdate(struct linerradon3d & par)
{ 
    par.df=1.0/par.dz/par.nf;
    par.data.zeros(par.nx,par.ny,par.nz);
    par.datafP.zeros(par.npx,par.npy,par.nf);
    //par.datafx.zeros(par.nx,par.ny,par.nf);
    par.realdataTP.zeros(par.npx,par.npy,par.nz);
    //par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    //par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
    par.p_power.zeros(par.npx,par.npy);
    par.nf1=1;
    par.rulef1=1;
    par.nf2=int(par.nf2/par.df);
    par.rulef2=int(par.rulef2/par.df);
    if(par.nf2>par.nf/2){
        par.nf2=par.nf/2;
    }
    if(par.rulef2>par.nf2/2){
        par.rulef2=par.nf2/2;
    }

    par.px_coord.zeros(par.npx,1);
    par.py_coord.zeros(par.npy,1);
    
    int k;
    par.coordx0=(-par.dx*par.nx/2),par.coordy0=(-par.dy*par.ny/2);
    //par.p0x=(-par.dpx*par.npx/2),par.p0y=(-par.dpy*par.npy/2);
    for(k=0;k<par.npx;k++){
        par.px_coord(k,0)=k*par.dpx+par.p0x;
    }
    for(k=0;k<par.npy;k++){
        par.py_coord(k,0)=k*par.dpy+par.p0y;
    }
    for(k=0;k<par.nx;k++){
        par.x_coord(k,0)=k*par.dx+par.coordx0;
    }
    for(k=0;k<par.ny;k++){
        par.y_coord(k,0)=k*par.dy+par.coordy0;
    }
    if(par.dpx==0){
        par.dpx=0.00001;
    }
    if(par.dpy==0){
        par.dpy=par.dpx;
    }
}

void beamforming_cleardata(struct linerradon3d & par)
{

}

////////////beamforminginv3d use c++ thread/////////////////
struct beamforminginv3d
{
    fmat digw_fmat_npxnpy;
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
    for(k=0;k<ncpu;k++){
        parinv.hessianinv_cxfmat_p1ncpu_npxynpxy[k].zeros(n,n);
    }
}

//////////////////////////////////////////////////////////////////////////////////////
//kernel function of diagonally weighted matrix W
void beamforminginv3d_getdigfmat(struct linerradon3d & par,\
 fmat &digw_fmat_npxnpy)
{
    digw_fmat_npxnpy.fill(0.0);
    int k,i;
    float maxpower,minpower(0);
    cx_fmat cxmat(par.npx,par.npy);
    for(k=par.rulef1;k<par.rulef2;k++){
        digw_fmat_npxnpy+=abs(cxmat=par.datafP.slice(k));
    } 

    maxpower=((digw_fmat_npxnpy.max()));
    minpower=((digw_fmat_npxnpy.min()));
    digw_fmat_npxnpy=(digw_fmat_npxnpy-minpower)/\
        (maxpower-minpower);
    digw_fmat_npxnpy+=sum(sum(digw_fmat_npxnpy))/\
        digw_fmat_npxnpy.n_elem;
    digw_fmat_npxnpy=1.0/digw_fmat_npxnpy;
    maxpower=((digw_fmat_npxnpy.max()));
    minpower=((digw_fmat_npxnpy.min()));
    digw_fmat_npxnpy=(digw_fmat_npxnpy-minpower)/\
        (maxpower-minpower);

    par.p_power=digw_fmat_npxnpy;
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
            //datawrite(par.p_power,par.npx,par.npy,"dig.bin");
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
            //datawrite(digline=par.p_power.col(k),par.npx,1,"dig.bin");
        }
    }
    //get diagonally weighted matrix W
    digw_fmat_npxnpy=par.dig_n2+par.dig_n1*par.p_power;
    par.p_power=digw_fmat_npxnpy;
}

//get mat LTL
void beamforminginv3d_hessianget_thread(struct linerradon3d * par,\
 struct beamforminginv3d * par2, int pncpu, int pnf, bool regularization,\
 bool doinv=true)
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
            fmat hess;

    kf=pnf;
    if(kf<par[0].nf2){
        w=2.0*pi*df*(kf); 
        if(regularization)
        {
            for(kpx=0;kpx<npx;kpx++){
            for(kpy=0;kpy<npy;kpy++){
                i=kpx*npy+kpy;
                float fpx1(par[0].px_coord(kpx,0));
                float fpy1(par[0].py_coord(kpy,0));
            for(kpx2=0;kpx2<npx;kpx2++){
            for(kpy2=0;kpy2<npy;kpy2++){
                j=kpx2*npy+kpy2;
                float fpx2(par[0].px_coord(kpx2,0));
                float fpy2(par[0].py_coord(kpy2,0));
                a=cx_hessianelement2d(w,nx,kx1,dx,fpx1-fpx2,\
                    ny,ky1,dy,fpy1-fpy2);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,j)=a(0,0);
                //par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](j,i)=a(0,0); 
            }}        
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,i)\
                +=par->p_power(kpx,kpy);
            kpx2=npx;kpy2=npy;
            }}
        }
        else
        {
            cx_fmat hessian_small(npx*2-1,npy*2-1);
            cx_fcube basey3d;
            cx_fmat basex,basey;
            cx_fmat basepy0,basepx0;
            cx_fmat basedpy,basedpx;
            float fpx1,fpy1,fpx2,fpy2,fpxmax,fpxmin,fpymax,fpymin;
            int nkpx,nkpy;
            fpxmax=par[0].px_coord.max()-par[0].px_coord.min();
            fpxmin=par[0].px_coord.min()-par[0].px_coord.max();
            fpymax=par[0].py_coord.max()-par[0].py_coord.min();
            fpymin=par[0].py_coord.min()-par[0].py_coord.max();

            basey.zeros(nx,ny);basex.zeros(nx,ny);
            basey3d.zeros(nx,ny,2*npy);
            basepy0.zeros(nx,ny);basepx0.zeros(nx,ny);
            basedpy.zeros(nx,ny);basedpx.zeros(nx,ny);
            basepx0.set_imag(w*par->coordx3d*(fpxmin-dpx));  
            basepy0.set_imag(w*par->coordy3d*(fpymin));
            basedpx.set_imag(w*par->coordx3d*dpx);  
            basedpy.set_imag(w*par->coordy3d*dpy);
            basedpy=exp(basedpy);basedpx=exp(basedpx);
            basepy0=exp(basepy0);basepx0=exp(basepx0);

            basex=basepx0;
            for(fpx2=fpxmin;fpx2<=fpxmax;fpx2+=dpx){
                kpx2=round(fpx2/dpx)+npx-1;
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    basex(kx,ky)=basex(kx,ky)*basedpx(kx,ky);
                }}
                if(fpx2==fpxmin)
                {basey3d.slice(0)=basepy0;}
            for(fpy2=fpymin;fpy2<=fpymax;fpy2+=dpy){
                kpy2=round(fpy2/dpy)+npy-1;
                a=cx_fmatmul(basex,basey=basey3d.slice(kpy2));
                hessian_small(kpx2,kpy2)=a(0,0);
                if(fpx2==fpxmin){
                    for(kx=0;kx<nx;kx++){
                    for(ky=0;ky<ny;ky++){
                    basey3d(kx,ky,kpy2+1)=basey3d(kx,ky,kpy2)*basedpy(kx,ky);
                }}}
            }}

            for(kpx=0;kpx<npx;kpx++){
            for(kpy=0;kpy<npy;kpy++){
                i=kpx*npy+kpy;
                fpx1=(par[0].px_coord(kpx,0));
                fpy1=(par[0].py_coord(kpy,0));
            for(kpx2=0;kpx2<npx;kpx2++){
            for(kpy2=0;kpy2<npy;kpy2++){
                j=kpx2*npy+kpy2;
                fpx2=(par[0].px_coord(kpx2,0));
                fpy2=(par[0].py_coord(kpy2,0));
                nkpx=round((fpx1-fpx2)/dpx)+npx-1;
                nkpy=round((fpy1-fpy2)/dpy)+npy-1;
                a(0,0)=hessian_small(nkpx,nkpy);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,j)=a(0,0);
                //par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](j,i)=a(0,0);
            }}
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,i)\
                +=par->p_power(kpx,kpy);
            }}
        }
        /*
        if(kf==30){
            datawrite(hess=real(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]),
                npx*npy,npx*npy,"hess2.bin");
        }*/

        if(doinv){
        par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]=\
            inv(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]);
        }
    }
}

void beamforminginv3d_beamget_thread(struct linerradon3d * par,\
 struct beamforminginv3d * par2, int pncpu, int pnf,\
 bool regularization,bool *finish_thread)
{
    beamforminginv3d_hessianget_thread(\
        par,par2,pncpu,pnf,regularization,true);
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
    cout<<"now is inv: "<<kf<<endl;
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
    finish_thread[0]=true;
}

//inv function of radon
void beamforminginv3d(struct linerradon3d & par, bool regularization=true)
{
    int numf(par.nf2-par.nf1),ncpu(par.numthread),kf,kcpu;
    regularization=par.regularization;

    struct beamforminginv3d par2;
    struct linerradon3d * ppar;
    ppar=&par;
    beamforminginv3d_parset(ppar[0], par2);
    beamforminginv3d_getdigfmat(ppar[0], par2.digw_fmat_npxnpy);

    bool* finish_thread;
    finish_thread=new bool[ncpu];
    thread *pcal;
    pcal=new thread[ncpu];
    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(beamforminginv3d_beamget_thread,&par,&par2,\
                kcpu,kcpu+par.nf1,regularization,&(finish_thread[kcpu]));
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable()&&kf<par.nf2&&finish_thread[kcpu]){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(beamforminginv3d_beamget_thread,&par,&par2,\
                kcpu,kf,regularization,&(finish_thread[kcpu]));
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }

    delete[] pcal;
    delete[] finish_thread;

    fx2tx_3d_thread(par.realdataTP,par.datafP,par.numthread);
}

///////////////linerradon 3D transform: slant stack///////////////
void linerradon_fthread(struct linerradon3d * par,int pnf1, int pnf2,\
 bool *finish_thread)
{
    int kf,kpx,kpy,kx,ky,i,j;//cout<<"ok"<<endl;
    float fx,fy,fpx,fpy;
    float w,pi(3.1415926);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy);

    cx_fmat a(1,1),A,B;
    cx_fcube basey3d;
    cx_fmat basex,basey;
    cx_fmat basepy0,basepx0;
    cx_fmat basedpy,basedpx;

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        if(par->regularization){
            A.zeros(nx,ny);B.zeros(nx,ny);
            basey.zeros(ny,npy);basex.zeros(nx,npx);
            basey.set_imag(w*(par->y_coord*par->py_coord.t()));  
            basex.set_imag(w*(par->x_coord*par->px_coord.t()));  
            basey=exp(basey);basex=exp(basex);
            B=par->datafx.slice(kf); 
            for(kpx=0;kpx<npx;kpx++){
                for(kpy=0;kpy<npy;kpy++){
                    A=basex.col(kpx)*basey.col(kpy).st();
                    a=cx_fmatmul(A,B);
                    par->datafP(kpx,kpy,kf)=a(0,0);
                }
            }
        }
        else{
            basey3d.set_size(nx,ny,npy+1);
            basepy0.zeros(nx,ny);basepx0.zeros(nx,ny);
            basedpy.zeros(nx,ny);basedpx.zeros(nx,ny);
            basepx0.set_imag(w*par->coordx3d*(par->px_coord(0,0)));  
            basepy0.set_imag(w*par->coordy3d*(par->py_coord(0,0)));
            basedpx.set_imag(w*par->coordx3d*par->dpx);  
            basedpy.set_imag(w*par->coordy3d*par->dpy);
            basedpy=exp(basedpy);basedpx=exp(basedpx);
            basepy0=exp(basepy0);basepx0=exp(basepx0);
            B.set_size(nx,ny);
            B=par->datafx.slice(kf); 

            kpx=0;{
                {basey3d.slice(0)=basepy0;}
            kpy=0;{
                cx_float fzero;
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    B(kx,ky)*=basepx0(kx,ky);
                    fzero+=B(kx,ky)*basey3d(kx,ky,kpy);
                    {basey3d(kx,ky,kpy+1)=basey3d(kx,ky,kpy)*basedpy(kx,ky);}
                }}
                par->datafP(kpx,kpy,kf)=fzero;
                }
                for(kpy=1;kpy<npy-1;kpy++){
                cx_float fzero;
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    fzero+=B(kx,ky)*basey3d(kx,ky,kpy);
                    {basey3d(kx,ky,kpy+1)=basey3d(kx,ky,kpy)*basedpy(kx,ky);}
                }}
                par->datafP(kpx,kpy,kf)=fzero;
                }
                kpy=npy-1;
                if(kpy>0){
                cx_float fzero;
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    fzero+=B(kx,ky)*basey3d(kx,ky,kpy);
                }}
                par->datafP(kpx,kpy,kf)=fzero;
                }
            }

            for(kpx=1;kpx<npx;kpx++){
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    B(kx,ky)*=basedpx(kx,ky);
                }}
            for(kpy=0;kpy<npy;kpy++){
                cx_float fzero;
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    fzero+=B(kx,ky)*basey3d(kx,ky,kpy);
                }}
                par->datafP(kpx,kpy,kf)=fzero;
            }}
        }
    }
    finish_thread[0]=true;
}

void linerradon(struct linerradon3d & par)
{
    par.datafx.set_size(par.nx,par.ny,par.nf);
    tx2fx_3d_thread(par.datafx,par.data,par.numthread);
    
    int ncpu(par.numthread),pnf1,pnf2,k,kcpu,kf;
    float dnf;
    bool* finish_thread;
    finish_thread=new bool[ncpu];
    thread *pcal;
    pcal=new thread[ncpu];
    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(linerradon_fthread,\
            &par,kcpu+par.nf1,kcpu+par.nf1+1,&(finish_thread[kcpu]));
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable() && kf<par.nf2 && finish_thread[kcpu]){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(linerradon_fthread,\
                    &par,kf,kf+1,&(finish_thread[kcpu]));
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }
    
    //fp_to_tp3d_linerradon3d(par);
    fx2tx_3d_thread(par.realdataTP,par.datafP,par.numthread);
    par.datafx.set_size(1,1,1);

    delete[] pcal;
    delete[] finish_thread;
}

////////////////////rebuild 3D data/////////////////////////
void rebuildsignal_fthread(struct linerradon3d * par,int pnf1, int pnf2,\
 bool *finish_thread)
{
    int kf,kpx,kpy,kx,ky;//cout<<"ok"<<endl;
    float w,pi(3.1415926);
    int nx(par->nx),npx(par->npx),nf(par->nf),\
        ny(par->ny),npy(par->npy);

    cx_fmat A,a(1,1),B;
    cx_fcube basey3d;
    cx_fmat basey,basex;
    cx_fmat basepy0,basepx0;
    cx_fmat basedpy,basedpx;

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        w=2.0*pi*par->df*kf; 
        if(par->regularization){
            A.zeros(npx,npy);B.zeros(npx,npy);
            basey.zeros(ny,npy);basex.zeros(nx,npx);
            basey.set_imag(w*(par->y_coord*par->py_coord.t()));  
            basex.set_imag(w*(par->x_coord*par->px_coord.t()));  
            basey=exp(basey);basex=exp(basex);
            B=par->datafP.slice(kf);
            for(kx=0;kx<nx;kx++){
            for(ky=0;ky<ny;ky++){
                A=basex.row(kx).st()*basey.row(ky);
                A.set_imag(-imag(A));
                a=cx_fmatmul(A,B);
                par->rebuildfx(kx,ky,kf)=a(0,0);
            }}
        }
        else{
            basey3d.set_size(nx,ny,npy+1);
            basex.set_size(nx,ny);
            basepy0.zeros(nx,ny);basepx0.zeros(nx,ny);
            basedpy.zeros(nx,ny);basedpx.zeros(nx,ny);
            basepx0.set_imag(-w*par->coordx3d*(par->px_coord(0,0)));  
            basepy0.set_imag(-w*par->coordy3d*(par->py_coord(0,0)));
            basedpx.set_imag(-w*par->coordx3d*par->dpx);  
            basedpy.set_imag(-w*par->coordy3d*par->dpy);
            basedpy=exp(basedpy);basedpx=exp(basedpx);
            basepy0=exp(basepy0);basepx0=exp(basepx0);
            B.set_size(npx,npy);
            A.zeros(nx,ny);
            B=par->datafP.slice(kf);
            kpx=0;{
                kpy=0;{
                basey3d.slice(0)=basepy0;
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    basey3d(kx,ky,kpy)*=basepx0(kx,ky);
                    A(kx,ky)+=(B(kpx,kpy)*basey3d(kx,ky,kpy));
                    {basey3d(kx,ky,kpy+1)=basey3d(kx,ky,kpy)*basedpy(kx,ky);}
                }}}
                for(kpy=1;kpy<npy-1;kpy++){
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    A(kx,ky)+=(B(kpx,kpy)*basey3d(kx,ky,kpy));
                    {basey3d(kx,ky,kpy+1)=basey3d(kx,ky,kpy)*basedpy(kx,ky);}
                }}}
                kpy=npy-1;
                if(kpy>0){
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    A(kx,ky)+=(B(kpx,kpy)*basey3d(kx,ky,kpy));
                }}}
            }
            basex=basedpx;
            for(kpx=1;kpx<npx-1;kpx++){
                for(kpy=0;kpy<npy;kpy++){
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    A(kx,ky)+=(B(kpx,kpy)*basex(kx,ky)*basey3d(kx,ky,kpy));
                }}}
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    basex(kx,ky)*=basedpx(kx,ky);
                }}
            }
            kpx=npx-1;{
                for(kpy=0;kpy<npy;kpy++){
                for(kx=0;kx<nx;kx++){
                for(ky=0;ky<ny;ky++){
                    A(kx,ky)+=(B(kpx,kpy)*basex(kx,ky)*basey3d(kx,ky,kpy));
                }}}
            }
            par->rebuildfx.slice(kf)=A;
        }
    }
    finish_thread[0]=true;
}

void rebuildsignal(struct linerradon3d & par)
{
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    tx2fx_3d_thread(par.datafP,par.realdataTP,par.numthread);

    int ncpu(par.numthread),pnf1,pnf2,k,kf,kcpu;
    float dnf;
    bool* finish_thread;
    finish_thread=new bool[ncpu];
    thread *pcal;
    pcal=new thread[ncpu];
    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(rebuildsignal_fthread,\
            &par,kcpu+par.nf1,kcpu+par.nf1+1,&(finish_thread[kcpu]));
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable()&&kf<par.nf2&&finish_thread[kcpu]){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(rebuildsignal_fthread,\
                    &par,kf,kf+1,&(finish_thread[kcpu]));
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }
    par.realrebuildtx.set_size(par.nx,par.ny,par.nz);
    fx2tx_3d_thread(par.realrebuildtx,par.rebuildfx,par.numthread);
    par.realrebuildtx=4*par.realrebuildtx;
    par.rebuildfx.set_size(1,1,1);

    delete[] pcal;
    delete[] finish_thread;
}

///////////////////////////////////////////////////////////
inline cx_fmat cx_fmatmul(cx_fmat & mat1, cx_fmat & mat2)
{
    int nz,nx;
    nz=mat1.n_rows;
    nx=mat1.n_cols;

    cx_fmat a(1,1,fill::zeros);
    int i,j;
    for(i=0;i<nz;i++){
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
/////////////beamformingCG3d: !/////////////
int getfftnum(int num){
    int n(1);
    while(n<num){
        n*=2;
    }
    return n;
}
cx_fmat beamforminginv3d_CG_hessianget_thread(struct linerradon3d * par,\
 cx_fmat &hessiancg_cxfmat_p1ncpu_2npx2npy_small, int pncpu, int pnf)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,i,j;//cout<<"ok"<<endl;
    float w,pi(3.1415926);
    float df(par[0].df),dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy),nfft,k;
    nfft=getfftnum(2*npx-1);
    cx_fmat a(1,1),A(nx,ny),B(nx,ny);

    kf=pnf;
    if(kf<par[0].nf2){
        w=2.0*pi*df*(kf); 
        cx_fcube basey3d;
        cx_fmat basex,basey;
        cx_fmat basepy0,basepx0;
        cx_fmat basedpy,basedpx;
        float fpx1,fpy1,fpx2,fpy2,fpxmax,fpxmin,fpymax,fpymin;
        int nkpx,nkpy;
        fpxmax=par[0].px_coord.max()-par[0].px_coord.min();
        fpxmin=par[0].px_coord.min()-par[0].px_coord.max();
        fpymax=par[0].py_coord.max()-par[0].py_coord.min();
        fpymin=par[0].py_coord.min()-par[0].py_coord.max();

        basey.zeros(nx,ny);basex.zeros(nx,ny);
        basey3d.set_size(nx,ny,2*npy-1);
        basepy0.zeros(nx,ny);basepx0.zeros(nx,ny);
        basedpy.zeros(nx,ny);basedpx.zeros(nx,ny);
        basepx0.set_imag(w*par->coordx3d*(fpxmin-dpx));  
        basepy0.set_imag(w*par->coordy3d*(fpymin));
        basedpx.set_imag(w*par->coordx3d*dpx);  
        basedpy.set_imag(w*par->coordy3d*dpy);
        basedpy=exp(basedpy);basedpx=exp(basedpx);
        basepy0=exp(basepy0);basepx0=exp(basepx0);

        basex=basepx0;
        basey3d.slice(0)=basepy0;
        for(kpy2=0;kpy2<(2*npy-2);kpy2++){
        for(kx=0;kx<nx;kx++){
            for(ky=0;ky<ny;ky++){
            basey3d(kx,ky,kpy2+1)=basey3d(kx,ky,kpy2)*basedpy(kx,ky);
        }}}
        for(kpx2=0;kpx2<=(2*npx-2);kpx2++){
            for(kx=0;kx<nx;kx++){
            for(ky=0;ky<ny;ky++){
                basex(kx,ky)=basex(kx,ky)*basedpx(kx,ky);
            }}
        for(kpy2=0;kpy2<=(2*npy-2);kpy2++){
            a=cx_fmatmul(basex,basey=basey3d.slice(kpy2));
            hessiancg_cxfmat_p1ncpu_2npx2npy_small(kpx2,kpy2)\
                =a(0,0);
        }}
    }
    cx_fmat hessfft(nfft, 2*npy-1);
    for(kpy=0;kpy<npy;kpy++){
        float fpy2=(par[0].py_coord(kpy,0));
    for(kpy2=0;kpy2<npy;kpy2++){
    float fpy1=(par[0].py_coord(kpy2,0));
    int nkpy=round((fpy1-fpy2)/dpy)+npy-1;
    cx_fmat a1(nfft,1,fill::zeros);
    kpx=0;
    float fpx1=(par[0].px_coord(kpx,0));
    for(kpx2=0;kpx2<npx;kpx2++){
        float fpx2=(par[0].px_coord(kpx2,0));
        int nkpx=round((fpx1-fpx2)/dpx)+npx-1;
        a1(kpx2,0)=hessiancg_cxfmat_p1ncpu_2npx2npy_small(nkpx,nkpy);
    }
    for(k=1;k<npx;k++){
        a1(nfft-k,0)=a1(k,0);
    }

    kpx2=0;
    float fpx2=(par[0].px_coord(kpx2,0));
    for(kpx=0;kpx<npx;kpx++){
        float fpx1=(par[0].px_coord(kpx,0));
        int nkpx=round((fpx1-fpx2)/dpx)+npx-1;
        a1(kpx,0)=hessiancg_cxfmat_p1ncpu_2npx2npy_small(nkpx,nkpy);
    }
    hessfft.col(nkpy)=fft(a1.col(0));
    }}
    return hessfft;
}

void cxfmatget_Ap_small_fft(cx_fmat & Ap, cx_fmat & hessfft,cx_fmat & p,int pncpu,\
 struct linerradon3d * par)
{
    Ap.fill(0.0);
    int kf,kpx,kpy,kpx2,kpy2,kx,ky,k,nfft;//cout<<"ok"<<endl;
    float df(par[0].df),dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy);
    nfft=hessfft.n_rows;

for(kpy=0;kpy<npy;kpy++){
    float fpy2=(par[0].py_coord(kpy,0));
    cx_fmat dfft(nfft,1,fill::zeros),d(nfft,1,fill::zeros);
    for(kpx=0;kpx<npx;kpx++){
        d(kpx,0)=p(kpx,kpy);
    }
    dfft.col(0)=fft(d.col(0));

for(kpy2=0;kpy2<npy;kpy2++){
    float fpy1=(par[0].py_coord(kpy2,0));
    int nkpy=round((fpy1-fpy2)/dpy)+npy-1;
    cx_fmat a1fft(nfft,1,fill::zeros),a1ifft(nfft,1,fill::zeros);

    a1fft.col(0)=hessfft.col(nkpy);
    for(k=0;k<nfft;k++){
        a1fft(k,0)*=dfft(k,0);
    }
    a1ifft.col(0)=ifft(a1fft.col(0));
    for(kpx=0;kpx<npx;kpx++){
        Ap(kpx,kpy2)+=a1ifft(kpx,0);
    }
}}

    for(kpx=0;kpx<npx;kpx++){
    for(kpy=0;kpy<npy;kpy++){
        Ap(kpx,kpy)+=par->p_power(kpx,kpy)*p(kpx,kpy);
    }}
}
void cxfmatget_Ap_small(cx_fmat & Ap, cx_fmat & A,cx_fmat & p,int pncpu,\
 struct linerradon3d * par)
{
    int kf,kpx,kpy,kpx2,kpy2,kx,ky;//cout<<"ok"<<endl;
    float df(par[0].df),dpx(par[0].dpx),dpy(par[0].dpy),p0x(par[0].p0x);
    int nx(par[0].nx),npx(par[0].npx),nf(par[0].nf),\
        ny(par[0].ny),npy(par[0].npy);

    for(kpx=0;kpx<npx;kpx++){
    for(kpy=0;kpy<npy;kpy++){
        float fpx1=(par[0].px_coord(kpx,0));
        float fpy1=(par[0].py_coord(kpy,0));
        cx_float a;
    for(kpx2=0;kpx2<npx;kpx2++){
    for(kpy2=0;kpy2<npy;kpy2++){
        float fpx2=(par[0].px_coord(kpx2,0));
        float fpy2=(par[0].py_coord(kpy2,0));
        int nkpx=round((fpx1-fpx2)/dpx)+npx-1;
        int nkpy=round((fpy1-fpy2)/dpy)+npy-1;
        a+=A(nkpx,nkpy)*p(kpx2,kpy2);
        if(kpx==kpx2 && kpy==kpy2){
             a+=par->p_power(kpx,kpy)*p(kpx2,kpy2);
        }
    }}
    Ap(kpx,kpy)=a;
    }}
}
inline cx_fmat cx_fmatmul_CG(cx_fmat & mat1, cx_fmat & mat2)
{
    int nz,nx;
    nz=mat1.n_rows;
    nx=mat1.n_cols;

    cx_fmat a(1,1,fill::zeros);
    int i,j;
    for(i=0;i<nz;i++){
        a+=mat1.row(i)*mat2.row(i).t();
    }
    return a;
}
void beamformingCG3d_fthread(struct linerradon3d * par,\
 bool *convergence,int kcpu, int kf,bool regularization,\
 int iterations_num, float residual_ratio,bool *finish_thread)
{
    cx_fmat hess_A(par[0].npx*2-1,par[0].npy*2-1),hessfft;
    hessfft=beamforminginv3d_CG_hessianget_thread(par,hess_A,kcpu,kf);

    int ip,jp,in,jn,k,iter(0);
    int np1(par->datafP.n_rows),np2(par->datafP.n_cols);
    cx_fmat gradient_rk,gradient_rk_1,\
        gradient_cg_pk,gradient_cg_pk_1,\
        datatp_k,datatp_k_1,\
        recoverdatatx_uk,A_gradient_cg_pk,A_datatp_k;
    fmat sum_num(1,1),residual_pow(1,1),residual_k(1,1);
    cx_fmat beta_k(1,1),alpha_k(1,1);
    datatp_k.zeros(np1,np2);
    datatp_k_1.copy_size(datatp_k);
    gradient_rk.copy_size(datatp_k);
    gradient_rk_1.copy_size(datatp_k);
    gradient_cg_pk.copy_size(datatp_k);
    gradient_cg_pk_1.copy_size(datatp_k); 
    A_gradient_cg_pk.copy_size(datatp_k);
    A_datatp_k.copy_size(datatp_k);

    iter=0;
    datatp_k=par[0].datafP.slice(kf);
    //datatp_k.fill(0.0);
    datatp_k.set_real(real(datatp_k)/par[0].nx/par[0].ny);
    datatp_k.set_imag(imag(datatp_k)/par[0].nx/par[0].ny);
    sum_num=sum(sum(sum(abs(datatp_k))));
    residual_pow=sum_num(0,0);
    residual_pow*=residual_ratio;

    //cxfmatget_Ap_small(A_datatp_k,hess_A,datatp_k,kcpu,par);
    cxfmatget_Ap_small_fft(A_datatp_k,hessfft,datatp_k,kcpu,par);
    gradient_rk=par[0].datafP.slice(kf)-A_datatp_k;
    gradient_cg_pk=gradient_rk;

//    cxfmatget_Ap_small(A_gradient_cg_pk,hess_A,gradient_cg_pk,kcpu,par);
    cxfmatget_Ap_small_fft(A_gradient_cg_pk,hessfft,gradient_cg_pk,kcpu,par);
        
    alpha_k=cx_fmatmul_CG(gradient_rk,gradient_rk);
    alpha_k=alpha_k/cx_fmatmul_CG(gradient_cg_pk,A_gradient_cg_pk);
    //get new solution
    datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
    gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;

    beta_k=cx_fmatmul_CG(gradient_rk_1,gradient_rk_1);
    beta_k=beta_k/cx_fmatmul_CG(gradient_rk,gradient_rk);
    gradient_cg_pk_1=gradient_rk_1+beta_k(0,0)*gradient_cg_pk;
    //updata
    datatp_k=datatp_k_1;
    gradient_rk=gradient_rk_1;
    gradient_cg_pk=gradient_cg_pk_1;
    //cal residual_pow
    sum_num=sum(sum(sum(abs(gradient_rk))));
    residual_k=sum_num;

    while(iter<iterations_num && residual_k(0,0)>residual_pow(0,0)){
        iter++;

//        cxfmatget_Ap_small(A_gradient_cg_pk,hess_A,gradient_cg_pk,kcpu,par);
        cxfmatget_Ap_small_fft(A_gradient_cg_pk,hessfft,gradient_cg_pk,kcpu,par);
        
        alpha_k=cx_fmatmul_CG(gradient_rk,gradient_rk);
        alpha_k=alpha_k/cx_fmatmul_CG(gradient_cg_pk,A_gradient_cg_pk);
        //get new solution
        datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
        gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;

        beta_k=cx_fmatmul_CG(gradient_rk_1,gradient_rk_1);
        beta_k=beta_k/cx_fmatmul_CG(gradient_rk,gradient_rk);
        gradient_cg_pk_1=gradient_rk_1+real(beta_k(0,0))*gradient_cg_pk;
        //updata
        datatp_k=datatp_k_1;
        gradient_rk=gradient_rk_1;
        gradient_cg_pk=gradient_cg_pk_1;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(gradient_rk))));
        residual_k=sum_num;
    }
    //cout<<kf<<"||"<<iter<<"||"<<residual_k/residual_pow<<endl;

    if((residual_k(0,0)/residual_pow(0,0))>(0.5/residual_ratio)\
        ||isnan(residual_k(0,0))){
            par[0].datafP.slice(kf).fill(0.0);
            convergence[0]=false;
        }
        else{
            par[0].datafP.slice(kf)=datatp_k;
            convergence[0]=true;
        }
    finish_thread[0]=true;
}

void beamformingCG3d_redo_fthread(struct linerradon3d * par,\
 bool *convergence,int kcpu, int kf,bool regularization,\
 int iterations_num, float residual_ratio,bool *finish_thread)
{
    if(convergence[0]){
        finish_thread[0]=true;
    }
    else{
    cx_fmat hess_A(par[0].npx*2-1,par[0].npy*2-1),hessfft;
    hessfft=beamforminginv3d_CG_hessianget_thread(par,hess_A,kcpu,kf);

    int ip,jp,in,jn,k,iter(0);
    int np1(par->datafP.n_rows),np2(par->datafP.n_cols);
    cx_fmat gradient_rk,gradient_rk_1,\
        gradient_cg_pk,gradient_cg_pk_1,\
        datatp_k,datatp_k_1,\
        recoverdatatx_uk,A_gradient_cg_pk,A_datatp_k;
    fmat sum_num(1,1),residual_pow(1,1),residual_k(1,1);
    cx_fmat beta_k(1,1),alpha_k(1,1);
    datatp_k.zeros(np1,np2);
    datatp_k_1.copy_size(datatp_k);
    gradient_rk.copy_size(datatp_k);
    gradient_rk_1.copy_size(datatp_k);
    gradient_cg_pk.copy_size(datatp_k);
    gradient_cg_pk_1.copy_size(datatp_k); 
    A_gradient_cg_pk.copy_size(datatp_k);
    A_datatp_k.copy_size(datatp_k);

    iter=0;
    datatp_k=par[0].datafP.slice(kf);
    sum_num=sum(sum(sum(abs(datatp_k))));
    residual_pow=sum_num(0,0);
    residual_pow*=residual_ratio;
    //datatp_k.fill(0.0);

//    cxfmatget_Ap_small(A_datatp_k,hess_A,datatp_k,kcpu,par);
    cxfmatget_Ap_small_fft(A_datatp_k,hessfft,datatp_k,kcpu,par);
    gradient_rk=par[0].datafP.slice(kf)-A_datatp_k;
    gradient_cg_pk=gradient_rk;

//    cxfmatget_Ap_small(A_gradient_cg_pk,hess_A,gradient_cg_pk,kcpu,par);
    cxfmatget_Ap_small_fft(A_gradient_cg_pk,hessfft,gradient_cg_pk,kcpu,par);
        
    alpha_k=cx_fmatmul_CG(gradient_rk,gradient_rk);
    alpha_k=alpha_k/cx_fmatmul_CG(gradient_cg_pk,A_gradient_cg_pk);
    //get new solution
    datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
    gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;

    beta_k=cx_fmatmul_CG(gradient_rk_1,gradient_rk_1);
    beta_k=beta_k/cx_fmatmul_CG(gradient_rk,gradient_rk);
    gradient_cg_pk_1=gradient_rk_1+beta_k(0,0)*gradient_cg_pk;
    //updata
    datatp_k=datatp_k_1;
    gradient_rk=gradient_rk_1;
    gradient_cg_pk=gradient_cg_pk_1;
    //cal residual_pow
    sum_num=sum(sum(sum(abs(gradient_rk))));
    residual_k=sum_num;

    while(iter<iterations_num && residual_k(0,0)>residual_pow(0,0)){
        //cout<<iter<<"||"<<residual_k/residual_pow<<endl;
        iter++;

//        cxfmatget_Ap_small(A_gradient_cg_pk,hess_A,gradient_cg_pk,kcpu,par);
        cxfmatget_Ap_small_fft(A_gradient_cg_pk,hessfft,gradient_cg_pk,kcpu,par);
        
        alpha_k=cx_fmatmul_CG(gradient_rk,gradient_rk);
        alpha_k=alpha_k/cx_fmatmul_CG(gradient_cg_pk,A_gradient_cg_pk);
        //get new solution
        datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
        gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;

        beta_k=cx_fmatmul_CG(gradient_rk_1,gradient_rk_1);
        beta_k=beta_k/cx_fmatmul_CG(gradient_rk,gradient_rk);
        gradient_cg_pk_1=gradient_rk_1+real(beta_k(0,0))*gradient_cg_pk;
        //updata
        datatp_k=datatp_k_1;
        gradient_rk=gradient_rk_1;
        gradient_cg_pk=gradient_cg_pk_1;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(gradient_rk))));
        residual_k=sum_num;
    }
    if((residual_k(0,0)/residual_pow(0,0))>(0.5/residual_ratio)\
        ||isnan(residual_k(0,0))){
            par[0].datafP.slice(kf).fill(0.0);
            convergence[0]=false;
            cout<<"kf="<<kf<<" ||iteration times:"<<iter<<" ||err level:"<<\
            residual_k(0,0)/residual_pow(0,0)<<" ("<<convergence[0]<<endl;
        }
        else{
            par[0].datafP.slice(kf)=datatp_k;
            convergence[0]=true;
        }
    finish_thread[0]=true;
    }
}

void beamformingCG3d(struct linerradon3d & par,\
int iterations_num=45, float residual_ratio=0.1,\
bool regularization=true)
{
    int ncpu(par.numthread),numf(par.nf2-par.nf1-ncpu);
    regularization=par.regularization;
    int kcpu,kf,k,i,j;//cout<<"ok"<<endl;

    struct linerradon3d * ppar;
    ppar=&par;
    fmat dig_w_fmat(par.npx,par.npy);
    beamforminginv3d_getdigfmat(ppar[0], dig_w_fmat);
    datawrite(dig_w_fmat,par.npx,par.npy,"dig.bin");
    
    bool *finish_thread,*convergence;
    convergence=new bool[par.nf2];
    finish_thread=new bool[ncpu];
    thread *pcal;
    pcal=new thread[ncpu];
    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(beamformingCG3d_fthread,&par,\
            &(convergence[kcpu+par.nf1]),kcpu,kcpu+par.nf1,regularization,\
            iterations_num,residual_ratio,&(finish_thread[kcpu]));
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable()&&kf<par.nf2&&finish_thread[kcpu]){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(beamformingCG3d_fthread,&par,\
                    &(convergence[kf]),kcpu,kf,regularization,\
                    iterations_num,residual_ratio,&(finish_thread[kcpu]));
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }
cout<<"==========cg has finished=========="<<endl;
    cx_fmat s1,s2;
    bool *convergence2;
    convergence2=new bool[par.nf2];
    int kk,k1(0),k2(0);

for(kk=0;kk<iterations_num;kk++){
    k2=0;
    for(kf=par.nf1;kf<par.nf2;kf++){
        convergence2[kf]=convergence[kf];
        if(!convergence[kf]){
            k2++;
        }
    }
    if(k2==k1){
        break;
    }
    else{
        k1=k2;
    }
    for(kf=par.nf2-2;kf>=par.nf1+1;kf--){
        if(!convergence2[kf]&&convergence2[kf-1]&&convergence2[kf+1]){
            s1=par.datafP.slice(kf-1);
            s1+=par.datafP.slice(kf+1);
            par.datafP.slice(kf).set_real(real(s1)/2.0);
            par.datafP.slice(kf).set_imag(imag(s1)/2.0);
            convergence2[kf]=true;
        }
        else if(!convergence2[kf]&&convergence2[kf+1]){
            par.datafP.slice(kf)=par.datafP.slice(kf+1);
            convergence2[kf]=true;
        }
    }

    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(beamformingCG3d_redo_fthread,&par,\
            &(convergence[kcpu+par.nf1]),kcpu,kcpu+par.nf1,regularization,\
            iterations_num,residual_ratio,&(finish_thread[kcpu]));
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable()&&kf<par.nf2&&finish_thread[kcpu]){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(beamformingCG3d_redo_fthread,&par,\
                    &(convergence[kf]),kcpu,kf,regularization,\
                    iterations_num,residual_ratio,&(finish_thread[kcpu]));
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }
    cout<<"==========cg redo "<<kk<<" finished=========="<<endl;
}

    for(kf=par.nf1+1;kf<par.nf2-1;kf++){
        if(!convergence[kf]&&convergence[kf-1]&&convergence[kf+1]){
            s1=par.datafP.slice(kf-1);
            s1+=par.datafP.slice(kf+1);
            par.datafP.slice(kf).set_real(real(s1)/2.0);
            par.datafP.slice(kf).set_imag(imag(s1)/2.0);
            convergence[kf]=true;
        }
        else if(!convergence[kf]&&convergence[kf-1]){
            par.datafP.slice(kf)=par.datafP.slice(kf-1);
            convergence[kf]=true;
        }
    }

    fx2tx_3d_thread(par.realdataTP,par.datafP,par.numthread);
    delete[] pcal;
    delete[] finish_thread;
    delete[] convergence;
    delete[] convergence2;
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
// 2-D Linear Radon Transform in the frequence - space domain.
int Beamforming_LS_2D(float **trace, int ntrace, int ns, float dt,\
 float *coor, float x0, float **tauppanel, int npsr,float psrmin, float dpsr,\
 float fmax=150,float frule=50,int ncpu=1, bool regularization=false,\
 float factor_L2=0.1, int iterations_num=45, float residual_ratio=0.1)
{
    struct linerradon3d par;

    int i,j,k;
    int nz2(ns),nz(ns),nx(ntrace),ny(1),nf(ns);
//////////////////////////radon par-set////////////////////////////
    beamforming_parset(nx,ny,nz,par);
    par.dpx=dpsr;
    par.dpy=0;
    par.dz=dt;

//The default px of central channel is zero
    par.npx=npsr;
    par.p0x=psrmin;
    par.npy=1;
    par.p0y=0.0;   

//ncpu
    par.numthread=ncpu;

//regularization parameter
    par.dig_n2=nx*ny*factor_L2;  //L2, Tikhonov 
//Seismic trace coordinates
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            par.coordx3d(i,j)=coor[i]-x0;
            par.coordy3d(i,j)=0;
        }}
    par.regularization=regularization;
//Parameters updated
    beamforming_parupdate(par);
//Frequency calculation range (number)
    par.nf2=int(fmax/par.df); 
    par.nf1=1; 
//Low frequency constraint range (number)
    par.rulef2=int(frule/par.df);
    par.rulef1=1;

//////////////////////////////////////////////////////////////////
//data input, (row,col,slice) of par.data is (nx,ny,nz)
//nx: Spatial sampling nummber
//ny: trace nummber
//nz: time sampling nummber
for(i=0;i<ntrace;i++){
    for(j=0;j<ns;j++){
        par.data(i,0,j)=trace[i][j];
    }
}
//////////////////////////////////////////////////
    linerradon(par); 
    //inv radon
    //beamforming_cleardata(par);
    //beamforminginv3d(par,0);
    beamformingCG3d(par,iterations_num,residual_ratio);
    //recover data
    rebuildsignal(par);

///////////////////output tau-p///////////////////
//output tau-p, (row,col,slice) of par.realdataTP is (npx,npy,nz)
//npx: X Ray parameters sampling nummber
//npy: Y Ray parameters sampling nummber
//nz: time sampling nummber
for(i=0;i<npsr;i++){
    for(j=0;j<ns;j++){
        tauppanel[i][j]=par.realdataTP(i,0,j);
    }
}
////////////////////output recover data/////////////////////
//output recover data, (row,col,slice) of par.realrebuildtx is (nx,ny,nz)
//npx: X Ray parameters sampling nummber
//npy: Y Ray parameters sampling nummber
//nz: time sampling nummber

    return 0;
}
int Beamforming_CG_2D(fmat &tauppanel,fmat &recoverdata,fmat &recovererr,\
 fmat trace, fmat coor, int ns, int ntrace, float dt,int npsr,float psrmin, float dpsr,\
 float fmax=150,float frule=50,int ncpu=1, float factor_L2=0.1,float factor_L1=1,\
 int iterations_num=45, float residual_ratio=0.1, bool dorecover=false)
{
    struct linerradon3d par;
    int i,j,k;
    int nz(ns),nx(ntrace),ny(1),nf(ns);
//////////////////////////radon par-set////////////////////////////
    beamforming_parset(nx,ny,nz,par);
    par.dpx=dpsr;
    par.dpy=dpsr;
    par.dz=dt;

//The default px of central channel is zero
    par.npx=npsr;
    par.p0x=psrmin;
    par.npy=1;
    par.p0y=0.0;  

//ncpu
    par.numthread=ncpu;

//regularization parameter
    par.dig_n2=nx*ny*factor_L2;  //L2, Tikhonov 
    par.dig_n1=nx*ny*factor_L1;  //L1,
//Seismic trace coordinates
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            par.coordx3d(i,j)=coor(i,j);
            par.coordy3d(i,j)=0;
        }}
    par.regularization=false;
//Parameters updated
    beamforming_parupdate(par);
//Frequency calculation range (number)
    par.nf2=int(fmax/par.df); 
    par.nf1=1; 
//Low frequency constraint range (number)
    par.rulef2=int(frule/par.df);
    par.rulef1=1;

    par.data.col(0)=trace;
//////////////////////////////////////////////////
    linerradon(par); 

    beamformingCG3d(par,iterations_num,residual_ratio);
    //beamforminginv3d(par);
    tauppanel=par.realdataTP.col(0);

    if(dorecover){
        rebuildsignal(par);

        recoverdata=par.realrebuildtx.col(0);
        recovererr=recoverdata-trace;
        //recover data
    }

    return 0;
}
int Beamforming_recoverdata_2D(fmat &recoverdata,fmat &tauppanel,\
 fmat coor, int ns, int ntrace, float dt,int npsr,float psrmin,\
 float dpsr,float fmax=150,float frule=50,int ncpu=1)
{
    struct linerradon3d par;
    int i,j,k;
    int nz(ns),nx(ntrace),ny(1),nf(ns);
//////////////////////////radon par-set////////////////////////////
    beamforming_parset(nx,ny,nz,par);
    par.dpx=dpsr;
    par.dpy=dpsr;
    par.dz=dt;

//The default px of central channel is zero
    par.npx=npsr;
    par.p0x=psrmin;
    par.npy=1;
    par.p0y=0.0;  

//ncpu
    par.numthread=ncpu;

//regularization parameter

//Seismic trace coordinates
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            par.coordx3d(i,j)=coor(i,j);
            par.coordy3d(i,j)=0;
        }}
    par.regularization=false;
//Parameters updated
    beamforming_parupdate(par);
//Frequency calculation range (number)
    par.nf2=int(fmax/par.df); 
    par.nf1=1; 
//Low frequency constraint range (number)
    par.rulef2=int(frule/par.df);
    par.rulef1=1;

//////////////////////////////////////////////////
    par.realdataTP.col(0)=tauppanel;

    rebuildsignal(par);

    recoverdata=par.realrebuildtx.col(0);
    //recover data

    return 0;
}

int Beamforming_CG_3D(fcube &tauppanel,fcube &recoverdata,fcube &recovererr,\
 fcube& trace, fmat coordx,fmat coordy, int ns, int ntrace,int nline,float dt,\
 int npx,float pxmin, float dpx,int npy,float pymin, float dpy,\
 float fmax=150,float frule=50,int ncpu=1, float factor_L2=0.1,float factor_L1=1,\
 int iterations_num=45, float residual_ratio=0.1, bool dorecover=false)
{
    struct linerradon3d par;
    int i,j,k;
    int nz(ns),nx(ntrace),ny(nline),nf(ns);
//////////////////////////radon par-set////////////////////////////
    beamforming_parset(nx,ny,nz,par);
    par.dpx=dpx;
    par.dpy=dpy;
    par.dz=dt;

//The default px of central channel is zero
    par.npx=npx;
    par.p0x=pxmin;
    par.npy=npy;
    par.p0y=pymin;  

//ncpu
    par.numthread=ncpu;

//regularization parameter
    par.dig_n2=nx*ny*factor_L2;  //L2, Tikhonov 
    par.dig_n1=nx*ny*factor_L1;  //L1,
//Seismic trace coordinates
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            par.coordx3d(i,j)=coordx(i,j);
            par.coordy3d(i,j)=coordy(i,j);
        }}
    par.regularization=false;
//Parameters updated
    beamforming_parupdate(par);
//Frequency calculation range (number)
    par.nf2=int(fmax/par.df); 
    par.nf1=1; 
//Low frequency constraint range (number)
    par.rulef2=int(frule/par.df);
    par.rulef1=1;

    par.data=trace;
    trace.set_size(1,1,1);
    tauppanel.set_size(1,1,1);
//////////////////////////////////////////////////
    linerradon(par); 

    beamformingCG3d(par,iterations_num,residual_ratio);
    //beamforminginv3d(par);
    tauppanel=par.realdataTP;
    trace=par.data;

    if(dorecover){
        //recover data
        rebuildsignal(par);
        recoverdata=par.realrebuildtx;
        recovererr=recoverdata-trace;
    }
    
    return 0;
}

int Beamforming_recoverdata_3D(fcube &recoverdata,fcube &tauppanel,\
 fmat coordx,fmat coordy, int ns, int ntrace,int nline,float dt,\
 int npx,float pxmin, float dpx,int npy,float pymin, float dpy,\
 float fmax=150,float frule=50,int ncpu=1)
{
    struct linerradon3d par;
    int i,j,k;
    int nz(ns),nx(ntrace),ny(nline),nf(ns);
//////////////////////////radon par-set////////////////////////////
    beamforming_parset(nx,ny,nz,par);
    par.dpx=dpx;
    par.dpy=dpy;
    par.dz=dt;

//The default px of central channel is zero
    par.npx=npx;
    par.p0x=pxmin;
    par.npy=npy;
    par.p0y=pymin;  

//ncpu
    par.numthread=ncpu;

//Seismic trace coordinates
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            par.coordx3d(i,j)=coordx(i,j);
            par.coordy3d(i,j)=coordy(i,j);
        }}
    par.regularization=false;
//Parameters updated
    beamforming_parupdate(par);
//Frequency calculation range (number)
    par.nf2=int(fmax/par.df); 
    par.nf1=1; 
//Low frequency constraint range (number)
    par.rulef2=int(frule/par.df);
    par.rulef1=1;

//////////////////////////////////////////////////

    par.realdataTP=tauppanel;
    //recover data
    rebuildsignal(par);
    recoverdata=par.realrebuildtx;

    
    return 0;
}

#endif
