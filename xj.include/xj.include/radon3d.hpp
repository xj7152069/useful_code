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
 float *coor, float x0, float **tauppanel, int npsr,float psrmin, float dpsr);
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
/*****************************************************/
    fcube data,realdataTP,realrebuildtx;
    cx_fcube datafx,rebuildfx,datafP,rebuildfp,dataTP,rebuildtx;
/*****************************************************/
    fmat p_power,py_coord,px_coord,x_coord,y_coord;
    char allAreal[99],allAimag[99];
    float dig_n1,dig_n2,kerpar1;
    float x_center,y_center;
    bool lsinvmat;
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
    //par.p0x=-par.dpx*int(par.npx/2)+par.px_center;
    //par.p0y=-par.dpy*int(par.npy/2)+par.py_center;
    par.nf1=0;par.nf2=par.nf/2;
    par.rulef1=0;par.rulef2=par.nf/2;
    par.numthread=1;
    par.dig_n=(1);
    par.dig_n1=(1);
    par.dig_n2=(1);
    par.allAreal[0]='\0';
    strcat(par.allAreal,"./allAreal.bin");
    par.allAimag[0]='\0';
    strcat(par.allAimag,"./allAimag.bin");
    par.lsinvmat=false;
    par.dig_n=par.dig_n*par.nx*par.ny;
    par.dig_n1=par.dig_n1*par.nx*par.ny;
    par.dig_n2=par.dig_n2*par.nx*par.ny;
    par.x_coord.zeros(par.nx,1);
    par.y_coord.zeros(par.ny,1);
    int k;
    for(k=0;k<par.nx;k++)
    {
        par.x_coord(k,0)=k*par.dx-(int(par.nx/2)*par.dx)+par.x_center;
    }
    for(k=0;k<par.ny;k++)
    {
        par.y_coord(k,0)=k*par.dy-(int(par.ny/2)*par.dy)+par.y_center;
    }
    par.dig_n1=0.0;  //L1
    par.kerpar1=0.0;     //kernel function par
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
    //cout<<"nf = "<<nf<<endl;
    beamforming_parset(nx,ny,nz,nf,ppar[0]);
}

void beamforming_parupdate(struct linerradon3d & par)
{ 
    par.df=1.0/par.dz/par.nf;
    //par.p0x=-par.dpx*int(par.npx/2)+par.px_center;
    //par.p0y=-par.dpy*int(par.npy/2)+par.py_center;
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
    par.nf1=1;
    par.rulef1=1;
    par.nf2=int(par.nf2/par.df);
    par.rulef2= par.nf2;
    if(par.nf2>par.nf/2)
    {
        par.nf2=par.nf/2;
    }
    if(par.rulef2>par.nf2/2)
    {
        par.rulef2=par.nf2/2;
    }
    par.allAreal[0]='\0';
    strcat(par.allAreal,"./allAreal.bin");
    par.allAimag[0]='\0';
    strcat(par.allAimag,"./allAimag.bin");
    par.lsinvmat=false;

    par.px_coord.zeros(par.npx,1);
    par.py_coord.zeros(par.npy,1);
    
    int k;
    for(k=0;k<par.npx;k++)
    {
        par.px_coord(k,0)=k*par.dpx+par.p0x;
    }
    for(k=0;k<par.npy;k++)
    {
        par.py_coord(k,0)=k*par.dpy+par.p0y;
    }

    /*
    ofstream outf1;
    outf1.open(par.allAreal);
    outf1.close();
    outf1.open(par.allAimag);
    outf1.close();*/
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
    for(k=par.rulef1;k<par.rulef2;k++){
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
                par.p_power(i,k)=1.0/(1.0+exp(10*(0.5-par.p_power(i,k))));
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
            //datawrite(digline=par.p_power.col(k),par.npx,1,"dig.bin");
        }
    }

    //get diagonally weighted matrix W
    parinv.digw_fmat_npxnpy=par.dig_n2+par.dig_n1*par.p_power;
    par.p_power=parinv.digw_fmat_npxnpy;
}

//get mat LS-LTL
void beamforminginv3d_LS_hessianget_thread(struct linerradon3d * par,\
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
    ifstream inf1,inf2;
    ofstream outf1,outf2;

    float ky1(-(int(ny/2)*dy));
    float kx1(-(int(nx/2)*dx));
    cx_fmat a(1,1),A(nx,ny),B(nx,ny);
    fmat hess(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows,\
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols);

    kf=pnf;
    if(kf<par[0].nf2) {
        if(!par[0].lsinvmat){
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

            outf1.open(par[0].allAreal,ios::app);
            outf1.seekp((kf-par[0].nf1)\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols\
                *sizeof(float),ios::beg);
            datawrite(hess=real(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]),\
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows,\
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols,outf1);
            outf1.close();

            outf1.open(par[0].allAimag,ios::app);
            outf1.seekp((kf-par[0].nf1)\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols\
                *sizeof(float),ios::beg);
            datawrite(hess=imag(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu]),\
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows,\
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols,outf1);
            outf1.close();

        }
        else if(par[0].lsinvmat){
            inf1.open(par[0].allAreal);
            inf1.seekg((kf-par[0].nf1)\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols\
                *sizeof(float),ios::beg);
            hess=dataread(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows,\
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols,inf1);
            inf1.close();
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].set_real(hess);

            inf1.open(par[0].allAimag);
            inf1.seekg((kf-par[0].nf1)\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows\
                *par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols\
                *sizeof(float),ios::beg);
            hess=dataread(par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_rows,\
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].n_cols,inf1);
            inf1.close();
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu].set_imag(hess);
        }
    }
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
    cx_fmat basey(ny,npy,fill::zeros),basex(nx,npx,fill::zeros); //

    kf=pnf;
    if(kf<par[0].nf2) {
        w=2.0*pi*df*(kf); 
        if(!regularization){
        basey.zeros(ny,npy);basex.zeros(nx,npx); //
        basey.set_imag(w*(par->y_coord*par->py_coord.t())); //
        basex.set_imag(w*(par->x_coord*par->px_coord.t())); //
        basey=exp(basey);basex=exp(basex); //
        }
        for(kpx=0;kpx<npx;kpx++){
        for(kpy=0;kpy<npy;kpy++){
            i=kpx*npy+kpy;
            if(!regularization){
            A=basex.col(kpx)*basey.col(kpy).st();
            }
            float fpx1(par[0].px_coord(kpx,0));
            float fpy1(par[0].py_coord(kpy,0));
            for(kpx2=0;kpx2<npx;kpx2++){
            for(kpy2=0;kpy2<npy;kpy2++){
                j=kpx2*npy+kpy2;
                if(!regularization){
                B=basex.col(kpx2)*basey.col(kpy2).st();
                B.set_imag(-imag(B));
                a=cx_fmatmul(A,B);
                }
                else{
                float fpx2(par[0].px_coord(kpx2,0));
                float fpy2(par[0].py_coord(kpy2,0));
                a=cx_hessianelement2d(w,nx,kx1,dx,fpx1-fpx2,\
                    ny,ky1,dy,fpy1-fpy2);
                }
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,j)=a(0,0);
                par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](j,i)=a(0,0); 
            }}        
            par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[pncpu](i,i)+=par->p_power(kpx,kpy);
            kpx2=npx;kpy2=npy;
        }}

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

    struct beamforminginv3d par2;
    struct linerradon3d * ppar;
    ppar=&par;
    beamforminginv3d_parset(ppar[0], par2);
    beamforminginv3d_getdigfmat(ppar[0], par2);

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

    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
}

//LS-inv function of radon
void beamformingLSinv3d(struct linerradon3d & par)
{
    int numf(par.nf2-par.nf1),ncpu(par.numthread);
    int kf,kcpu,i,j;//cout<<"ok"<<endl;

    struct beamforminginv3d par2;
    struct linerradon3d * ppar;
    ppar=&par;
    par.dig_n1=0.0;
    beamforminginv3d_parset(ppar[0], par2);
    if(!par.lsinvmat)
        {beamforminginv3d_getdigfmat(ppar[0], par2);}

    bool* finish_thread;
    finish_thread=new bool[ncpu];
    thread *pcal;
    pcal=new thread[ncpu];
    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(beamforminginv3d_beamget_thread,&par,&par2,\
                kcpu,kcpu+par.nf1,true,&(finish_thread[kcpu]));
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable()&&kf<par.nf2&&(finish_thread[kcpu])){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(beamforminginv3d_beamget_thread,&par,&par2,\
                kcpu,kf,true,&(finish_thread[kcpu]));
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }

    par.lsinvmat=true;
    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
    delete[] pcal;
    delete[] finish_thread;
}

///////////////linerradon 3D transform: slant stack///////////////
void linerradon_fthread(struct linerradon3d * par,int pnf1, int pnf2,\
 bool *finish_thread)
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
    finish_thread[0]=true;
}

void linerradon(struct linerradon3d & par)
{
    tx_to_fx3d_radon3d_thread(par);
    
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
    
    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
    delete[] pcal;
    delete[] finish_thread;
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
void rebuildsignal_fthread(struct linerradon3d * par,int pnf1, int pnf2,\
 bool *finish_thread)
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
    finish_thread[0]=true;
}

void rebuildsignal(struct linerradon3d & par)
{
    tp_to_fp3d_linerradon3d(par);

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
   
    //fx_to_tx3d_linerradon3d(par);
    fx_to_tx3d_radon3d_thread(par);
    par.realrebuildtx=real(par.rebuildtx);
    delete[] pcal;
    delete[] finish_thread;
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
    //datawrite(datapower,par.nf,1,"./fpower.bin");
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
void cxfmatget_Ap(cx_fmat & Ap,cx_fmat & A,cx_fmat & p)
{
    int npx(Ap.n_rows),npy(Ap.n_cols);
    int i,j,kpx,kpx2,kpy,kpy2;
    cx_fmat a(1,1),b(1,1);

    for(kpx=0;kpx<npx;kpx++){
    for(kpy=0;kpy<npy;kpy++){
        i=kpx*npy+kpy;
        a(0,0).real(0.0);
        a(0,0).imag(0.0);
        for(kpx2=0;kpx2<npx;kpx2++){
        for(kpy2=0;kpy2<npy;kpy2++){
            j=kpx2*npy+kpy2;
            a(0,0)+=A(i,j)*p(kpx2,kpy2);
        }}
        Ap(kpx,kpy)=a(0,0);
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
 struct beamforminginv3d * par2,bool *finish_thread,\
 int kcpu, int kf,bool regularization,\
 int iterations_num=25, float residual_ratio=0.5)
{
    beamforminginv3d_hessianget_thread(\
        par,par2,kcpu,kf,regularization,false);
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
    datatp_k.set_real(2*real(datatp_k)/par[0].nx/par[0].ny);
    datatp_k.set_imag(2*imag(datatp_k)/par[0].nx/par[0].ny);

    cxfmatget_Ap(A_datatp_k,par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu],\
        datatp_k);
    gradient_rk=par[0].datafP.slice(kf)-A_datatp_k;
    gradient_cg_pk=gradient_rk;

    //cal residual_pow
    sum_num=sum(sum(sum(abs(gradient_rk))));
    residual_k=sum_num;
    cxfmatget_Ap(A_gradient_cg_pk,par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu],\
        gradient_cg_pk);
        
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

    while(iter<iterations_num && residual_k(0,0)>residual_pow(0,0)){
        //cout<<iter<<"||"<<residual_k/residual_pow<<endl;
        iter++;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(gradient_rk))));
        residual_k=sum_num;

        cxfmatget_Ap(A_gradient_cg_pk,par2[0].hessianinv_cxfmat_p1ncpu_npxynpxy[kcpu],\
            gradient_cg_pk);
        
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
    }
    par[0].datafP.slice(kf)=datatp_k;
    finish_thread[0]=true;
    cout<<"kf="<<kf<<"||iteration times: "<<iter<<"||residual level:"<<\
        residual_k/residual_pow<<endl;
}

void beamformingCG3d(struct linerradon3d & par,\
int iterations_num=45, float residual_ratio=0.1,\
bool regularization=true)
{
    int ncpu(par.numthread),numf(par.nf2-par.nf1-ncpu);
    int kcpu,kf,k,i,j;//cout<<"ok"<<endl;
    float pi(3.1415926);

    struct beamforminginv3d par2;
    struct linerradon3d * ppar;
    ppar=&par;
    beamforminginv3d_parset(ppar[0], par2);
    beamforminginv3d_getdigfmat(ppar[0], par2);
    
    bool* finish_thread;
    finish_thread=new bool[ncpu];
    thread *pcal;
    pcal=new thread[ncpu];
    for(kcpu=0;kcpu<ncpu;kcpu++){
        finish_thread[kcpu]=false;
        pcal[kcpu]=thread(beamformingCG3d_fthread,&par,&par2,\
            &(finish_thread[kcpu]),kcpu,kcpu+par.nf1,regularization,\
            iterations_num,residual_ratio);
    }
    kf=ncpu+par.nf1;
    while(kf<par.nf2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
            if(pcal[kcpu].joinable()&&kf<par.nf2&&finish_thread[kcpu]){
                pcal[kcpu].join();
                finish_thread[kcpu]=false;
                pcal[kcpu]=thread(beamformingCG3d_fthread,&par,&par2,\
                    &(finish_thread[kcpu]),\
                    kcpu,kf,regularization,iterations_num,residual_ratio);
                kf++;
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }

    fp_to_tp3d_linerradon3d(par);
    par.realdataTP=real(par.dataTP);
    delete[] pcal;
    delete[] finish_thread;
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
 float fmax=150, int ncpu=1)
{
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
    struct linerradon3d par;

    int i,j,k;
    int nz2(ns),nz(ns),nx(ntrace),ny(1),nf(0);
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

//Frequency calculation range (number)
    par.nf2=fmax; 
//Low frequency constraint range (number)
    par.rulef2=fmax/2;
//regularization parameter
    par.dig_n2=nx*ny*0.1;  //L2, Tikhonov 
//Seismic trace coordinates
    for(k=0;k<par.nx;k++){
        par.x_coord(k,0)=coor[k]-x0;
    }
//Parameters updated
    beamforming_parupdate(par);

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
    beamforming_cleardata(par);
    beamforminginv3d(par,0);
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

#endif
