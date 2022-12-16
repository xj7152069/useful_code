/************update in 2021.01.07**************/

    
    
/*********************************************/

#ifndef WIENER3D_HPP
#define WIENER3D_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "../xjc.h"
#include <thread>
#include <future>
using namespace std;
using namespace arma;

inline fmat invfmat(fmat mat1, int n);
inline cx_fmat invcxfmat(cx_fmat mat1, int n);
inline cx_fmat solvecxfmat(cx_fmat mata, cx_fmat matb, int n);
inline fmat solvefmat(fmat mata, fmat matb, int n);

struct wiener3d{
    int nz,nx,ny,nf,nf1,nf2,nwx,nwy,narx,nary,\
        halfnarx,halfnary,nmovex,nmovey,ncpu,covermax;
    float dz,dy,dx,df,dig_n,dig_n2,normal_con;
    fcube data,realrebuildtx,dataagc;
    cx_fcube rebuildfx,rebuildtx,datafx;
    cx_fcube fx3dInterpolation;
    fcube tx3dInterpolation;
    int use_way,agcwx,agcwy,agcwz;
    bool use_anti_alias,normal_bool;
};

//设置一组常用的初始化参数
//Mat(nx,ny,nz); Mat(nx,ny,nf)
void wiener3d_parset(int nx, int ny, int nz,\
    struct wiener3d & par)
{
    par.normal_bool=false,par.normal_con=1.0;
    par.use_way=1;
    par.nx=nx,par.ny=ny,par.nz=nz;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.dig_n=0.001;par.dig_n2=0.001;
    par.nwx=1+floor(nx/2.0+0.01);par.nwy=1+floor(ny/2.0+0.01);          
    par.narx=floor(par.nwx/2.0+0.01);par.nary=floor(par.nwy/2.0+0.01);
    if(par.nwx>21){par.nwx=21;}
    if(par.nwy>21){par.nwy=21;}
    if(par.narx>5){par.narx=5;}
    if(par.nary>5){par.nary=5;}
    par.nmovex=par.nwx;par.nmovey=par.nwy;
    par.covermax=1;
    par.ncpu=1;par.nf=1;
    //par.df=1.0/par.dz/par.nf;
    par.halfnarx=int(par.narx/2);
    par.narx=2*par.halfnarx+1;
    par.halfnary=int(par.nary/2);
    par.nary=2*par.halfnary+1;
    par.nf=par.nz;
    /*
    while (par.nf<par.nz){
        par.nf*=2;
    }*/
    cout<<"nf = "<<par.nf<<endl;;
    par.nf1=0,par.nf2=par.nf/2;
    par.df=1.0/par.dz/par.nf;
    if(par.nx<par.ny){
        cout<<"error in wiener3d: find {nx < ny};\
         the data size {nx should >= ny}"<<endl;
    }
    par.use_anti_alias=false;
    //par.datafx.zeros(par.nx,par.ny,par.nf);
    par.data.zeros(par.nx,par.ny,par.nz);
    par.dataagc.zeros(par.nx,par.ny,par.nz);
    par.dataagc.fill(1.0);
    //par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nz);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
}

void wiener3d_parupdate(struct wiener3d & par){
    int k;
    par.nwx=min(1+(k=floor(par.nx/2.0+0.01)),par.nwx);
    par.nwy=min(1+(k=floor(par.ny/2.0+0.01)),par.nwy);          
    par.narx=floor(par.nwx/2.0+0.01);
    par.nary=floor(par.nwy/2.0+0.01);
    if(par.nwx>21){par.nwx=21;}
    if(par.nwy>21){par.nwy=21;}
    if(par.narx>5){par.narx=5;}
    if(par.nary>5){par.nary=5;}
    par.nmovex=par.nwx;par.nmovey=par.nwy;
    par.agcwx=min(par.narx,par.agcwx);
    par.agcwy=min(par.nary,par.agcwy);
    par.agcwz=min(par.agcwz,par.nz/5);

    par.halfnarx=floor(0.01+par.narx/2.0);
    par.narx=2*par.halfnarx+1;
    par.halfnary=floor(0.01+par.nary/2.0);
    par.nary=2*par.halfnary+1;
    //par.agcwx=par.halfnarx;
    //par.agcwy=par.halfnary;
    par.agcwx=2*par.narx;
    par.agcwy=2*par.nary;
    par.agcwz=50;
    //par.data.zeros(par.nx,par.ny,par.nz);
    par.dataagc.zeros(par.nx,par.ny,par.nz);
    par.dataagc.fill(1.0);
    //par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nz);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
    par.datafx.zeros(par.nx,par.ny,par.nf);
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    if(par.ncpu>par.nx){
        par.ncpu=par.nx;
    }
}

void wiener3d_cleardata(struct wiener3d & par){
    par.datafx.zeros(1,1,1);
    par.rebuildfx.zeros(1,1,1);
}

//void getagc_wiener3d_thread(struct wiener3d & par)
void getagc_wiener3d(struct wiener3d & par){
    int i,j,k,i1,j1,k1;
    float agcpow;
    for(i=par.agcwx;i<par.nx-par.agcwx;i++){
    for(j=par.agcwy;j<par.ny-par.agcwy;j++){
    for(k=par.agcwz;k<par.nz-par.agcwz;k++){
        agcpow=0;
        for(i1=-par.agcwx;i1<=par.agcwx;i1++){
        for(j1=-par.agcwy;j1<=par.agcwy;j1++){
        for(k1=-par.agcwz;k1<=par.agcwz;k1++){
            agcpow+=abs(par.data(i+i1,j+j1,k+k1));
        }}}
        par.dataagc(i,j,k)=agcpow;
        agcpow=0;
    }}}
    for(i=0;i<par.agcwx;i++){
        par.dataagc.row(i)=par.dataagc.row(par.agcwx);
        par.dataagc.row(par.nx-1-i)=par.dataagc.row(par.nx-1-par.agcwx);
    }
    for(i=0;i<par.agcwy;i++){
        par.dataagc.col(i)=par.dataagc.col(par.agcwy);
        par.dataagc.col(par.ny-1-i)=par.dataagc.col(par.ny-1-par.agcwy);
    }
    for(i=0;i<par.agcwz;i++){
        par.dataagc.slice(i)=par.dataagc.slice(par.agcwz);
        par.dataagc.slice(par.nz-1-i)=par.dataagc.slice(par.nz-1-par.agcwz);
    }
    par.dataagc=par.dataagc/(par.dataagc.max()-par.dataagc.min());
    cout<<par.dataagc.max()<<"|"<<par.dataagc.min()<<endl;
    par.dataagc=par.dataagc+0.001*(par.dataagc.max()-par.dataagc.min());
    par.dataagc=1.0/par.dataagc;
    for(i=0;i<par.nx;i++){
    for(j=0;j<par.ny;j++){
    for(k=0;k<par.nz;k++){
        par.data(i,j,k)=par.dataagc(i,j,k)*par.data(i,j,k);
    }}}
}
void deagc_wiener3d(struct wiener3d & par)
{
    int i,j,k;
    for(i=0;i<par.realrebuildtx.n_rows;i++){
    for(j=0;j<par.realrebuildtx.n_cols;j++){
    for(k=0;k<par.realrebuildtx.n_slices;k++){
        par.realrebuildtx(i,j,k)=par.realrebuildtx(i,j,k)/par.dataagc(i,j,k);
        //par.data(i,j,k)=par.data(i,j,k)/par.dataagc(i,j,k);
    }}}
}

void tx_to_fx3d_wiener3d(struct wiener3d & par){
    int i,j;
    fmat data(par.ny,par.nz);
    cx_fmat datafx(par.ny,par.nf);
    for(i=0;i<par.nx;i++){
        data=par.data.row(i);
        for(j=0;j<par.ny;j++){
            datafx.row(j)=fft(data.row(j),par.nf);
        }
        par.datafx.row(i)=datafx;
    } 
}

void fx_to_tx3d_wiener3d(struct wiener3d & par){
    int i,j;
    cx_fmat data(par.ny,par.nf),datafx(par.ny,par.nf),data2(par.ny,par.nz);
    for(i=0;i<par.nx;i++){
        datafx=par.rebuildfx.row(i);
        for(j=0;j<par.ny;j++){
            data.row(j)=ifft(datafx.row(j),par.nf);
        }
        for(j=0;j<par.nz;j++){
            data2.col(j)=data.col(j);
        }
        par.rebuildtx.row(i)=data2;
        par.realrebuildtx.row(i)=real(data2);
    } 
}

void tx_to_fx3d_wiener3d_pthread(int i, struct wiener3d * par){
    int j;
    cx_fmat datafx(par[0].ny,par[0].nf);
    fmat data(par[0].ny,par[0].nz);
    data=par[0].data.row(i);
    for(j=0;j<par[0].ny;j++){
        datafx.row(j)=fft(data.row(j),par[0].nf);
    }
    par[0].datafx.row(i)=datafx;
}
void tx_to_fx3d_wiener3d_thread(struct wiener3d & par){
    int i,k;
    thread *pcal;
    pcal=new thread[par.ncpu];
    
    for(i=0;i<(par.nx-par.ncpu);i+=par.ncpu){
    for(k=0;k<par.ncpu;k++){
        pcal[k]=thread(tx_to_fx3d_wiener3d_pthread,i+k,&par);
    }
    for(k=0;k<par.ncpu;k++){
        pcal[k].join();
    }
    } 
    for(i=(par.nx-par.ncpu);i<(par.nx);i++){
        k=i-(par.nx-par.ncpu);
        pcal[k]=thread(tx_to_fx3d_wiener3d_pthread,i,&par);
    }
    for(i=(par.nx-par.ncpu);i<(par.nx);i++){
        k=i-(par.nx-par.ncpu);
        pcal[k].join();
    }
}

void fx_to_tx3d_wiener3d_pthread(int i, struct wiener3d * par){
    cx_fmat data(par[0].ny,par[0].nf),datafx(par[0].ny,par[0].nf),\
        data2(par[0].ny,par[0].nz);
    int j;
    datafx=par[0].rebuildfx.row(i);
    for(j=0;j<par[0].ny;j++){
        data.row(j)=ifft(datafx.row(j),par[0].nf);
    }
    for(j=0;j<par[0].nz;j++){
        data2.col(j)=data.col(j);
    }
    par[0].rebuildtx.row(i)=data2;
    par[0].realrebuildtx.row(i)=real(data2);
}

void fx_to_tx3d_wiener3d_thread(struct wiener3d & par){
    int i,k;
    thread *pcal;
    pcal=new thread[par.ncpu];
    for(i=0;i<(par.nx-par.ncpu);i+=par.ncpu){
    for(k=0;k<par.ncpu;k++){
        pcal[k]=thread(fx_to_tx3d_wiener3d_pthread,i+k,&par);
    }
    for(k=0;k<par.ncpu;k++){
        pcal[k].join();
    }
    } 
    for(i=(par.nx-par.ncpu);i<(par.nx);i++){
        k=i-(par.nx-par.ncpu);
        pcal[k]=thread(fx_to_tx3d_wiener3d_pthread,i,&par);
    }
    for(i=(par.nx-par.ncpu);i<(par.nx);i++){
        k=i-(par.nx-par.ncpu);
        pcal[k].join();
    }
}

///////////////////thread use wiener3dCG//////////////////
struct wienerCG3d
{
    int nx,nz,npx,npz,numi,mini;
    float *a;
    cx_fmat *b_cxfmatp1ncpunpxy;
    cx_fmat **r_k1_cxfmatp1ncpunpxy;
    cx_fmat **r_k2_cxfmatp1ncpunpxy;
    cx_fmat **p_k1_cxfmatp1ncpunpxy;
    cx_fmat **p_k2_cxfmatp1ncpunpxy;
    cx_fmat **x_k1_cxfmatp1ncpunpxy;
    cx_fmat **x_k2_cxfmatp1ncpunpxy;
};
void wienerGC3d_parset(struct wiener3d & par, struct wienerCG3d & cg,\
    int numi=25, int mini=9);
void wienerGC3d_getfilter(cx_fmat & X, cx_fmat A, cx_fmat B, \
    struct wienerCG3d & cg, int kcpu);
void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu,\
    cx_fmat& arFilter);
void filter_xboundary(int kwinx, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu);
void filter_yboundary(int kwinx, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu);
void filter_xpoint(int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu);

void wiener3d_mid_pthread(struct wiener3d * par, int pnf1, int pnf2,\
    struct wienerCG3d * cg, int kcpu, cx_fmat* arFilter)
{
    int i,j,k;//cout<<"ok"<<endl;
    int nx(par[0].nx),nz(par[0].nz),ny(par[0].ny),nf(par[0].nf),\
        halfarx(par[0].halfnarx),halfary(par[0].halfnary),\
        nwx(par[0].nwx),nwy(par[0].nwy),movex(par[0].nmovex),\
        movey(par[0].nmovey),narx(par[0].narx),nary(par[0].nary);
    int nar(narx*nary),nw(nwx*nwy);
    fmat ar_weight(nx,ny);
    int kwinx,kwiny,kf,kout(0);

    for(kf=pnf1;kf<pnf2;kf++){
        ar_weight.fill(0.0);kout=0;
        //cout<<"now is run kf = "<<kf<<endl;
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex){
        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey){
            filter_mid(kwinx,kwiny,kf,par[0],ar_weight,cg[0],kcpu,\
                arFilter[0]);
            //filter_forword_back(kwinx,kwiny,kf,ppar[0],ar_weight);
            kout++;
        }
        filter_mid(kwinx,ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu,\
            arFilter[0]);
        kout++;
        //filter_forword_back(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight);
        }
        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey){
            filter_mid(nx-halfarx-nwx,kwiny,kf,par[0],ar_weight,cg[0],kcpu,\
                arFilter[0]);
            //filter_forword_back(nx-halfarx-nwx,kwiny,kf,ppar[0],ar_weight);
        }
        filter_mid(nx-halfarx-nwx,ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu,\
            arFilter[0]);
        //filter_forword_back(nx-halfarx-nwx,ny-halfary-nwy,kf,ppar[0],ar_weight);

        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey){
            filter_xboundary(kwiny,kf,par[0],ar_weight,cg[0],kcpu);
        }
        filter_xboundary(ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu);
        
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex){
            filter_yboundary(kwinx,kf,par[0],ar_weight,cg[0],kcpu);
        }
        filter_yboundary(nx-halfarx-nwx,kf,par[0],ar_weight,cg[0],kcpu);
        filter_xpoint(kf,par[0],ar_weight,cg[0],kcpu);

        for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            par[0].rebuildfx(i,j,kf).imag(2*imag(par[0].rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000000001));
            par[0].rebuildfx(i,j,kf).real(2*real(par[0].rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000000001));
        }}
    }
}

void wiener3d_mid_thread(struct wiener3d & par,\
    int midnumcpu, int numi=25, int mini=9)
{ 
    struct wiener3d * ppar;
    struct wienerCG3d cg;
    ppar=&par;
    getagc_wiener3d(ppar[0]);
    tx_to_fx3d_wiener3d_thread(ppar[0]);

    int ncpuout=par.ncpu;
    par.ncpu=midnumcpu;

    wienerGC3d_parset(ppar[0],cg,numi,mini);
    int numf(par.nf2-par.nf1);
    int pn(par.ncpu),pnf1,pnf2,k,kf;
    float dnf;
    thread *pcal;
    pcal=new thread[pn];
    dnf=float(par.nf2-par.nf1)/pn;
    
//Use for test, and ncpu=1 should be obey//
    //ofstream outfarreal,outfarimag,outfr;
    cx_fmat * arFilter;
    arFilter=new cx_fmat[pn];
///////////////////////////////////////////

    for(k=0;k<pn-1;k++){
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(wiener3d_mid_pthread,&par,pnf1,pnf2,\
            &cg, k, &arFilter[k]);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(wiener3d_mid_pthread,&par,pnf1,pnf2,\
        &cg, k, &arFilter[k]);

    for(k=0;k<pn;k++){
        pcal[k].join();
    }

    par.ncpu=ncpuout;
    fx_to_tx3d_wiener3d_thread(ppar[0]);
    deagc_wiener3d(ppar[0]);
}
 
/////////////////////////////////////////////////////////////////////////////
void wienerGC3d_parset(struct wiener3d & par, struct wienerCG3d & cg,\
    int numi, int mini)
{
    int i,j;//cout<<"ok"<<endl;
    int narx(par.narx),nary(par.nary);
    int nar(narx*nary-1),ncpu(par.ncpu);
    cg.numi=numi;cg.mini=mini;
    cg.r_k1_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.r_k2_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.p_k1_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.p_k2_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.x_k1_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    cg.x_k2_cxfmatp1ncpunpxy=new cx_fmat*[ncpu];
    for(i=0;i<ncpu;i++){
        cg.r_k1_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.r_k2_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.p_k1_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.p_k2_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.x_k1_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.x_k2_cxfmatp1ncpunpxy[i]=new cx_fmat[1];
        cg.r_k1_cxfmatp1ncpunpxy[i][0].zeros(nar,1);
        cg.r_k2_cxfmatp1ncpunpxy[i][0].zeros(nar,1);
        cg.p_k1_cxfmatp1ncpunpxy[i][0].zeros(nar,1);
        cg.p_k2_cxfmatp1ncpunpxy[i][0].zeros(nar,1);
        cg.x_k1_cxfmatp1ncpunpxy[i][0].zeros(nar,1);
        cg.x_k2_cxfmatp1ncpunpxy[i][0].zeros(nar,1);
    }
}
 
void wienerGC3d_getfilter(cx_fmat & X, cx_fmat A, cx_fmat B, \
    struct wienerCG3d & cg, int kcpu)
{
    int k,nar(A.n_cols);
    float errcg;
    cx_fmat meo(1,1),me(1,1),ak(1,1),bk(1,1),memat(nar,1);
    cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]=B-A*X;
    cg.p_k1_cxfmatp1ncpunpxy[kcpu][0]=\
        cg.r_k1_cxfmatp1ncpunpxy[kcpu][0];
    cg.x_k1_cxfmatp1ncpunpxy[kcpu][0]=X;
    errcg=0.000001*abs(sum(sum(A)))/A.n_elem;
    //errcg=0.000001;

    for(k=0;k<cg.numi;k++){ 
        memat = A*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
        me = cg.p_k1_cxfmatp1ncpunpxy[kcpu][0].st()*(memat.t()).st();
        ak = cg.r_k1_cxfmatp1ncpunpxy[kcpu][0].t()*\
            cg.r_k1_cxfmatp1ncpunpxy[kcpu][0];
        //cout<<"|"<<abs(me(0,0));
        ak(0,0)= ak(0,0)/me(0,0);
        if(k>cg.mini && abs(me(0,0))<errcg){
            //cout<<k<<"|";
            break;
        }
        //else{cout<<"|"<<real(ak(0,0));}
        cg.x_k2_cxfmatp1ncpunpxy[kcpu][0]=\
            cg.x_k1_cxfmatp1ncpunpxy[kcpu][0]+\
            ak(0,0)*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
        
        cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]=\
            cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]-ak(0,0)*memat;
        bk=(cg.r_k2_cxfmatp1ncpunpxy[kcpu][0].t()*\
            cg.r_k2_cxfmatp1ncpunpxy[kcpu][0])/\
            (cg.r_k1_cxfmatp1ncpunpxy[kcpu][0].t()*\
            cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]);
        cg.p_k2_cxfmatp1ncpunpxy[kcpu][0]=\
            cg.r_k2_cxfmatp1ncpunpxy[kcpu][0]+\
            bk(0,0)*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
        
        meo=me;
        memat = A*cg.p_k2_cxfmatp1ncpunpxy[kcpu][0];
        me = cg.p_k2_cxfmatp1ncpunpxy[kcpu][0].st()*(memat.t()).st();
        if(abs(me(0,0))<abs(meo(0,0))){  
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
            X=cg.x_k1_cxfmatp1ncpunpxy[kcpu][0];
        }
        else{
            X=cg.x_k2_cxfmatp1ncpunpxy[kcpu][0];
            cg.r_k1_cxfmatp1ncpunpxy[kcpu][0]=B-A*X;
            cg.p_k1_cxfmatp1ncpunpxy[kcpu][0]=\
                cg.r_k1_cxfmatp1ncpunpxy[kcpu][0];
            cg.x_k1_cxfmatp1ncpunpxy[kcpu][0]=X;
        }
    }   
}

void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu,\
    cx_fmat& arFilter)
{
    int i,j,kw(0),kar,kax,kay;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int nar(narx*nary-1),nw(nwx*nwy);
    float dig_n(par.dig_n);
    cx_fmat datafx(nx,ny);
    cx_fmat hankel_D(nw,nar),ar_S(nw,1),\
        filter_A(nar,1),R(nar,nar);
    fmat digw(nar,nar,fill::zeros),digw2(nar,nar,fill::zeros),\
        ardata(nar,1),realR(nar,nar);
    digw.diag()+=1; digw2.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);
    
    hankel_D.zeros(nw,nar),ar_S.zeros(nw,1),\
    R.zeros(nar,nar);
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        ar_S(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<=i+halfarx;kax++){
        for(kay=j-halfary;kay<=j+halfary;kay++){
            if((kax==i) && (kay==j))
            {continue;}
            else{
                hankel_D(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    R=hankel_D.t()*hankel_D;
    //dig_n=dig_n*(trace(abs(R)))/nar;
    dig_n=dig_n*max(max(abs(R)));
    digw=digw*dig_n;
 
    //filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    //void wienerGC3d_getfilter(cx_fmat & X, cx_fmat A, cx_fmat B, \
    struct wienerCG3d & cg, int kcpu)
    if(par.use_way==1){
        filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_A.zeros(nar,1),pcg=&cg;
        wienerGC3d_getfilter(filter_A, R+digw+digw2,\
            hankel_D.t()*ar_S, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_A=solvecxfmat((R+digw+digw2),hankel_D.t()*ar_S,nar);
    }
    else if(par.use_way==4){
        filter_A=invcxfmat(R+digw+digw2,nar)*hankel_D.t()*ar_S;
    }
    //filter_A=solve(R+digw,hankel_D.t()*ar_S);
    //output the AR-Wiener filter
    arFilter=filter_A;
    ar_S=hankel_D*filter_A;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_S(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}
}

void filter_xboundary(int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu)
{
    int i,j,kw(0),kar,kax,kay;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int narf((2*halfarx+1)*nary-1),\
        narb((2*halfarx+1)*nary-1),nw(nwx*nwy);
    float dig_n(par.dig_n);
    cx_fmat datafx(nx,ny);
    cx_fmat hankel_Df(nw,narf,fill::zeros),ar_Sf(nw,1,fill::zeros),\
        filter_Af(narf,1,fill::zeros),Rf(narf,narf,fill::zeros);
    cx_fmat hankel_Db(nw,narb,fill::zeros),ar_Sb(nw,1,fill::zeros),\
        filter_Ab(narb,1,fill::zeros),Rb(narb,narb,fill::zeros);
    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag()+=1;digw2f.diag(par.dig_n2); 
    fmat digwb(narb,narb,fill::zeros),digw2b(narb,narb,fill::zeros);
    digwb.diag()+=1;digw2b.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    kw=0;
    for(i=nx-nwx;i<nx;i++){
    for(j=kwiny;j<kwiny+nwy;j++)  {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++){
        for(kay=j-halfary;kay<=j+halfary;kay++){
            if((kax!=i) || (kay!=j)){
            hankel_Df(kw,kar)=datafx(kax,kay);
            kar++;
            }
        }}
        kw++;
    }}
    kw=0;
    for(i=0;i<nwx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        ar_Sb(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++){
        for(kay=j-halfary;kay<=j+halfary;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Db(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1){
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Af.zeros(narf,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
            hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4){
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}

    Rb=hankel_Db.t()*hankel_Db;
    float dig_nb(dig_n*max(max(abs(Rb))));
    digwb=digwb*dig_nb;
    //filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    if(par.use_way==1){
        filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Ab.zeros(narb,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Ab, Rb+digwb+digw2b,\
            hankel_Db.t()*ar_Sb, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Ab=solvecxfmat((Rb+digwb+digw2b),\
            hankel_Db.t()*ar_Sb,narb);
    }
    else if(par.use_way==4){
        filter_Ab=invcxfmat(Rb+digwb+digw2b,narb)\
            *hankel_Db.t()*ar_Sb;
    }

    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=0;i<nwx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        //if(i<halfarx)
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}
}

void filter_yboundary(int kwinx, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu)
{
    int i,j,kw(0),kar,kax,kay;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int narf((2*halfary+1)*narx-1),\
        narb((2*halfary+1)*narx-1),nw(nwy*nwx);
    float dig_n(par.dig_n);
    cx_fmat datafx(nx,ny);
    cx_fmat hankel_Df(nw,narf,fill::zeros),ar_Sf(nw,1,fill::zeros),\
        filter_Af(narf,1,fill::zeros),Rf(narf,narf,fill::zeros);
    cx_fmat hankel_Db(nw,narb,fill::zeros),ar_Sb(nw,1,fill::zeros),\
        filter_Ab(narb,1,fill::zeros),Rb(narb,narb,fill::zeros);
    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag()+=1;digw2f.diag(par.dig_n2); 
    fmat digwb(narb,narb,fill::zeros),digw2b(narb,narb,fill::zeros);
    digwb.diag()+=1;digw2b.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=ny-nwy;j<ny;j++){
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<=i+halfarx;kax++){
        for(kay=j-2*halfary;kay<=j;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Df(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=0;j<nwy;j++){
        ar_Sb(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<=i+halfarx;kax++){
        for(kay=j;kay<=j+2*halfary;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Db(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1){
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Af.zeros(narf,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
            hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4){
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=ny-nwy;j<ny;j++){
        //if(j>=ny-halfary)
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}

    Rb=hankel_Db.t()*hankel_Db;
    float dig_nb(dig_n*max(max(abs(Rb))));
    digwb=digwb*dig_nb;
    //filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    if(par.use_way==1){
        filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Ab.zeros(narb,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Ab, Rb+digwb+digw2b,\
            hankel_Db.t()*ar_Sb, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Ab=solvecxfmat((Rb+digwb+digw2b),\
            hankel_Db.t()*ar_Sb,narb);
    }
    else if(par.use_way==4){
        filter_Ab=invcxfmat(Rb+digwb+digw2b,narb)\
            *hankel_Db.t()*ar_Sb;
    }

    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=0;j<nwy;j++){
        //if(j<halfary)
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}
}

void filter_xpoint(int kf,\
    struct wiener3d & par, fmat & weightmat,\
    struct wienerCG3d & cg, int kcpu)
{
    int i,j,kw(0),kar,kax,kay;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int narf((2*halfarx+1)*nary-1),\
        narb((2*halfarx+1)*nary-1),nw(nwx*nwy);
    float dig_n(par.dig_n);
    cx_fmat datafx(nx,ny);
    cx_fmat hankel_Df(nw,narf),ar_Sf(nw,1),\
        filter_Af(narf,1),Rf(narf,narf);

    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag()+=1;digw2f.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=nx-nwx;i<nx;i++){
    for(j=ny-nwy;j<ny;j++){
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++){
        for(kay=j-2*halfary;kay<=j;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Df(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1){
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Af.zeros(narf,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
            hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4){
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }
    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++){
    for(j=ny-nwy;j<ny;j++){
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=nx-nwx;i<nx;i++){
    for(j=0;j<nwy;j++){
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++){
        for(kay=j;kay<=j+2*halfary;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Df(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1){
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Af.zeros(narf,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
            hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4){
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++){
    for(j=0;j<nwy;j++){
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=0;i<nwx;i++){
    for(j=ny-nwy;j<ny;j++){
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++){
        for(kay=j-2*halfary;kay<=j;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Df(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1){
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Af.zeros(narf,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
            hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4){
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=0;i<nwx;i++){
    for(j=ny-nwy;j<ny;j++) 
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=0;i<nwx;i++){
    for(j=0;j<nwy;j++){
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++){
        for(kay=j;kay<=j+2*halfary;kay++){
            if((kax!=i) || (kay!=j)){
                hankel_Df(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }}
        kw++;
    }}
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1){
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_Af.zeros(narf,1),pcg=&cg;
        wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
            hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4){
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=0;i<nwx;i++){
    for(j=0;j<nwy;j++){
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax){
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }}
}

inline cx_fmat invcxfmat(cx_fmat mat1, int n)
{
    fmat matA(n,n), matB(n,n);
    fmat matAinv(n,n), matABinv(n,n);
    matA=real(mat1);
    matB=imag(mat1);
    matAinv=invfmat(matA,n);
    matABinv=invfmat(matA+matB*matAinv*matB,n);

    cx_fmat mat2(n,n);
    mat2.set_real(matABinv);
    mat2.set_imag(-matAinv*matB*matABinv);
    return mat2;
}

inline fmat invfmat(fmat mat1, int n)
{
    fmat mat2(n,n,fill::zeros);
    mat2.diag()+=1.0;
    int k1,k2;
    float xs,max_mat_elem;
    max_mat_elem=abs(mat1).max();
    mat1=mat1/max_mat_elem;
    //cout<<max_mat_elem<<'|';
    for(k2=0;k2<n;k2++){
    for(k1=1+k2;k1<n;k1++){
        if(abs(mat1(k1,k2))>0.000001){
            xs=mat1(k2,k2)/mat1(k1,k2);
            mat1.row(k1)=xs*mat1.row(k1)-mat1.row(k2);
            mat2.row(k1)=xs*mat2.row(k1)-mat2.row(k2);
        }
    }}
    for(k2=0;k2<n;k2++){
    for(k1=0;k1<k2;k1++){
        if(abs(mat1(k1,k2))>0.000001){
            xs=mat1(k2,k2)/mat1(k1,k2);
            mat1.row(k1)=xs*mat1.row(k1)-mat1.row(k2);
            mat2.row(k1)=xs*mat2.row(k1)-mat2.row(k2);
        }
    }}
    for(k2=0;k2<n;k2++){
        mat2.row(k2)=mat2.row(k2)/mat1(k2,k2);
    }
    mat2=mat2/max_mat_elem;
    //cout<<mat2.has_inf()<<","<<mat2.has_nan()<<"|";
    
    for(k2=0;k2<n;k2++){
    for(k1=0;k1<n;k1++){
        if(isnan(mat2(k1,k2)))
            {mat2(k1,k2)=0;}
        if(isinf(mat2(k1,k2)))
            {mat2(k1,k2)=0;}
    }}
    
    return mat2;
}

inline cx_fmat solvecxfmat(cx_fmat mata, cx_fmat matb, int n)
{
    cx_fmat mat2(n,1,fill::zeros);
    int k1,k2;
    cx_float xs;
    float max_mata_ele;
    max_mata_ele=0.000001*abs(mata).max();

    for(k2=0;k2<n;k2++){
    for(k1=1+k2;k1<n;k1++){
        if(abs(mata(k1,k2))>max_mata_ele){
            xs=mata(k2,k2)/mata(k1,k2);
            mata.row(k1)=xs*mata.row(k1)-mata.row(k2);
            matb.row(k1)=xs*matb.row(k1)-matb.row(k2);
        }
    }}

    for(k2=0;k2<n;k2++){
        if(isnan(abs(matb(k2,0)))){
            matb(k2,0).imag(0);
            matb(k2,0).real(0);
        }
        if(isinf(abs(matb(k2,0)))){
            matb(k2,0).imag(0);
            matb(k2,0).real(0);
        }
        for(k1=0;k1<n;k1++){
            if(isnan(abs(mata(k1,k2)))){
                mata(k1,k2).imag(0);
                mata(k1,k2).real(0);
            }
            if(isinf(abs(mata(k1,k2)))){
                mata(k1,k2).imag(0);
                mata(k1,k2).real(0);
            }
        }
    }

    for(k2=n-1;k2>=0;k2--){
        xs=matb(k2,0);
        for(k1=n-1;k1>k2;k1--){
            xs=xs-mat2(k1,0)*mata(k2,k1);
        }
        mat2(k2,0)=xs/mata(k2,k2);
    }

    for(k2=0;k2<n;k2++){
        if(isnan(abs(mat2(k2,0)))){
            mat2(k2,0).imag(0);
            mat2(k2,0).real(0);
        }
        if(isinf(abs(mat2(k2,0)))){
            mat2(k2,0).imag(0);
            mat2(k2,0).real(0);
        }
    }
    return mat2;
}

inline fmat solvefmat(fmat mata, fmat matb, int n)
{
    fmat mat2(n,1,fill::zeros);
    int k1,k2;
    float xs;
    float max_mata_ele;
    max_mata_ele=0.000001*abs(mata).max();

    for(k2=0;k2<n;k2++){
    for(k1=1+k2;k1<n;k1++){
        if(abs(mata(k1,k2))>max_mata_ele){
            xs=mata(k2,k2)/mata(k1,k2);
            mata.row(k1)=xs*mata.row(k1)-mata.row(k2);
            matb.row(k1)=xs*matb.row(k1)-matb.row(k2);
        }
    }}

    for(k2=0;k2<n;k2++){
        if(isnan(abs(matb(k2,0)))){
            matb(k2,0)=(0);
        }
        if(isinf(abs(matb(k2,0)))){
            matb(k2,0)=(0);
        }
        for(k1=0;k1<n;k1++){
            if(isnan(abs(mata(k1,k2)))){
                mata(k1,k2)=(0);
            }
            if(isinf(abs(mata(k1,k2)))){
                mata(k1,k2)=(0);
            }
        }
    }

    for(k2=n-1;k2>=0;k2--){
        xs=matb(k2,0);
        for(k1=n-1;k1>k2;k1--){
            xs=xs-mat2(k1,0)*mata(k2,k1);
        }
        mat2(k2,0)=xs/mata(k2,k2);
    }

    for(k2=0;k2<n;k2++){
        if(isnan(abs(mat2(k2,0)))){
            mat2(k2,0)=(0);
        }
        if(isinf(abs(mat2(k2,0)))){
            mat2(k2,0)=(0);
        }
    }
    return mat2;
}

void fxWiener3dFilteRandomNoise(fcube& data3dOrig,\
    fcube& data3dResult,\
    float dt, float fmin, float fmax, int ncpu=1,\
    int nwx=21, int nwy=21, int narx=5, int nary=5,\
    float TikhonovReg=0.001,float digNum=0.001)
{
    int nx(data3dOrig.n_rows),ny(data3dOrig.n_cols),nz(data3dOrig.n_slices);
    struct wiener3d par;           //结构体，包含数据和参数
    wiener3d_parset(nx,ny,nz,par); //结构体初始化函数,并设置一组初始化参数
    par.nwx=nwx;par.nwy=nwy;          //数据窗求解误差逼近泛函
    par.narx=narx;par.nary=nary;         //wiener滤波器系数阶数
    par.dig_n=TikhonovReg;               //法方程求解对角稳定系数
    par.dig_n2=digNum;
    par.nmovex=par.nwx;            //滑动滤波的滑动步长，
    par.nmovey=par.nwy;            //这里设置为等于数据窗大小，避免多次覆盖滤波
    par.covermax=1;
    float df=1.0/nz/dt;
    par.nf1=fmin/df,par.nf2=fmax/df;    //滤波频率范围
    par.ncpu=ncpu;
    wiener3d_parupdate(par);       //调整参数后更新结构体
    par.data=data3dOrig;

    par.use_way=1;
    wiener3d_mid_thread(par,ncpu); 
    wiener3d_cleardata(par);
    data3dResult=par.realrebuildtx;
}

///////////////////////////////////////////////////////////////
void interpolationMidByRow(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat,fmat & weightmatInter,\
    struct wienerCG3d & cg, int kcpu,\
    cx_fmat& arFilter)
{
    int i,j,kw(0),kar,kax,kay;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int nar((narx-1)*nary),nw(nwx*nwy);
    float dig_n(par.dig_n);
    cx_fmat datafx(nx,ny),datafxInter(nx,ny);
    cx_fmat hankel_D(nw,nar),hankelInter(nw,nar),\
        ar_S(nw,1),arInter(nw,1),\
        filter_A(nar,1),R(nar,nar);
    fmat digw(nar,nar,fill::zeros),digw2(nar,nar,fill::zeros),\
        ardata(nar,1),realR(nar,nar);
    digw.diag()+=1; 

    datafx=par.datafx.slice(kf);
    datafxInter=par.datafx.slice(2*kf);
    hankel_D.zeros(nw,nar),ar_S.zeros(nw,1),\
    R.zeros(nar,nar);
    
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        ar_S(kw,0)=datafx(i,j);
        arInter(kw,0)=datafxInter(i,j);
        kar=0;
        for(kax=i-2*halfarx+1;kax<=i+2*halfarx-1;kax+=2){
            int yloop;
            if(j+2*halfary>=0){yloop=j+2*halfary;}
            else{yloop=j+halfary;}
        for(kay=max(j-2*halfary,0);kay<=yloop;kay+=2){
            //if((kax==i) && (kay==j))
            //{continue;}
            //else{
                hankel_D(kw,kar)=datafx(kax,kay);
                kar++;
            //}
        }}
        kar=0;
        for(kax=i-halfarx+1;kax<=i+halfarx;kax+=1){
        for(kay=j-halfary;kay<=j+halfary;kay+=1){
            //if((kax==i) && (kay==j))
            //{continue;}
            //else{
                hankelInter(kw,kar)=datafxInter(kax,kay);
                kar++;
            //}
        }}
        kw++;
    }}
    R=hankel_D.t()*hankel_D;
    //dig_n=dig_n*(trace(abs(R)))/nar;
    dig_n=dig_n*max(max(abs(R)));
    digw=digw*dig_n;
    if(dig_n<par.dig_n2){
        digw2.diag(par.dig_n2);
    }
    //filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    //void wienerGC3d_getfilter(cx_fmat & X, cx_fmat A, cx_fmat B, \
    struct wienerCG3d & cg, int kcpu)
    if(par.use_way==1){
        filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    }
    else if(par.use_way==2){
        struct wienerCG3d * pcg;
        filter_A.zeros(nar,1),pcg=&cg;
        wienerGC3d_getfilter(filter_A, R+digw+digw2,\
            hankel_D.t()*ar_S, pcg[0], kcpu);
    }
    else if(par.use_way==3){
        filter_A=solvecxfmat((R+digw+digw2),hankel_D.t()*ar_S,nar);
    }
    else if(par.use_way==4){
        filter_A=invcxfmat(R+digw+digw2,nar)*hankel_D.t()*ar_S;
    }
    //filter_A=solve(R+digw,hankel_D.t()*ar_S);
    //output the AR-Wiener filter
    arFilter=filter_A;
    ar_S=hankel_D*filter_A;
    arInter=hankelInter*filter_A;
    kw=0;

    for(i=kwinx;i<kwinx+nwx;i++){
    for(j=kwiny;j<kwiny+nwy;j++){
        //if(weightmat(i,j)<par.covermax){
            par.fx3dInterpolation((2*i+1),j,2*kf)+=arInter(kw,0);
            weightmatInter((2*i+1),j)+=1;
        //}
        kw++;
    }}
}

void wiener3dMidInterpolationByRow_pthread(\
    struct wiener3d * par, int pnf1, int pnf2,\
    struct wienerCG3d * cg, int kcpu, cx_fmat* arFilter)
{
    int i,j,k;//cout<<"ok"<<endl;
    int nx(par[0].nx),nz(par[0].nz),ny(par[0].ny),nf(par[0].nf),\
        halfarx(par[0].halfnarx),halfary(par[0].halfnary),\
        nwx(par[0].nwx),nwy(par[0].nwy),movex(par[0].nmovex),\
        movey(par[0].nmovey),narx(par[0].narx),nary(par[0].nary);
    int nar(narx*nary),nw(nwx*nwy);
    fmat ar_weight(par->fx3dInterpolation.n_rows,ny);
    fmat ar_weight_inter(par->fx3dInterpolation.n_rows,ny);
    int kwinx,kwiny,kf,kout(0);
    arFilter[0].zeros(nar,1);

    for(kf=pnf1;kf<pnf2;kf++){
        ar_weight.fill(0.0);
        ar_weight_inter.fill(0.0);
        kout=0;
        cout<<"now is run kf = "<<kf<<endl;
        for(kwinx=2*halfarx-1;kwinx<nx-2*halfarx-nwx-1;kwinx+=movex){
            int yloop;
            if(ny-2*halfary-nwy>0){yloop=ny-2*halfary-nwy;}
            else{yloop=ny-halfary-nwy;}
            //cout<<yloop<<endl;
        for(kwiny=max(2*halfary,0);kwiny<yloop;kwiny+=movey){
            interpolationMidByRow(kwinx,kwiny,kf,par[0],ar_weight,\
                ar_weight_inter,cg[0], kcpu, arFilter[0]);
            //filter_forword_back(kwinx,kwiny,kf,ppar[0],ar_weight);
            kout++;
        }
        interpolationMidByRow(kwinx,yloop,kf,par[0],ar_weight,\
                ar_weight_inter,cg[0], kcpu, arFilter[0]);
        kout++;
        //filter_forword_back(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight);
        }
        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey){
            filter_mid(nx-halfarx-nwx,kwiny,kf,par[0],ar_weight,cg[0],kcpu,\
                arFilter[0]);
            //filter_forword_back(nx-halfarx-nwx,kwiny,kf,ppar[0],ar_weight);
        }
        filter_mid(nx-halfarx-nwx,ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu,\
            arFilter[0]);
        //filter_forword_back(nx-halfarx-nwx,ny-halfary-nwy,kf,ppar[0],ar_weight);

        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey){
            filter_xboundary(kwiny,kf,par[0],ar_weight,cg[0],kcpu);
        }
        filter_xboundary(ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu);
        
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex){
            filter_yboundary(kwinx,kf,par[0],ar_weight,cg[0],kcpu);
        }
        filter_yboundary(nx-halfarx-nwx,kf,par[0],ar_weight,cg[0],kcpu);
        filter_xpoint(kf,par[0],ar_weight,cg[0],kcpu);
        
        for(i=1;i<nx;i++){
        for(j=0;j<ny;j++){
            par[0].fx3dInterpolation(2*i-1,j,2*kf).imag\
                (1*imag(par[0].fx3dInterpolation(2*i-1,j,2*kf))\
                /(ar_weight_inter(2*i-1,j)+0.000000001));
            par[0].fx3dInterpolation(2*i-1,j,2*kf).real\
                (1*real(par[0].fx3dInterpolation(2*i-1,j,2*kf))\
                /(ar_weight_inter(2*i-1,j)+0.000000001));
        }}
    }
}

void wiener3dMidInterpolationByRowThread(struct wiener3d & par,\
    int midnumcpu, int numi=25, int mini=9)
{ 
    struct wiener3d * ppar;
    struct wienerCG3d cg;
    ppar=&par;
    getagc_wiener3d(ppar[0]);
    tx_to_fx3d_wiener3d_thread(ppar[0]);
    //par.fx3dInterpolation.zeros(par.datafx.n_rows*2-1,\
        par.datafx.n_cols,par.datafx.n_slices);
    par.fx3dInterpolation=par.datafx;

    cxfcubeLinearInterpolation3dByRow(par.fx3dInterpolation);
    cout<<par.nx<<"ok"<<par.fx3dInterpolation.n_rows<<endl;

    for(int k=1;k<par.nx;k++){
        par.fx3dInterpolation.row(2*k-1).fill(0.0);
    }

    int ncpuout=par.ncpu;
    par.ncpu=midnumcpu;

    wienerGC3d_parset(ppar[0],cg,numi,mini);
    int pn(par.ncpu),pnf1,pnf2,k,kf;
    float dnf;
    thread *pcal;
    pcal=new thread[pn];
    dnf=float(par.nf2-par.nf1)/pn;
    
//Use for test, and ncpu=1 should be obey//
    //ofstream outfarreal,outfarimag,outfr;
    cx_fmat * arFilter;
    arFilter=new cx_fmat[pn];
///////////////////////////////////////////

    for(k=0;k<pn-1;k++){
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(wiener3dMidInterpolationByRow_pthread,\
            &par,pnf1,pnf2,&cg, k, &arFilter[k]);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(wiener3dMidInterpolationByRow_pthread,\
        &par,pnf1,pnf2,&cg, k, &arFilter[k]);

    for(k=0;k<pn;k++){
        pcal[k].join();
    }
    par.rebuildfx.zeros(par.fx3dInterpolation.n_rows,\
        par.fx3dInterpolation.n_cols,\
        par.fx3dInterpolation.n_slices/2);
    for(int kf=par.nf1;kf<par.nf2;kf++){
        par.rebuildfx.slice(kf)=par.fx3dInterpolation.slice(2*kf);
    }

    par.ncpu=ncpuout;
    par.realrebuildtx.copy_size(par.rebuildfx);
    fx2tx_3d_thread(par.realrebuildtx,par.rebuildfx,par.ncpu);
    fcubeLinearInterpolation3dByRow(par.dataagc);
    par.realrebuildtx*=2.0;
    cout<<"finished wiener"<<endl;
    deagc_wiener3d(ppar[0]);
    cout<<"finished deagc"<<endl;
}

void fxWiener3dInterpolationByRow(fcube& data3dOrig,\
    fcube& data3dResult,\
    float dt, float fmin, float fmax, int ncpu=1,\
    int nwx=21, int nwy=21, int narx=5, int nary=5,\
    float TikhonovReg=0.000001,float digNum=0.000001)
{
    int nx(data3dOrig.n_rows),ny(data3dOrig.n_cols),\
        nz(data3dOrig.n_slices*2);
    struct wiener3d par;           //结构体，包含数据和参数
    wiener3d_parset(nx,ny,nz,par); //结构体初始化函数,并设置一组初始化参数
    par.nwx=nwx;par.nwy=nwy;          //数据窗求解误差逼近泛函
    par.narx=narx;par.nary=nary;         //wiener滤波器系数阶数
    par.dig_n=TikhonovReg;               //法方程求解对角稳定系数
    par.dig_n2=digNum;
    par.nmovex=par.nwx;            //滑动滤波的滑动步长，
    par.nmovey=par.nwy;            //这里设置为等于数据窗大小，避免多次覆盖滤波
    par.covermax=1;
    float df=1.0/nz/dt;
    par.nf1=fmin/df,par.nf2=fmax/df;    //滤波频率范围
    par.ncpu=ncpu;
    wiener3d_parupdate(par);       //调整参数后更新结构体
    for(int k=0;k<data3dOrig.n_slices;k++){
        par.data.slice(k)=data3dOrig.slice(k);
    }
////////////////////////////////
    par.use_way=1;
    wiener3dMidInterpolationByRowThread(par,ncpu); 
    wiener3d_cleardata(par);
///////////////////////////////
    data3dResult=par.realrebuildtx;
}

#endif
