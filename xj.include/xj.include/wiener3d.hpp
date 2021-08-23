/************update in 2021.01.07**************/

    
    
/*********************************************/

#ifndef WIENER3D_HPP
#define WIENER3D_HPP

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

struct wiener3d
{
    int nz,nx,ny,nf,nf1,nf2,nwx,nwy,narx,nary,\
        halfnarx,halfnary,nbackx,nbacky;
    float dz,dy,dx,df,dig_n;
    fcube data,realrebuildtx;
    cx_fcube rebuildfx,rebuildtx,datafx;
};

//设置一组常用的初始化参数
//Mat(nx,ny,nz); Mat(nx,ny,nf)
void wiener3d_parset(int nx, int ny, int nz,int nf,\
    struct wiener3d & par)
{
    par.nx=nx,par.ny=ny,par.nz=nz,par.nf=nf;
    par.nf1=0,par.nf2=nf/2;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.df=1.0/par.dz/par.nf;
    par.dig_n=0.01;
    //par.datafx.zeros(par.nx,par.ny,par.nf);
    par.data.zeros(par.nx,par.ny,par.nz);
    //par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nf);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nf);
}

void wiener3d_parupdate(struct wiener3d & par)
{
    par.df=1.0/par.dz/par.nf;
    par.halfnarx=int(par.narx/2);
    par.narx=2*par.halfnarx+1;
    par.halfnary=int(par.nary/2);
    par.nary=2*par.halfnary+1;
    par.datafx.zeros(par.nx,par.ny,par.nf);
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
}

void wiener3d_cleardata(struct wiener3d & par)
{
    par.data.zeros(par.nx,par.ny,par.nz);
    par.rebuildtx.zeros(par.nx,par.ny,par.nf);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nf);
}

void tx_to_fx3d_wiener3d(struct wiener3d & par)
{
    int i,j;
    fmat data(par.ny,par.nz);
    cx_fmat datafx(par.ny,par.nf);
    for(i=0;i<par.nx;i++)
    {
        data=par.data.row(i);
        for(j=0;j<par.ny;j++)
        {
            datafx.row(j)=fft(data.row(j),par.nf);
        }
        par.datafx.row(i)=datafx;
    } 
}

void fx_to_tx3d_wiener3d(struct wiener3d & par)
{
    int i,j;
    cx_fmat data(par.ny,par.nf),datafx(par.ny,par.nf);
    for(i=0;i<par.nx;i++)
    {
        datafx=par.rebuildfx.row(i);
        for(j=0;j<par.ny;j++)
        {
            data.row(j)=ifft(datafx.row(j),par.nf);
        }
        par.rebuildtx.row(i)=data;
        par.realrebuildtx.row(i)=real(data);
    } 
}

void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par)
{
    cout<<"now is run kf = "<<kf<<endl;
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
    fmat digw(nar,nar,fill::eye); 
    datafx=par.datafx.slice(kf);

    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        ar_S(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<=i+halfarx;kax++)
        {
        for(kay=j-halfary;kay<=j+halfary;kay++)  
        {
            if(kax!=i && kay!=j)
            {
                hankel_D(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }
        }
        kw++;
    }
    }
    R=hankel_D.t()*hankel_D;
    dig_n=dig_n*(R(0,0).real());
    digw=digw*dig_n;
    filter_A=inv(R+digw)*hankel_D.t()*ar_S;
    ar_S=hankel_D*filter_A;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        par.rebuildfx(i,j,kf)=ar_S(kw,0);
        kw++;
    }
    }
}

void wiener3d_mid(struct wiener3d & par)
{
    struct wiener3d * ppar;
    ppar=&par;
    tx_to_fx3d_wiener3d(ppar[0]);
    int i,j,k;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int nar(narx*nary),nw(nwx*nwy);

    int kwinx,kwiny,kf;
    for(kf=par.nf1;kf<par.nf2;kf++)
    {
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=nwx)
        {
        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=nwy)
        {
            filter_mid(kwinx,kwiny,kf,ppar[0]);
        }
        }

    }
    fx_to_tx3d_wiener3d(ppar[0]);
}

#endif
