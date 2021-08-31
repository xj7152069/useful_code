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
        halfnarx,halfnary,nmovex,nmovey,ncpu;
    float dz,dy,dx,df,dig_n,dig_n2;
    fcube data,realrebuildtx;
    cx_fcube rebuildfx,rebuildtx,datafx;
};

//设置一组常用的初始化参数
//Mat(nx,ny,nz); Mat(nx,ny,nf)
void wiener3d_parset(int nx, int ny, int nz,\
    struct wiener3d & par)
{
    par.nx=nx,par.ny=ny,par.nz=nz;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.dig_n=0.001;par.dig_n2=0.00;
    par.nwx=nx/2;par.nwy=ny/2;          
    par.narx=nx/4;par.nary=ny/4;
    par.nmovex=1;par.nmovey=1;
    //par.datafx.zeros(par.nx,par.ny,par.nf);
    par.data.zeros(par.nx,par.ny,par.nz);
    //par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nz);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
    par.ncpu=1;par.nf=1;
    par.df=1.0/par.dz/par.nf;
    par.halfnarx=int(par.narx/2);
    par.narx=2*par.halfnarx+1;
    par.halfnary=int(par.nary/2);
    par.nary=2*par.halfnary+1;
    while (par.nf<par.nz)
    {
        par.nf*=2;
    }
    cout<<"nf = "<<par.nf<<endl;;
    par.nf1=0,par.nf2=par.nf/2;
    par.df=1.0/par.dz/par.nf;
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
    cx_fmat data(par.ny,par.nf),datafx(par.ny,par.nf),data2(par.ny,par.nz);
    for(i=0;i<par.nx;i++)
    {
        datafx=par.rebuildfx.row(i);
        for(j=0;j<par.ny;j++)
        {
            data.row(j)=ifft(datafx.row(j),par.nf);
        }
        for(j=0;j<par.nz;j++)
        {
            data2.col(j)=data.col(j);
        }
        par.rebuildtx.row(i)=data2;
        par.realrebuildtx.row(i)=real(data2);
    } 
}

void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat);
void filter_forword_back(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat);
void filter_xboundary(int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat);
void filter_yboundary(int kwinx, int kf,\
    struct wiener3d & par, fmat & weightmat);

void wiener3d_mid(struct wiener3d & par)
{
    struct wiener3d * ppar;
    ppar=&par;
    tx_to_fx3d_wiener3d(ppar[0]);
    int i,j,k;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),movex(par.nmovex),\
        movey(par.nmovey),narx(par.narx),nary(par.nary);
    int nar(narx*nary),nw(nwx*nwy);
    fmat ar_weight(nx,ny);

    int kwinx,kwiny,kf;
    for(kf=par.nf1;kf<par.nf2;kf++)
    {
        ar_weight.fill(0.0);
        cout<<"now is run kf = "<<kf<<endl;
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex)
        {
            for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
            {
                filter_mid(kwinx,kwiny,kf,ppar[0],ar_weight);
                //filter_forword_back(kwinx,kwiny,kf,ppar[0],ar_weight);
            }
            filter_mid(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight);
            //filter_forword_back(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight);
        }
        filter_mid(nx-halfarx-nwx,ny-halfary-nwy,kf,ppar[0],ar_weight);
        //filter_forword_back(nx-halfarx-nwx,ny-halfary-nwy,kf,ppar[0],ar_weight);

        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
        {
            filter_xboundary(kwiny,kf,ppar[0],ar_weight);
        }
        filter_xboundary(ny-halfary-nwy,kf,ppar[0],ar_weight);
        
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex)
        {
            filter_yboundary(kwinx,kf,ppar[0],ar_weight);
        }
        filter_yboundary(nx-halfarx-nwx,kf,ppar[0],ar_weight);

        for(i=0;i<nx;i++)
        {
        for(j=0;j<ny;j++)  
        {
            par.rebuildfx(i,j,kf).imag(2*imag(par.rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000001));
            par.rebuildfx(i,j,kf).real(2*real(par.rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000001));
        }
        }
    }
    fx_to_tx3d_wiener3d(ppar[0]);
}

/////////////////////////////////////////////////////////////////////////////
void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat)
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
    fmat digw(nar,nar,fill::zeros),digw2(nar,nar,fill::zeros);
    digw.diag(1.0); digw2.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);
    
    kw=0;
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
            if((kax==i) && (kay==j))
            {continue;}
            else
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
    //dig_n=dig_n*(trace(abs(R)))/nar;
    dig_n=dig_n*max(max(abs(R)));
    digw=digw*dig_n;
    filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    //filter_A=solve(R+digw,hankel_D.t()*ar_S);
    ar_S=hankel_D*filter_A;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        par.rebuildfx(i,j,kf)+=ar_S(kw,0);
        weightmat(i,j)+=1;
        kw++;
    }
    }
}

void filter_forword_back(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat)
{
    int i,j,kw(0),kar,kax,kay;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),\
        narx(par.narx),nary(par.nary);
    int narf(halfarx*nary),narb((halfarx+1)*nary-1),nw(nwx*nwy);
    float dig_n(par.dig_n);
    cx_fmat datafx(nx,ny);
    cx_fmat hankel_Df(nw,narf),ar_Sf(nw,1),\
        filter_Af(narf,1),Rf(narf,narf);
    cx_fmat hankel_Db(nw,narb),ar_Sb(nw,1),\
        filter_Ab(narb,1),Rb(narb,narb);
    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag(1.0);digw2f.diag(par.dig_n2); 
    fmat digwb(narb,narb,fill::zeros),digw2b(narb,narb,fill::zeros);
    digwb.diag(1.0);digw2b.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<i;kax++)
        {
        for(kay=j-halfary;kay<=j+halfary;kay++)  
        {
            hankel_Df(kw,kar)=datafx(kax,kay);
            kar++;
        }
        }
        ar_Sb(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+halfarx;kax++)
        {
        for(kay=j-halfary;kay<=j+halfary;kay++)  
        {
            if((kax!=i) || (kay!=j))
            {
                hankel_Db(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }
        }
        kw++;
    }
    }
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
        weightmat(i,j)+=1;
        kw++;
    }
    }

    Rb=hankel_Db.t()*hankel_Db;
    float dig_nb(dig_n*max(max(abs(Rb))));
    digwb=digwb*dig_nb;
    filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
        weightmat(i,j)+=1;
        kw++;
    }
    }
}

void filter_xboundary(int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat)
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
    cx_fmat hankel_Db(nw,narb),ar_Sb(nw,1),\
        filter_Ab(narb,1),Rb(narb,narb);
    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag(1.0);digw2f.diag(par.dig_n2); 
    fmat digwb(narb,narb,fill::zeros),digw2b(narb,narb,fill::zeros);
    digwb.diag(1.0);digw2b.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++)
        {
        for(kay=j-halfary;kay<=j+halfary;kay++)  
        {
            if((kax!=i) || (kay!=j))
            {
            hankel_Df(kw,kar)=datafx(kax,kay);
            kar++;
            }
        }
        }
        kw++;
    }}
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        ar_Sb(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++)
        {
        for(kay=j-halfary;kay<=j+halfary;kay++)  
        {
            if((kax!=i) || (kay!=j))
            {
                hankel_Db(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }
        }
        kw++;
    }
    }
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<1)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    Rb=hankel_Db.t()*hankel_Db;
    float dig_nb(dig_n*max(max(abs(Rb))));
    digwb=digwb*dig_nb;
    filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        //if(i<halfarx)
        if(weightmat(i,j)<1)
        {
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }
}

void filter_yboundary(int kwinx, int kf,\
    struct wiener3d & par, fmat & weightmat)
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
    cx_fmat hankel_Df(nw,narf),ar_Sf(nw,1),\
        filter_Af(narf,1),Rf(narf,narf);
    cx_fmat hankel_Db(nw,narb),ar_Sb(nw,1),\
        filter_Ab(narb,1),Rb(narb,narb);
    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag(1.0);digw2f.diag(par.dig_n2); 
    fmat digwb(narb,narb,fill::zeros),digw2b(narb,narb,fill::zeros);
    digwb.diag(1.0);digw2b.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++)
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<=i+halfarx;kax++)
        {
        for(kay=j-2*halfary;kay<=j;kay++)  
        {
            if((kax!=i) || (kay!=j))
            {
            hankel_Df(kw,kar)=datafx(kax,kay);
            kar++;
            }
        }
        }
        kw++;
    }}
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=0;j<nwy;j++)
    {
        ar_Sb(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-halfarx;kax<=i+halfarx;kax++)
        {
        for(kay=j;kay<=j+2*halfary;kay++) 
        {
            if((kax!=i) || (kay!=j))
            {
                hankel_Db(kw,kar)=datafx(kax,kay);
                kar++;
            }
        }
        }
        kw++;
    }
    }
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++)
    {
        //if(j>=ny-halfary)
        if(weightmat(i,j)<1)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    Rb=hankel_Db.t()*hankel_Db;
    float dig_nb(dig_n*max(max(abs(Rb))));
    digwb=digwb*dig_nb;
    filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=0;j<nwy;j++)
    {
        //if(j<halfary)
        if(weightmat(i,j)<1)
        {
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }
}

#endif
