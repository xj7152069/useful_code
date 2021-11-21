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

struct wiener3d
{
    int nz,nx,ny,nf,nf1,nf2,nwx,nwy,narx,nary,\
        halfnarx,halfnary,nmovex,nmovey,ncpu,covermax;
    float dz,dy,dx,df,dig_n,dig_n2,normal_con;
    fcube data,realrebuildtx,dataagc;
    cx_fcube rebuildfx,rebuildtx,datafx;
    int use_way,agcwx,agcwy,agcwz;
    bool use_anti_alias,normal_bool;
};

//设置一组常用的初始化参数
//Mat(nx,ny,nz); Mat(nx,ny,nf)
void wiener3d_parset(int nx, int ny, int nz,\
    struct wiener3d & par)
{
    par.agcwx=2,par.agcwy=2,par.agcwz=5;
    par.normal_bool=false,par.normal_con=1.0;
    par.use_way=1;
    par.nx=nx,par.ny=ny,par.nz=nz;
    par.dx=10,par.dy=10,par.dz=0.001;
    par.dig_n=0.001;par.dig_n2=0.00;
    par.nwx=1+nx/2;par.nwy=1+ny/2;          
    par.narx=par.nwx/2;par.nary=par.nwy/2;
    if(par.nwx>21){par.nwx=21;}
    if(par.nwy>21){par.nwy=21;}
    if(par.narx>5){par.narx=5;}
    if(par.nary>3){par.nary=3;}
    par.nmovex=par.nwx;par.nmovey=par.nwy;
    par.covermax=1;
    par.ncpu=1;par.nf=1;
    //par.df=1.0/par.dz/par.nf;
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
    if(par.nx<par.ny)
    {
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

void wiener3d_parupdate(struct wiener3d & par)
{
    par.halfnarx=int(par.narx/2);
    par.narx=2*par.halfnarx+1;
    par.halfnary=int(par.nary/2);
    par.nary=2*par.halfnary+1;
    par.data.zeros(par.nx,par.ny,par.nz);
    par.dataagc.zeros(par.nx,par.ny,par.nz);
    par.dataagc.fill(1.0);
    //par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    par.rebuildtx.zeros(par.nx,par.ny,par.nz);
    par.realrebuildtx.zeros(par.nx,par.ny,par.nz);
    par.datafx.zeros(par.nx,par.ny,par.nf);
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
    if(par.ncpu>par.nx)
    {
        par.ncpu=par.nx;
    }
}

void wiener3d_cleardata(struct wiener3d & par)
{
    par.datafx.zeros(par.nx,par.ny,par.nf);
    par.rebuildfx.zeros(par.nx,par.ny,par.nf);
}

//void getagc_wiener3d_thread(struct wiener3d & par)
void getagc_wiener3d(struct wiener3d & par)
{
    int i,j,k,i1,j1,k1;
    float agcpow;
    for(i=par.agcwx;i<par.nx-par.agcwx;i++)
    {
    for(j=par.agcwy;j<par.ny-par.agcwy;j++)
    {
    for(k=par.agcwz;k<par.nz-par.agcwz;k++)
    {
        agcpow=0;
        for(i1=-par.agcwx;i1<=par.agcwx;i1++)
        {
        for(j1=-par.agcwy;j1<=par.agcwy;j1++)
        {
        for(k1=-par.agcwz;k1<=par.agcwz;k1++)
        {
            agcpow+=abs(par.data(i+i1,j+j1,k+k1));
        }}}
        par.dataagc(i,j,k)=agcpow;
        agcpow=0;
    }}}
    for(i=0;i<par.agcwx;i++)
    {
        par.dataagc.row(i)=par.dataagc.row(par.agcwx);
        par.dataagc.row(par.nx-1-i)=par.dataagc.row(par.nx-1-par.agcwx);
    }
    for(i=0;i<par.agcwy;i++)
    {
        par.dataagc.col(i)=par.dataagc.col(par.agcwy);
        par.dataagc.col(par.ny-1-i)=par.dataagc.col(par.ny-1-par.agcwy);
    }
    for(i=0;i<par.agcwz;i++)
    {
        par.dataagc.slice(i)=par.dataagc.slice(par.agcwz);
        par.dataagc.slice(par.nz-1-i)=par.dataagc.slice(par.nz-1-par.agcwz);
    }
    par.dataagc=par.dataagc/(par.dataagc.max()-par.dataagc.min());
    cout<<par.dataagc.max()<<"|"<<par.dataagc.min()<<endl;
    par.dataagc=par.dataagc+0.001*(par.dataagc.max()-par.dataagc.min());
    par.dataagc=1.0/par.dataagc;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
    for(k=0;k<par.nz;k++)
    {
        par.data(i,j,k)=par.dataagc(i,j,k)*par.data(i,j,k);
    }}}
}
void deagc_wiener3d(struct wiener3d & par)
{
    int i,j,k;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
    for(k=0;k<par.nz;k++)
    {
    par.realrebuildtx(i,j,k)=par.realrebuildtx(i,j,k)/par.dataagc(i,j,k);
    par.data(i,j,k)=par.data(i,j,k)/par.dataagc(i,j,k);
    }}}
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

void tx_to_fx3d_wiener3d_pthread(int i, struct wiener3d * par)
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
void tx_to_fx3d_wiener3d_thread(struct wiener3d & par)
{
    int i,k;
    thread *pcal;
    pcal=new thread[par.ncpu];
    
    for(i=0;i<(par.nx-par.ncpu);i+=par.ncpu)
    {
        for(k=0;k<par.ncpu;k++)
        {
            pcal[k]=thread(tx_to_fx3d_wiener3d_pthread,i+k,&par);
        }
        for(k=0;k<par.ncpu;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(par.nx-par.ncpu);i<(par.nx);i++)
    {
        k=i-(par.nx-par.ncpu);
        pcal[k]=thread(tx_to_fx3d_wiener3d_pthread,i,&par);
    }
    for(i=(par.nx-par.ncpu);i<(par.nx);i++)
    {
        k=i-(par.nx-par.ncpu);
        pcal[k].join();
    }
}

void fx_to_tx3d_wiener3d_pthread(int i, struct wiener3d * par)
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

void fx_to_tx3d_wiener3d_thread(struct wiener3d & par)
{
    int i,k;
    thread *pcal;
    pcal=new thread[par.ncpu];
    
    for(i=0;i<(par.nx-par.ncpu);i+=par.ncpu)
    {
        for(k=0;k<par.ncpu;k++)
        {
            pcal[k]=thread(fx_to_tx3d_wiener3d_pthread,i+k,&par);
        }
        for(k=0;k<par.ncpu;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(par.nx-par.ncpu);i<(par.nx);i++)
    {
        k=i-(par.nx-par.ncpu);
        pcal[k]=thread(fx_to_tx3d_wiener3d_pthread,i,&par);
    }
    for(i=(par.nx-par.ncpu);i<(par.nx);i++)
    {
        k=i-(par.nx-par.ncpu);
        pcal[k].join();
    }
}

void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    ofstream & outfarreal,ofstream & outfarimag,\
    ofstream & outfarr, int & outar_flag);
void filter_forword_back(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat);
void filter_xboundary(int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat);
void filter_yboundary(int kwinx, int kf,\
    struct wiener3d & par, fmat & weightmat);
void filter_xpoint(int kf,\
    struct wiener3d & par, fmat & weightmat);
//////////////////////////////////////////////
void wiener3d_mid(struct wiener3d & par, int midnumcpu=1)
{
    struct wiener3d * ppar;
    ppar=&par;
    //tx_to_fx3d_wiener3d(ppar[0]);
    tx_to_fx3d_wiener3d_thread(ppar[0]);

    int i,j,k;//cout<<"ok"<<endl;
    int nx(par.nx),nz(par.nz),ny(par.ny),nf(par.nf),\
        halfarx(par.halfnarx),halfary(par.halfnary),\
        nwx(par.nwx),nwy(par.nwy),movex(par.nmovex),\
        movey(par.nmovey),narx(par.narx),nary(par.nary);
    int nar(narx*nary-1),nw(nwx*nwy);
    fmat ar_weight(nx,ny);

    int kwinx,kwiny,kf,kout(0);
    ofstream outfarreal,outfarimag,outfr;

    for(kf=par.nf1;kf<par.nf2;kf++)
    {
        ar_weight.fill(0.0);kout=0;
        //cout<<"now is run kf = "<<kf<<endl;
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex)
        { 
            for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
            {
                filter_mid(kwinx,kwiny,kf,ppar[0],ar_weight,\
                    outfarreal,outfarimag,outfr,kout);
                //filter_forword_back(kwinx,kwiny,kf,ppar[0],ar_weight);
            }
            filter_mid(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight,\
                outfarreal,outfarimag,outfr,kout);
            //filter_forword_back(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight);
        }
        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
        {
            filter_mid(nx-halfarx-nwx,kwiny,kf,ppar[0],ar_weight,\
                outfarreal,outfarimag,outfr,kout);
            //filter_forword_back(nx-halfarx-nwx,kwiny,kf,ppar[0],ar_weight);
        }
        filter_mid(nx-halfarx-nwx,ny-halfary-nwy,kf,ppar[0],ar_weight,\
            outfarreal,outfarimag,outfr,kout);
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
        filter_xpoint(kf,ppar[0],ar_weight);

        for(i=0;i<nx;i++)
        {
        for(j=0;j<ny;j++)  
        {
            par.rebuildfx(i,j,kf).imag(2*imag(par.rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000000001));
            par.rebuildfx(i,j,kf).real(2*real(par.rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000000001));
        }
        }
    }

    //fx_to_tx3d_wiener3d(ppar[0]);
    fx_to_tx3d_wiener3d_thread(ppar[0]);
}

void filter_mid(int kwinx, int kwiny, int kf,\
    struct wiener3d & par, fmat & weightmat,\
    ofstream & outfarreal,ofstream & outfarimag,\
    ofstream & outfarr,int & outar_flag)
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
    struct wiener3d * ppar;ppar=&par;

    digw.diag()+=1; digw2.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);
    
    hankel_D.zeros(nw,nar),ar_S.zeros(nw,1),\
    filter_A.zeros(nar,1),R.zeros(nar,nar);
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
    cx_float cx_dig_n;
    cx_dig_n.real(1);
    cx_dig_n.imag(1);
    
    filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    //filter_A=solve(R+digw,hankel_D.t()*ar_S);
    ar_S=hankel_D*filter_A;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        if(weightmat(i,j)<par.covermax)
        {
        par.rebuildfx(i,j,kf)+=ar_S(kw,0);
        weightmat(i,j)+=1;
        }
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
        if(weightmat(i,j)<par.covermax)
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
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {   
        if(weightmat(i,j)<par.covermax)
        {
        par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
        weightmat(i,j)+=1;
        }
        kw++;
    }
    }
}

void filter_xpoint(int kf,\
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

    fmat digwf(narf,narf,fill::zeros),digw2f(narf,narf,fill::zeros);
    digwf.diag()+=1;digw2f.diag(par.dig_n2); 
    datafx=par.datafx.slice(kf);

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=ny-nwy;j<ny;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++)
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
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=ny-nwy;j<ny;j++) 
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++)
        {
        for(kay=j;kay<=j+2*halfary;kay++)  
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
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++)
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
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++) 
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++)
        {
        for(kay=j;kay<=j+2*halfary;kay++)  
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
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
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
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
     
    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
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
    //filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;

    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        //if(i<halfarx)
        if(weightmat(i,j)<par.covermax)
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
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++)
    {
        //if(j>=ny-halfary)
        if(weightmat(i,j)<par.covermax)
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
    //filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;

    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=0;j<nwy;j++)
    {
        //if(j<halfary)
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
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
    ofstream & outfarreal,ofstream & outfarimag,\
    ofstream & outfarr,int outar_flag);
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
    struct wienerCG3d * cg, int kcpu, \
    ofstream * outfarreal, ofstream * outfarimag, ofstream * outfr)
{
    int i,j,k;//cout<<"ok"<<endl;
    int nx(par[0].nx),nz(par[0].nz),ny(par[0].ny),nf(par[0].nf),\
        halfarx(par[0].halfnarx),halfary(par[0].halfnary),\
        nwx(par[0].nwx),nwy(par[0].nwy),movex(par[0].nmovex),\
        movey(par[0].nmovey),narx(par[0].narx),nary(par[0].nary);
    int nar(narx*nary),nw(nwx*nwy);
    fmat ar_weight(nx,ny);
    int kwinx,kwiny,kf,kout(0);

    for(kf=pnf1;kf<pnf2;kf++)  
    {
        ar_weight.fill(0.0);kout=0;
        //cout<<"now is run kf = "<<kf<<endl;
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex)
        {
            for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
            {
                filter_mid(kwinx,kwiny,kf,par[0],ar_weight,cg[0],kcpu,\
                    outfarreal[0],outfarimag[0],outfr[0],kout);
                //filter_forword_back(kwinx,kwiny,kf,ppar[0],ar_weight);
                kout++;
            }
            filter_mid(kwinx,ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu,\
                outfarreal[0],outfarimag[0],outfr[0],kout);
            kout++;
            //filter_forword_back(kwinx,ny-halfary-nwy,kf,ppar[0],ar_weight);
        }
        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
        {
            filter_mid(nx-halfarx-nwx,kwiny,kf,par[0],ar_weight,cg[0],kcpu,\
                outfarreal[0],outfarimag[0],outfr[0],kout);
            //filter_forword_back(nx-halfarx-nwx,kwiny,kf,ppar[0],ar_weight);
        }
        filter_mid(nx-halfarx-nwx,ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu,\
            outfarreal[0],outfarimag[0],outfr[0],kout);
        //filter_forword_back(nx-halfarx-nwx,ny-halfary-nwy,kf,ppar[0],ar_weight);

        for(kwiny=halfary;kwiny<ny-halfary-nwy;kwiny+=movey)
        {
            filter_xboundary(kwiny,kf,par[0],ar_weight,cg[0],kcpu);
        }
        filter_xboundary(ny-halfary-nwy,kf,par[0],ar_weight,cg[0],kcpu);
        
        for(kwinx=halfarx;kwinx<nx-halfarx-nwx;kwinx+=movex)
        {
            filter_yboundary(kwinx,kf,par[0],ar_weight,cg[0],kcpu);
        }
        filter_yboundary(nx-halfarx-nwx,kf,par[0],ar_weight,cg[0],kcpu);
        filter_xpoint(kf,par[0],ar_weight,cg[0],kcpu);

        for(i=0;i<nx;i++)
        {
        for(j=0;j<ny;j++)  
        {
            par[0].rebuildfx(i,j,kf).imag(2*imag(par[0].rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000000001));
            par[0].rebuildfx(i,j,kf).real(2*real(par[0].rebuildfx(i,j,kf))\
                /(ar_weight(i,j)+0.000000001));
        }
        }
    }
}

void wiener3d_mid_thread(struct wiener3d & par,int midnumcpu, int numi=25, int mini=9)
{ 
    struct wiener3d * ppar;
    struct wienerCG3d cg;
    ppar=&par;
    //getagc_wiener3d(ppar[0]);
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
    ofstream outfarreal,outfarimag,outfr;

///////////////////////////////////////////

    for(k=0;k<pn-1;k++)
    {
        pnf1=round(par.nf1+k*dnf);
        pnf2=round(par.nf1+(k+1)*dnf);
        pcal[k]=thread(wiener3d_mid_pthread,&par,pnf1,pnf2,\
            &cg, k, &outfarreal, &outfarimag, &outfr);
    }
    k=pn-1;
    pnf1=round(par.nf1+k*dnf);
    pnf2=par.nf2;
    pcal[k]=thread(wiener3d_mid_pthread,&par,pnf1,pnf2,\
        &cg, k, &outfarreal, &outfarimag, &outfr);

    for(k=0;k<pn;k++)
    {
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
    for(i=0;i<ncpu;i++)
    {
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

for(k=0;k<cg.numi;k++)
{ 
    memat = A*cg.p_k1_cxfmatp1ncpunpxy[kcpu][0];
    me = cg.p_k1_cxfmatp1ncpunpxy[kcpu][0].st()*(memat.t()).st();
    ak = cg.r_k1_cxfmatp1ncpunpxy[kcpu][0].t()*\
        cg.r_k1_cxfmatp1ncpunpxy[kcpu][0];
    //cout<<"|"<<abs(me(0,0));
    ak(0,0)= ak(0,0)/me(0,0);
    if(k>cg.mini && abs(me(0,0))<errcg)
    {
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
    if(abs(me(0,0))<abs(meo(0,0)))
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
        X=cg.x_k1_cxfmatp1ncpunpxy[kcpu][0];
    }
    else
    {
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
    ofstream & outfarreal,ofstream & outfarimag,\
    ofstream & outfarr,int outar_flag)
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
 
    //filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    //void wienerGC3d_getfilter(cx_fmat & X, cx_fmat A, cx_fmat B, \
    struct wienerCG3d & cg, int kcpu)
    if(par.use_way==1)
    {
        filter_A=inv(R+digw+digw2)*hankel_D.t()*ar_S;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_A.zeros(nar,1),pcg=&cg;
    wienerGC3d_getfilter(filter_A, R+digw+digw2,\
        hankel_D.t()*ar_S, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_A=solvecxfmat((R+digw+digw2),hankel_D.t()*ar_S,nar);
    }
    else if(par.use_way==4)
    {
        filter_A=invcxfmat(R+digw+digw2,nar)*hankel_D.t()*ar_S;
    }

    //filter_A=solve(R+digw,hankel_D.t()*ar_S);
    ar_S=hankel_D*filter_A;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        if(weightmat(i,j)<par.covermax)
        {
        par.rebuildfx(i,j,kf)+=ar_S(kw,0);
        weightmat(i,j)+=1;
        }
        kw++;
    }
    }
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
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1)
    {
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Af.zeros(narf,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
        hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4)
    {
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
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
    //filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    if(par.use_way==1)
    {
        filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Ab.zeros(narb,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Ab, Rb+digwb+digw2b,\
        hankel_Db.t()*ar_Sb, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Ab=solvecxfmat((Rb+digwb+digw2b),\
            hankel_Db.t()*ar_Sb,narb);
    }
    else if(par.use_way==4)
    {
        filter_Ab=invcxfmat(Rb+digwb+digw2b,narb)\
            *hankel_Db.t()*ar_Sb;
    }

    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=kwiny;j<kwiny+nwy;j++)  
    {
        //if(i<halfarx)
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }
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
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1)
    {
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Af.zeros(narf,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
        hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4)
    {
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++)
    {
        //if(j>=ny-halfary)
        if(weightmat(i,j)<par.covermax)
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
    //filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    //filter_Ab=solve(Rb+digwb,hankel_Db.t()*ar_Sb);
    if(par.use_way==1)
    {
        filter_Ab=inv(Rb+digwb+digw2b)*hankel_Db.t()*ar_Sb;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Ab.zeros(narb,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Ab, Rb+digwb+digw2b,\
        hankel_Db.t()*ar_Sb, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Ab=solvecxfmat((Rb+digwb+digw2b),\
            hankel_Db.t()*ar_Sb,narb);
    }
    else if(par.use_way==4)
    {
        filter_Ab=invcxfmat(Rb+digwb+digw2b,narb)\
            *hankel_Db.t()*ar_Sb;
    }

    ar_Sb=hankel_Db*filter_Ab;
    kw=0;
    for(i=kwinx;i<kwinx+nwx;i++)
    {
    for(j=0;j<nwy;j++)
    {
        //if(j<halfary)
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sb(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }
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
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=ny-nwy;j<ny;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++)
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
    Rf=hankel_Df.t()*hankel_Df;
    float dig_nf(dig_n*max(max(abs(Rf))));
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1)
    {
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Af.zeros(narf,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
        hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4)
    {
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }
    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=ny-nwy;j<ny;j++) 
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i-2*halfarx;kax<=i;kax++)
        {
        for(kay=j;kay<=j+2*halfary;kay++)  
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
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1)
    {
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Af.zeros(narf,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
        hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4)
    {
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=nx-nwx;i<nx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++)
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
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1)
    {
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Af.zeros(narf,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
        hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4)
    {
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=ny-nwy;j<ny;j++) 
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }

    hankel_Df.zeros(nw,narf),ar_Sf.zeros(nw,1),\
    filter_Af.zeros(narf,1),Rf.zeros(narf,narf);
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        ar_Sf(kw,0)=datafx(i,j);
        kar=0;
        for(kax=i;kax<=i+2*halfarx;kax++)
        {
        for(kay=j;kay<=j+2*halfary;kay++)  
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
    Rf=hankel_Df.t()*hankel_Df;
    dig_nf=(dig_n*max(max(abs(Rf))));
    digwf.fill(0.0);
    digwf.diag()+=1;
    digwf=digwf*dig_nf;
    //filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    //filter_Af=solve(Rf+digwf,hankel_Df.t()*ar_Sf);
    if(par.use_way==1)
    {
        filter_Af=inv(Rf+digwf+digw2f)*hankel_Df.t()*ar_Sf;
    }
    else if(par.use_way==2)
    {
    struct wienerCG3d * pcg;
    filter_Af.zeros(narf,1),pcg=&cg;
    wienerGC3d_getfilter(filter_Af, Rf+digwf+digw2f,\
        hankel_Df.t()*ar_Sf, pcg[0], kcpu);
    }
    else if(par.use_way==3)
    {
        filter_Af=solvecxfmat((Rf+digwf+digw2f),\
            hankel_Df.t()*ar_Sf,narf);
    }
    else if(par.use_way==4)
    {
        filter_Af=invcxfmat(Rf+digwf+digw2f,narf)\
            *hankel_Df.t()*ar_Sf;
    }

    ar_Sf=hankel_Df*filter_Af;
    kw=0;
    for(i=0;i<nwx;i++)
    {
    for(j=0;j<nwy;j++)  
    {
        //if(i>=(nx-halfarx))
        if(weightmat(i,j)<par.covermax)
        {
            par.rebuildfx(i,j,kf)+=ar_Sf(kw,0);
            weightmat(i,j)+=1;
        }
        kw++;
    }
    }
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
    for(k2=0;k2<n;k2++)
    {
        for(k1=1+k2;k1<n;k1++)
        {
            if(abs(mat1(k1,k2))>0.000001)
            {
                xs=mat1(k2,k2)/mat1(k1,k2);
                mat1.row(k1)=xs*mat1.row(k1)-mat1.row(k2);
                mat2.row(k1)=xs*mat2.row(k1)-mat2.row(k2);
            }
        }
    }
    for(k2=0;k2<n;k2++)
    {
        for(k1=0;k1<k2;k1++)
        {
            if(abs(mat1(k1,k2))>0.000001)
            {
                xs=mat1(k2,k2)/mat1(k1,k2);
                mat1.row(k1)=xs*mat1.row(k1)-mat1.row(k2);
                mat2.row(k1)=xs*mat2.row(k1)-mat2.row(k2);
            }
        }
    }
    for(k2=0;k2<n;k2++)
    {
        mat2.row(k2)=mat2.row(k2)/mat1(k2,k2);
    }
    mat2=mat2/max_mat_elem;
    //cout<<mat2.has_inf()<<","<<mat2.has_nan()<<"|";
    
    for(k2=0;k2<n;k2++)
    {
        for(k1=0;k1<n;k1++)
        {
            if(isnan(mat2(k1,k2)))
                {mat2(k1,k2)=0;}
            if(isinf(mat2(k1,k2)))
                {mat2(k1,k2)=0;}
        }
    }
    
    return mat2;
}

inline cx_fmat solvecxfmat(cx_fmat mata, cx_fmat matb, int n)
{
    cx_fmat mat2(n,1,fill::zeros);
    int k1,k2;
    cx_float xs;
    float max_mata_ele;
    max_mata_ele=0.000001*abs(mata).max();

    for(k2=0;k2<n;k2++)
    {
        for(k1=1+k2;k1<n;k1++)
        {
            if(abs(mata(k1,k2))>max_mata_ele)
            {
                xs=mata(k2,k2)/mata(k1,k2);
                mata.row(k1)=xs*mata.row(k1)-mata.row(k2);
                matb.row(k1)=xs*matb.row(k1)-matb.row(k2);
            }
        }
    }

    for(k2=0;k2<n;k2++)
    {
        if(isnan(abs(matb(k2,0))))
        {
            matb(k2,0).imag(0);
            matb(k2,0).real(0);
        }
        if(isinf(abs(matb(k2,0))))
        {
            matb(k2,0).imag(0);
            matb(k2,0).real(0);
        }
        for(k1=0;k1<n;k1++)
        {
            if(isnan(abs(mata(k1,k2))))
            {
                mata(k1,k2).imag(0);
                mata(k1,k2).real(0);
            }
            if(isinf(abs(mata(k1,k2))))
            {
                mata(k1,k2).imag(0);
                mata(k1,k2).real(0);
            }
        }
    }

    for(k2=n-1;k2>=0;k2--)
    {
        xs=matb(k2,0);
        for(k1=n-1;k1>k2;k1--)
        {
            xs=xs-mat2(k1,0)*mata(k2,k1);
        }
        mat2(k2,0)=xs/mata(k2,k2);
    }

    for(k2=0;k2<n;k2++)
    {
        if(isnan(abs(mat2(k2,0))))
        {
            mat2(k2,0).imag(0);
            mat2(k2,0).real(0);
        }
        if(isinf(abs(mat2(k2,0))))
        {
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

    for(k2=0;k2<n;k2++)
    {
        for(k1=1+k2;k1<n;k1++)
        {
            if(abs(mata(k1,k2))>max_mata_ele)
            {
                xs=mata(k2,k2)/mata(k1,k2);
                mata.row(k1)=xs*mata.row(k1)-mata.row(k2);
                matb.row(k1)=xs*matb.row(k1)-matb.row(k2);
            }
        }
    }

    for(k2=0;k2<n;k2++)
    {
        if(isnan(abs(matb(k2,0))))
        {
            matb(k2,0)=(0);
        }
        if(isinf(abs(matb(k2,0))))
        {
            matb(k2,0)=(0);
        }
        for(k1=0;k1<n;k1++)
        {
            if(isnan(abs(mata(k1,k2))))
            {
                mata(k1,k2)=(0);
            }
            if(isinf(abs(mata(k1,k2))))
            {
                mata(k1,k2)=(0);
            }
        }
    }

    for(k2=n-1;k2>=0;k2--)
    {
        xs=matb(k2,0);
        for(k1=n-1;k1>k2;k1--)
        {
            xs=xs-mat2(k1,0)*mata(k2,k1);
        }
        mat2(k2,0)=xs/mata(k2,k2);
    }

    for(k2=0;k2<n;k2++)
    {
        if(isnan(abs(mat2(k2,0))))
        {
            mat2(k2,0)=(0);
        }
        if(isinf(abs(mat2(k2,0))))
        {
            mat2(k2,0)=(0);
        }
    }
    return mat2;
}

#endif
