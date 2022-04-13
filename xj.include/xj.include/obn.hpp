#ifndef OBN_HPP
#define OBN_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <thread>
#include <future>
#include <time.h>
#include "../xjc.h"
using namespace std;
using namespace arma;

/////////////////////////////////////////////////////
struct demultiple2d
{
    int nx,ny;
    float dy,dx,dkx,dky;
    int agcwx,agcwy;
    fmat datap,datavz,datapagc,datavzagc;
    fmat datap2,datavz2;
    fmat dataup,datadown;

    fmat data2d, datadict, datah, datad, datahd;
    int&n1=nx;
    int&n2=ny;
    int nwp,nwph,nwpd,nwphd;
    float lmd;
    int n1w,d1w;
};
////////////////////////////////////////////////////////////////
void single_trace_dewave(fmat & dataup, fmat & datadown,\
    fmat & data, fmat & dict,\
    fmat &dataph, fmat &datapd, fmat &dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh, int ncpu);
void single_trace_dewave(struct demultiple2d & par, int ncpu=1)
{
    single_trace_dewave(par.dataup, par.datadown,\
    par.data2d, par.datadict,\
    par.datah, par.datad, par.datahd,\
    par.n1,par.n2,par.nwp,par.nwph,par.nwpd,par.nwphd,\
    par.n1w,par.d1w,par.lmd,ncpu);
}

void single_trace_dewave_pthread(fmat* dataup, fmat* datadown,\
    fmat* datavz, fmat* datap,\
    fmat*dataph, fmat* datapd, fmat*dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh,int nx1,int nx2);
void single_trace_dewave(fmat & dataup, fmat & datadown,\
    fmat & data, fmat & dict,\
    fmat &dataph, fmat &datapd, fmat &dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh, int ncpu)
{
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd),nx1,nx2,dnx;
    //fmat mat1(nwt,nw),mat0(nwt,nw),data2(nt,nx,fill::zeros),\
    //    matq(nw,1),mat2(nw,1),matd(nwt,1),matD(nw,nw),matq0(nw,1);
    //fmat digmat(nw,nw,fill::zeros);
    fmat datap(nt,nx),datavz(nt,nx);

////////////////////////////////////////////////////////////////////
    float xs;
    xs=data.max()-data.min();
    datavz=data/(xs);
    datap=dict/(dict.max()-dict.min());
    dataph=dataph/(dataph.max()-dataph.min());
    datapd=datapd/(datapd.max()-datapd.min());
    dataphd=dataphd/(dataphd.max()-dataphd.min());
////////////////////////////////////////////////////////////////////
    fmat* pdataup=(&dataup); fmat* pdatadown=(&datadown);
    fmat* pdata=(&datavz); fmat* pdict=(&datap);
    fmat* pdataph=(&dataph); fmat* pdatapd=(&datapd);
    fmat* pdataphd=(&dataphd);
    thread* pcal;
    pcal=new thread[ncpu];
    dnx=nx/ncpu;
    nx1=0;
    for(k=0;k<ncpu-1;k++){
        nx2=nx1+dnx;
        pcal[k]=thread(single_trace_dewave_pthread,pdataup, pdatadown,\
            pdata, pdict,pdataph, pdatapd, pdataphd,\
            nt, nx,nwp, nwph, nwpd, nwphd,nwt, dnwt,zzh,nx1,nx2);
        nx1=nx2;
    }
    pcal[ncpu-1]=thread(single_trace_dewave_pthread,pdataup, pdatadown,\
        pdata, pdict,pdataph, pdatapd, pdataphd,\
        nt, nx,nwp, nwph, nwpd, nwphd,nwt, dnwt,zzh,nx1,nx);
    for(k=0;k<ncpu;k++){
        pcal[k].join();
    }

/*
    for(k=0;k<nx;k++){
        cout<<"tracl = "<<k<<endl;
        for(kwt=0;kwt<=(nt-nwt);kwt+=dnwt){
            mat1.fill(0.0);

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwp;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=datap(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwp;j++){
                mat1.col(j-kwt)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwph;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=dataph(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwph;j++){
                mat1.col(j-kwt+nwp)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwpd;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=datapd(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwpd;j++){
                mat1.col(j-kwt+nwp+nwph)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwphd;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=dataphd(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwphd;j++){
                mat1.col(j-kwt+nwp+nwph+nwpd)=mat0.col(j-kwt);
            }
            ///////////
            
            for(i=0;i<nwt;i++){
                matd(i,0)=datavz(i+kwt,k);
            }
            matD=mat1.t()*mat1;
            
            digmat.fill(0);
            digmat.diag()+=(matD.max()*zzh+zzh);
            //cout<<digmat(1,1)<<endl;
            if(kwt==0){
                matq=inv(matD+digmat)*mat1.t()*matd;
                matq0=matq;
            }
            else{
                matq=matq0;
            }
            matd=mat1*matq;
            for(i=nw;i<nwt;i++){
                dataup(i+kwt,k)+=datavz(i+kwt,k)-matd(i,0);
                datadown(i+kwt,k)+=datavz(i+kwt,k)+matd(i,0);
            }
        }
    }
*/
////////////////////////////////////////////////////////////////////
    for(i=0;i<nt;i++){
        for(j=0;j<nx;j++){
            if(isnan(abs(dataup(i,j))))
                {dataup(i,j)=0;}
            if(isnan(abs(datadown(i,j))))
                {datadown(i,j)=0;}
        }
    }
    datadown*=xs;
    dataup*=xs;
}

void single_trace_dewave_pthread(fmat* dataup, fmat* datadown,\
    fmat* datavz, fmat* datap,\
    fmat*dataph, fmat* datapd, fmat*dataphd,\
     int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh,int nx1,int nx2)
{
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd);
    fmat mat1(nwt,nw),mat0(nwt,nw),data2(nt,nx,fill::zeros),\
        matq(nw,1),mat2(nw,1),matd(nwt,1),matD(nw,nw),matq0(nw,1);
    fmat digmat(nw,nw,fill::zeros);
    digmat.diag()+=1;

        for(k=nx1;k<nx2;k++){
        cout<<"tracl = "<<k<<endl;
        for(kwt=0;kwt<=(nt-nwt);kwt+=dnwt){
            mat1.fill(0.0);

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwp;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=datap[0](i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwp;j++){
                mat1.col(j-kwt)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwph;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=dataph[0](i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwph;j++){
                mat1.col(j-kwt+nwp)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwpd;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=datapd[0](i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwpd;j++){
                mat1.col(j-kwt+nwp+nwph)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwphd;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=dataphd[0](i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwphd;j++){
                mat1.col(j-kwt+nwp+nwph+nwpd)=mat0.col(j-kwt);
            }
            ///////////
            
            for(i=0;i<nwt;i++){
                matd(i,0)=datavz[0](i+kwt,k);
            }
            matD=mat1.t()*mat1;
            
            digmat.fill(0);
            digmat.diag()+=(matD.max()*zzh+zzh);
            //cout<<digmat(1,1)<<endl;
            if(kwt==0){
                matq=inv(matD+digmat)*mat1.t()*matd;
                matq0=matq;
            }
            else{
                matq=matq0;
            }
            matd=mat1*matq;
            for(i=nw;i<nwt;i++){
                dataup[0](i+kwt,k)+=datavz[0](i+kwt,k)-matd(i,0);
                datadown[0](i+kwt,k)+=datavz[0](i+kwt,k)+matd(i,0);
            }
        }
    }

}
//////////////////////////////////////////////////////////////////////
void multiple_code2d(fmat& u2, fmat& u1, fmat& green,\
 float df, int fn1=0, int fn2=2048, int ncpu=1);
void multiple_code2d_pthread(cx_fmat* u2cx, cx_fmat* u1cx, fmat* green,\
 float df, int fn1, int fn2);
void multiple_code2d(fmat& u2, fmat& u1, fmat& green,\
 float df, int fn1, int fn2, int ncpu)
{
    int n1(u1.n_rows),n2(u1.n_cols);
    int i,j,k,dfn,pfn1,pfn2;
    float w,pi(3.1415926);
    cx_fmat codemat(n2,n2,fill::zeros),onecol(n1,1),\
        u1cx(n2,n1),u2cx(n2,n1,fill::zeros);
    cx_float a;
    a.real(0.0);

    for(k=0;k<n2;k++){
        onecol.col(0)=fft(u1.col(k),n1);
        u1cx.row(k)=onecol.col(0).st();
        //u2cx.row(k)=fft(u2.col(k).t(),n1);
    }

    dfn=(fn2-fn1)/ncpu;
    pfn1=fn1;
    thread * pcal;
    pcal=new thread[ncpu];
    cx_fmat* pu2cx=(&u2cx);
    cx_fmat* pu1cx=(&u1cx);
    fmat* pgreen=(&green);
    for(k=0;k<ncpu-1;k++){
        pfn2=pfn1+dfn;
        pcal[k]=thread(multiple_code2d_pthread,pu2cx,pu1cx,pgreen,df,pfn1,pfn2);
        pfn1=pfn2;
    }
    pcal[ncpu-1]=thread(multiple_code2d_pthread,pu2cx,pu1cx,pgreen,df,pfn1,fn2);
    for(k=0;k<ncpu;k++){
        pcal[k].join();
    }
/*
    for(k=fn1;k<fn2;k++){
        w=-2.0*pi*df*k;
        if(k%100==0){
            cout<<"k="<<k<<" | "<<"w="<<w<<endl;
        }
        for(i=0;i<n2;i++){
            for(j=0;j<n2;j++){
                a.imag(w*green(i,j));
                codemat(i,j)=exp(a);
            }
        }
        u2cx.col(k)=codemat*u1cx.col(k);
    }
*/
    for(k=0;k<n2;k++){
        //u1cx.row(k)=fft(u1.col(k).t(),n1);
        onecol.col(0)=u2cx.row(k).st();
        u2.col(k)=real(ifft(onecol.col(0),n1));
    }
}

void multiple_code2d_pthread(cx_fmat* u2cx, cx_fmat* u1cx, fmat* green,\
 float df, int fn1, int fn2)
{
    
    int i,j,k;
    float w,pi(3.1415926);
    cx_float a;
    a.real(0.0);
    int n2(u1cx->n_rows);
    cx_fmat codemat(n2,n2,fill::zeros);

    for(k=fn1;k<fn2;k++){
        w=-2.0*pi*df*k;
        if(k%100==0){
            cout<<"k="<<k<<" | "<<"w="<<w<<endl;
        }
        for(i=0;i<n2;i++){
            for(j=0;j<n2;j++){
                a.imag(w*green[0](i,j));
                codemat(i,j)=exp(a);
            }
        }
        u2cx[0].col(k)=codemat*u1cx[0].col(k);
    }
    
}


void getagc_demultiple2d(struct demultiple2d & par, fmat & data,\
     fmat & dataagc, float par1)
{
    int i,j,i1,j1;
    float agcpow,maxpow(1.0/par1);
    for(i=par.agcwx;i<par.nx-par.agcwx;i++)
    {
    for(j=par.agcwy;j<par.ny-par.agcwy;j++)
    {
        agcpow=0;
        for(i1=-par.agcwx;i1<=par.agcwx;i1++)
        {
        for(j1=-par.agcwy;j1<=par.agcwy;j1++)
        {
            agcpow+=abs(data(i+i1,j+j1));
        }}
        dataagc(i,j)=agcpow;
        agcpow=0;
    }}
    for(i=0;i<par.agcwx;i++)
    {
        dataagc.row(i)=dataagc.row(par.agcwx);
        dataagc.row(par.nx-1-i)=dataagc.row(par.nx-1-par.agcwx);
    }
    for(i=0;i<par.agcwy;i++)
    {
        dataagc.col(i)=dataagc.col(par.agcwy);
        dataagc.col(par.ny-1-i)=dataagc.col(par.ny-1-par.agcwy);
    }

    dataagc=dataagc/(dataagc.max()-dataagc.min());
    cout<<dataagc.max()<<"|"<<dataagc.min()<<endl;
    dataagc=dataagc+maxpow*(dataagc.max()-dataagc.min());
    dataagc=1.0/dataagc;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
        data(i,j)=dataagc(i,j)*data(i,j);
    }}
}

void deagc_demultiple2d(struct demultiple2d & par, fmat & data, fmat & dataagc)
{
    int i,j;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
        data(i,j)=data(i,j)/dataagc(i,j);
    }}
}

void addagc_demultiple2d(struct demultiple2d & par, fmat & data, fmat & dataagc)
{
    int i,j;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
        data(i,j)=data(i,j)*dataagc(i,j);
    }}
}

#endif
