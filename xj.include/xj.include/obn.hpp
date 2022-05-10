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
    int nwt,int dnwt,float zzh, int ncpu, bool multiple_win);
void single_trace_dewave(struct demultiple2d & par, int ncpu=1, bool multiple_win=true)
{
    single_trace_dewave(par.dataup, par.datadown,\
    par.data2d, par.datadict,\
    par.datah, par.datad, par.datahd,\
    par.n1,par.n2,par.nwp,par.nwph,par.nwpd,par.nwphd,\
    par.n1w,par.d1w,par.lmd,ncpu,multiple_win);
}

void single_trace_dewave_pthread(fmat* dataup, fmat* datadown,\
    fmat* datavz, fmat* datap,\
    fmat*dataph, fmat* datapd, fmat*dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh,int nx1,int nx2, bool multiple_win);
void single_trace_dewave(fmat & dataup, fmat & datadown,\
    fmat & data, fmat & dict,\
    fmat &dataph, fmat &datapd, fmat &dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh, int ncpu, bool multiple_win=true)
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
            nt, nx,nwp, nwph, nwpd, nwphd,nwt, dnwt,zzh,nx1,nx2,multiple_win);
        nx1=nx2;
    }
    pcal[ncpu-1]=thread(single_trace_dewave_pthread,pdataup, pdatadown,\
        pdata, pdict,pdataph, pdatapd, pdataphd,\
        nt, nx,nwp, nwph, nwpd, nwphd,nwt, dnwt,zzh,nx1,nx,multiple_win);
    for(k=0;k<ncpu;k++){
        pcal[k].join();
    }
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
    int nwt,int dnwt,float zzh,int nx1,int nx2, bool multiple_win=true)
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
            if(kwt==0 || multiple_win){
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

//////////////////////////////////////////////////////////////////////////////////////
/*
void multiple_code3d_onepoint_onefrequence(cx_fmat* u2, cx_fmat* u1, fmat* green,\
 int i1, int j1, float f1)
{
    int n1(u1[0].n_rows),n2(u1[0].n_cols);
    int i,j;
    float w,pi(3.1415926);
    cx_float a;
    a.real(0.0);

    w=-2.0*pi*f1;
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            a.imag(w*green[0](i,j));
            a=exp(a);
            u2[0](i,j)+=a*u1[0](i1,j1);
        }
    }
}*/
void multiple_code3d_onepoint_onefrequence(cx_fcube* u2, cx_fcube* u1, fmat* green,\
 int i1, int j1, int nf1, float df)
{
    int n1(u1[0].n_rows),n2(u1[0].n_cols);
    int i,j;
    float w,pi(3.1415926);
    cx_float a;
    
    w=-2.0*pi*df*nf1;
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            a.real(0.0);
            a.imag(w*green[0](i,j));
            a=exp(a);
            u2[0](i,j,nf1)+=a*u1[0](i1,j1,nf1);
        }
    }

}
void multiple_code3d_onepoint(cx_fcube& u2, cx_fcube& u1, fmat& green,\
 int i1, int j1, float df, int fn1, int fn2, int ncpu)
{
    int n1(u1.n_rows),n2(u1.n_cols),n3(u1.n_slices);
    int i,j,k;
    float w,pi(3.1415926),f1;

    thread * pcal;
    pcal=new thread[ncpu];
    cx_fcube *pu2=(&u2);
    cx_fcube *pu1=(&u1);
    fmat* pgreen=(&green);

    for(i=fn1;i<(fn2);i+=ncpu){
        for(k=0;k<ncpu;k++){
            pcal[k]=thread(multiple_code3d_onepoint_onefrequence,\
                pu2,pu1, pgreen,i1, j1,(i+k),df);
        }
        for(k=0;k<ncpu;k++){
            pcal[k].join();
        }
    }
    i=i-ncpu;
    for(k=i;k<fn2;k++){
        pcal[k-i]=thread(multiple_code3d_onepoint_onefrequence,\
            pu2,pu1, pgreen,i1, j1,k,df);
    }
    for(k=i;k<fn2;k++){
        pcal[k-i].join();
    }
}

void multiple_code3d(cx_fcube& u2, cx_fcube& u1, fmat& seabase_depth,\
 int source_site_x, int source_site_y, float water_velocity,\
 float dx, float dy, float df, int fn1, int fn2, int ncpu)
{
    int n1(u1.n_rows),n2(u1.n_cols),n3(u1.n_slices);
    int i,j,k,i1,j1,i2,j2,ix,jy,sx(source_site_x),sy(source_site_y);
    float depth, half_offset, d11,d12,d21,d22, l11,l12,l21,l22,fi,fj;
    float l_min(0.000000001),l_sum(0.0),l_trace;
    fmat green(n1,n2);
    cx_fcube *pu2(&u2);
    cx_fcube *pu1(&u1);
for(ix=0;ix<n1;ix++){
    cout<<ix<<endl;
    sx=ix;
for(jy=0;jy<n2;jy++){
    sy=jy;
    for(i=0;i<n1;i++){      
        fi=((float(sx)+i)/2.0);
        i1=floor((sx+i)/2);
        i2=1+floor((sx+i)/2);
        i1=max(i1,0);
        i2=min(i2,n1-1);
        for(j=0;j<n2;j++){
            fj=((float(sy)+j)/2.0);
            j1=floor((sy+j)/2);
            j2=1+floor((sy+j)/2);
            j1=max(j1,0);
            j2=min(j2,n2-1);
            l11=dx*dx*(fi-i1)*(fi-i1)+dy*dy*(fj-j1)*(fj-j1)+l_min;
            l12=dx*dx*(fi-i1)*(fi-i1)+dy*dy*(fj-j2)*(fj-j2)+l_min;
            l21=dx*dx*(fi-i2)*(fi-i2)+dy*dy*(fj-j1)*(fj-j1)+l_min;
            l22=dx*dx*(fi-i2)*(fi-i2)+dy*dy*(fj-j2)*(fj-j2)+l_min;
            l11=1/l11;l12=1/l12;l21=1/l21;l22=1/l22;
            l_sum=l11+l12+l21+l22;
            l11/=l_sum;l12/=l_sum;l21/=l_sum;l22/=l_sum;

            d11=seabase_depth(i1,j1);
            d12=seabase_depth(i1,j2);
            d21=seabase_depth(i2,j1);
            d22=seabase_depth(i2,j2);
            depth=l11*d11+l21*d21+l12*d12+l22*d22;
            half_offset=sqrt(dx*dx*(fi-i)*(fi-i)+dy*dy*(fj-j)*(fj-j));
            l_trace=2*sqrt(half_offset*half_offset+depth*depth);
            green(i,j)=l_trace/water_velocity;
        }
    }
    multiple_code3d_onepoint(pu2[0],pu1[0],green,sx,sy,df,fn1,fn2,ncpu);
}}
}
////////////function: tx2fx or fx2tx/////////////
void tx2fx_3d_pthread(cx_fcube *data3d_out,fcube *data3d,int i)
{
    int j,n1(data3d[0].n_rows),n2(data3d[0].n_cols),n3(data3d[0].n_slices);
    cx_fmat datafx(n2,n3);
    fmat data(n2,n3);
    data=data3d[0].row(i);
    for(j=0;j<n2;j++)
    {
        datafx.row(j)=fft(data.row(j),n3);
    }
    data3d_out[0].row(i)=datafx;
}
void tx2fx_3d_thread(cx_fcube &data3d_out,fcube &data3d, int ncpu)
{
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices);
    cx_fcube *pdata3d_out=&data3d_out;
    fcube *pdata3d=&data3d;
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=0;i<(n1-ncpu);i+=ncpu){
        for(k=0;k<ncpu;k++)
        {
            pcal[k]=thread(tx2fx_3d_pthread,pdata3d_out,pdata3d,i+k);
        }
        for(k=0;k<ncpu;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(n1-ncpu);i<(n1);i++)
    {
        k=i-(n1-ncpu);
        pcal[k]=thread(tx2fx_3d_pthread,pdata3d_out,pdata3d,i);
    }
    for(i=(n1-ncpu);i<(n1);i++)
    {
        k=i-(n1-ncpu);
        pcal[k].join();
    }
}

void fx2tx_3d_pthread(fcube *data3d_out,cx_fcube *data3d,int i)
{
    int j,n1(data3d[0].n_rows),n2(data3d[0].n_cols),n3(data3d[0].n_slices);
    cx_fmat data(n2,n3),datafx(n2,n3);
    datafx=data3d[0].row(i);
    for(j=0;j<n2;j++)
    {
        data.row(j)=ifft(datafx.row(j),n3);
    }
    data3d_out[0].row(i)=real(data);
}
void fx2tx_3d_thread(fcube &data3d_out,cx_fcube &data3d, int ncpu)
{
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices);
    fcube *pdata3d_out=&data3d_out;
    cx_fcube *pdata3d=&data3d;
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=0;i<(n1-ncpu);i+=ncpu){
        for(k=0;k<ncpu;k++)
        {
            pcal[k]=thread(fx2tx_3d_pthread,pdata3d_out,pdata3d,i+k);
        }
        for(k=0;k<ncpu;k++)
        {
            pcal[k].join();
        }
    } 
    for(i=(n1-ncpu);i<(n1);i++)
    {
        k=i-(n1-ncpu);
        pcal[k]=thread(fx2tx_3d_pthread,pdata3d_out,pdata3d,i);
    }
    for(i=(n1-ncpu);i<(n1);i++)
    {
        k=i-(n1-ncpu);
        pcal[k].join();
    }
}


///////////////////////////////////////////////////////////////////////
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
