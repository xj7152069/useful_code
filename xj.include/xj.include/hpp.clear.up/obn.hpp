/*
    Multiple wave prediction and suppression of seismic data.
    Clear up in 2022.11.05, by Xiang Jian, WPI.

*/

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
#include <algorithm>
#include <omp.h>

#include "../xjc.h"
using namespace std;
using namespace arma;
/////////////////////////////////////////////////////
fmat invcg(fmat x,fmat hess,fmat d,\
    float residual_ratio,int iterations_num);
/////////////////////////////////////////////////////
struct demultiple2d
{
    int ivalNt,ivalNx;
    float fvalDt,fvalDx;
    int FilterOrder1,FilterOrder2,FilterOrder3,FilterOrder4;
    int TimeWinLength,TimeWinSlideLength;

    fmat DataOrigTX,DataToSuppress,DataToSuppressHilbert,\
        DataToSuppressDerivative,DataToSuppressHilbertDerivative;
    fmat dataToSuppressWinBegin,dataToSuppressWinEnd;
    fmat dataProcessingResults,dataBeRemoved;
    fmat &data2d=DataOrigTX,\
        &datadict=DataToSuppress,\
        &datah=DataToSuppressHilbert,\
        &datad=DataToSuppressDerivative,\
        &datahd=DataToSuppressHilbertDerivative;
    fmat &datawinbeg=dataToSuppressWinBegin,\
        &datawinend=dataToSuppressWinEnd;
    fmat &data_remove_wave=dataProcessingResults,\
        &data_antiremove_wave=dataBeRemoved;

    int &n1=ivalNt, &n2=ivalNx;
    int &nt=ivalNt, &nx=ivalNx;
    float&d1=fvalDt, &d2=fvalDx;
    float&dt=fvalDt, &dx=fvalDx;
    int &nwp=FilterOrder1,\
        &nwph=FilterOrder2,\
        &nwpd=FilterOrder3,\
        &nwphd=FilterOrder4;
    int &n1w=TimeWinLength,\
        &d1w=TimeWinSlideLength;
    float lmd;
};
////////// Time Domain Deconvolution 1d for UDD/DGD//////////
void TimeDomainDeconvolution1d_pthread(fmat *data, fmat *deconOperator,\
    fmat *initialDecon, int nt, int nx, int deconDataLength, int deconLength,\
    float lmd, int nx1,int nx2, bool *threadFinished, fmat *decon);
fmat TimeDomainDeconvolution1d(fmat & data, fmat & deconOperator,\
     fmat & initialDecon, int deconDataLength, float lmd=0.1, int ncpu=1)
{
    int icpu,k,nx1,nx2,dnx,nx(data.n_cols),nt(data.n_rows),\
        deconLength(initialDecon.n_rows);
    fmat deconResult=initialDecon;
////////////////////////////////////////////////////////////////////
    fmat *pdata=(&data); 
    fmat *pdeconOperator=(&deconOperator);
    fmat *pdecon=(&deconResult); 
    fmat *pInitial=(&initialDecon);
    ncpu=min(ncpu,nx);

    thread *pcal;
    bool *threadFinished;
    pcal=new thread[ncpu];
    threadFinished=new bool[ncpu];

    for(icpu=0;icpu<ncpu;icpu++){
        threadFinished[icpu]=false;
        pcal[icpu]=thread(TimeDomainDeconvolution1d_pthread,\
            pdata, pdeconOperator,pInitial,nt,nx,deconDataLength,\
            deconLength,lmd,icpu,icpu+1,&threadFinished[icpu],pdecon);
    }
k=ncpu;
while(k<nx){
    for(icpu=0;icpu<ncpu;icpu++){
        if(pcal[icpu].joinable() && threadFinished[icpu]){
            pcal[icpu].join();
            threadFinished[icpu]=false;
        pcal[icpu]=thread(TimeDomainDeconvolution1d_pthread,\
            pdata, pdeconOperator,pInitial,nt,nx,deconDataLength,\
            deconLength,lmd,k,k+1,&threadFinished[icpu],pdecon);
            k++;
        }
    }
}
    for(icpu=0;icpu<ncpu;icpu++){
    if(pcal[icpu].joinable()){
        pcal[icpu].join();
    }}

    delete [] pcal;
    delete [] threadFinished;
    int i,j;
    for(i=0;i<deconLength;i++){
        for(j=0;j<nx;j++){
            if(isnan(abs(deconResult(i,j))))
                {deconResult(i,j)=0;}
        }
    }
    return deconResult;
}

void TimeDomainDeconvolution1d_pthread(fmat *data, fmat *deconOperator,\
    fmat *initialDecon, int nt, int nx, int deconDataLength,int deconLength,\
    float lmd, int nx1,int nx2, bool *threadFinished, fmat *decon)
{
    int i,j,k,wt1,kwt,nw(deconLength),nwt(deconDataLength);
    fmat mat1(nwt,nw),mat0(nwt,nw),data2(nt,nx,fill::zeros),datal,\
        matq(nw,1),mat2(nw,1),matd(nwt,1),matD(nw,nw),matq0(nw,1);
    fmat digmat(nw,nw,fill::zeros);
    digmat.diag()+=lmd;

    for(k=nx1;k<nx2;k++){
        cout<<"now is Deconvolution trace "<<k<<endl;
        kwt=0;
        mat1.fill(0.0);
        mat0.fill(0.0);
        for(i=kwt;i<kwt+nwt;i++){
            for(j=kwt;j<min(kwt+nw,i+1);j++){
                mat0(i-kwt,j-kwt)=deconOperator[0](i-j+kwt,k);
            }
        }
        for(j=kwt;j<kwt+nw;j++){
            mat1.col(j-kwt)=mat0.col(j-kwt);
        }
        for(i=0;i<nwt;i++){
            matd(i,0)=data[0](i+kwt,k);
        }
        matD=mat1.t()*mat1;
        
        fmat dignum(1,1);
        dignum=sum(abs(matd)/nwt,0);
        digmat.fill(0);
        digmat.diag()+=(dignum(0,0)*lmd*nw+0.00001);
        //decon
        {
            matq=inv(matD+digmat)*mat1.t()*matd;
            //matq=invcg(datacol[0](span(0,nw-1),span(0,0)),(matD+digmat),mat1.t()*matd,0.1,nw);
            //matq=invcg(mat1.t()*matd,(matD+digmat),mat1.t()*matd,0.001,nw);
            decon[0](span(0,nw-1),span(k,k))=matq.col(0);
        }
    }
}

//////////////////////////////////////////////////////////////////////
void AdaptiveRemoveMultiple2d(fmat & dataup, fmat & datadown,\
    fmat & data, fmat & dict,\
    fmat &dataph, fmat &datapd, fmat &dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh, int ncpu, int half_winx);
void AdaptiveRemoveMultiple2d(struct demultiple2d & par, int half_winx, int ncpu=1)
{
    AdaptiveRemoveMultiple2d(par.data_remove_wave, par.data_antiremove_wave,\
    par.data2d, par.datadict,\
    par.datah, par.datad, par.datahd,\
    par.n1,par.n2,par.nwp,par.nwph,par.nwpd,par.nwphd,\
    par.n1w,par.d1w,par.lmd,ncpu,half_winx);
}

void AdaptiveRemoveMultiple2dPthread(fmat* dataup, fmat* datadown,\
 fmat* datavz, fmat* datap,fmat*dataph, fmat* datapd, fmat*dataphd,\
 int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,int nwt,int dnwt,\
 float zzh,int nx1,int nx2, int half_winx,bool *end_of_thread);
void AdaptiveRemoveMultiple2d(fmat & dataup, fmat & datadown,\
 fmat & data, fmat & dict,fmat &dataph, fmat &datapd, fmat &dataphd,\
 int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,int nwt,int dnwt,\
 float zzh, int ncpu, int half_winx)
 {
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd),nx1;
    fmat datap(nt,nx),datavz(nt,nx);
////////////////////////////////////////////////////////////////////
    float xs;
    xs=data.max()-data.min();
    datavz=data/(xs);
    //datavz=data;
    datap=dict/(dict.max()-dict.min());
    //datap=dict;
    dataph=dataph/(dataph.max()-dataph.min());
    datapd=datapd/(datapd.max()-datapd.min());
    dataphd=dataphd/(dataphd.max()-dataphd.min());
    dataup.fill(0.0);
////////////////////////////////////////////////////////////////////
    fmat* pdataup=(&dataup); fmat* pdatadown=(&datadown);
    fmat* pdata=(&datavz); fmat* pdict=(&datap);
    fmat* pdataph=(&dataph); fmat* pdatapd=(&datapd);
    fmat* pdataphd=(&dataphd);
    thread* pcal;bool* end_of_thread;
    pcal=new thread[ncpu];
    end_of_thread=new bool[ncpu];

    for(k=0;k<ncpu;k++){
        end_of_thread[k]=false;
        pcal[k]=thread(AdaptiveRemoveMultiple2dPthread,\
            pdataup, pdatadown,pdata, pdict,pdataph, pdatapd,\
            pdataphd,nt, nx,nwp, nwph, nwpd, nwphd,nwt, dnwt,\
            zzh,k,k+1,half_winx,&(end_of_thread[k]));
    }
    nx1=ncpu;
    while(nx1<nx){
    for(k=0;k<ncpu;k++){
    if(end_of_thread[k]&&pcal[k].joinable()&&nx1<nx){
        pcal[k].join();
        end_of_thread[k]=false;
        pcal[k]=thread(AdaptiveRemoveMultiple2dPthread,\
            pdataup, pdatadown,pdata, pdict,pdataph, pdatapd,\
            pdataphd,nt, nx,nwp, nwph, nwpd, nwphd,nwt, dnwt,\
            zzh,nx1,nx1+1,half_winx,&(end_of_thread[k]));
        nx1++;
    }}}
    for(k=0;k<ncpu;k++){
        if(pcal[k].joinable())
            pcal[k].join();
    }
    delete [] pcal;
    delete [] end_of_thread;

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

void AdaptiveRemoveMultiple2dPthread(fmat* dataup, fmat* datadown,\
 fmat* datavz, fmat* datap,fmat*dataph, fmat* datapd, fmat*dataphd,\
 int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,int nwt,int dnwt,\
 float zzh,int nx1,int nx2, int half_winx,bool *end_of_thread)
{
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd);
    fmat digmat(nw,nw,fill::zeros),datal,\
        matq(nw,1),mat2(nw,1),matD(nw,nw),matq0(nw,1);
    digmat.diag()+=1;
    int winbeg,winend,winnum,datanum;
    float *psubmat;

    datal.zeros(nt,1);
    for(k=nx1;k<nx2;k++){
        //if(k%10==0)
            //cout<<"tracl = "<<k<<endl;
        winbeg=max(0,k-half_winx);
        winend=min(int(datap[0].n_cols-1),k+half_winx);
        winnum=winend-winbeg+1;
        datanum=nwt*winnum;
        for(kwt=0;kwt<(nt-nwt);kwt+=dnwt){
            fmat submat(nwt,winnum);
            fmat mat1(datanum,nw),mat0(datanum,nw);
            mat1.fill(0.0);

            mat0.fill(0.0);
            submat=datap[0](span(kwt,kwt+nwt-1),span(winbeg,winend));
            psubmat=&(submat(0,0));
            fvec pdatap(psubmat, datanum); 
            for(i=0;i<datanum;i++){
                for(j=0;j<min(nwp,i+1);j++){
                    //if((i-j)>=0)
                        mat0(i,j)=pdatap(i-j);
                }
            }
            for(j=kwt;j<kwt+nwp;j++){
                mat1.col(j-kwt)=mat0.col(j-kwt);
        }
        if(nwph>0){
            mat0.fill(0.0);
            submat=dataph[0](span(kwt,kwt+nwt-1),span(winbeg,winend));
            psubmat=&(submat(0,0));
            fvec pdataph(psubmat, datanum); 
            for(i=0;i<datanum;i++){
                for(j=0;j<nwph;j++){
                    if((i-j)>=0)
                        mat0(i,j)=pdataph(i-j);
                }
            }
            for(j=kwt;j<kwt+nwph;j++){
                mat1.col(j-kwt+nwp)=mat0.col(j-kwt);
            }
        }
        if(nwpd>0){
            mat0.fill(0.0);
            submat=datapd[0](span(kwt,kwt+nwt-1),span(winbeg,winend));
            psubmat=&(submat(0,0));
            fvec pdatapd(psubmat, datanum); 
            for(i=0;i<datanum;i++){
                for(j=0;j<nwpd;j++){
                    if((i-j)>=0)
                        mat0(i,j)=pdatapd(i-j);
                }
            }
            for(j=kwt;j<kwt+nwpd;j++){
                mat1.col(j-kwt+nwp+nwph)=mat0.col(j-kwt);
            }
        }
        if(nwphd>0){
            mat0.fill(0.0);
            submat=dataphd[0](span(kwt,kwt+nwt-1),span(winbeg,winend));
            psubmat=&(submat(0,0));
            fvec pdataphd(psubmat, datanum); 
            for(i=0;i<datanum;i++){
                for(j=0;j<nwphd;j++){
                    if((i-j)>=0)
                        mat0(i,j)=pdataphd(i-j);
                }
            }
            for(j=kwt;j<kwt+nwphd;j++){
                mat1.col(j-kwt+nwp+nwph+nwpd)=mat0.col(j-kwt);
            }
        }
            submat=datavz[0](span(kwt,kwt+nwt-1),span(winbeg,winend));
            psubmat=&(submat(0,0));
            fmat matd(psubmat, datanum,1); 
            matD=mat1.t()*mat1;   
            fmat dignum(1,1);
            dignum=sum(abs(matd)/nwt/half_winx,0);
            digmat.fill(0);
            digmat.diag()+=(dignum(0,0)*zzh*nw+0.00001);
            //multiple data_win for slip fliter
        {
            matq=inv(matD+digmat)*mat1.t()*matd;
            //matq=invcg(mat1.t()*matd,(matD+digmat),mat1.t()*matd,0.001,nw);
        }
            matd=mat1*matq;
            for(i=nwp;i<nwt;i++){
                datal(i+kwt,0)+=1;
                dataup[0](i+kwt,k)+=(datavz[0](i+kwt,k)\
                    -matd(i+(k-winbeg)*nwt,0));
                datadown[0](i+kwt,k)+=datavz[0](i+kwt,k);
                //datadown[0](i+kwt,k)=datavz[0](i+kwt,k)-dataup[0](i+kwt,k);
            }
        }
        for(i=0;i<nt;i++){
            if(datal(i,0)>0.5){
            dataup[0](i,k)/=datal(i,0);
            datadown[0](i,k)/=datal(i,0);
        }}
        datal.fill(0);
    }
    end_of_thread[0]=true;
}

//////////////////////////////////////////////////////////////////////
void single_trace_dewave_withdatawin_pthread(fmat* dataup, fmat* datadown,\
    fmat* datavz, fmat* datap,fmat*dataph, fmat* datapd, fmat*dataphd,\
    fmat * datawinbeg, fmat * datawinend,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    float zzh,int nx1,int nx2);
void single_trace_dewave_withdatawin(fmat & dataup, fmat & datadown,\
    fmat & data, fmat & dict,fmat &dataph, fmat &datapd, fmat &dataphd,\
    fmat & datawinbeg, fmat & datawinend,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    float zzh, int ncpu)
{
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd),nx1,nx2,dnx;
    fmat datap(nt,nx),datavz(nt,nx);
    float xs;
    xs=data.max()-data.min();
    cout<<xs<<endl;
    datavz=data/(xs);
    datap=dict/(dict.max()-dict.min());
    dataph=dataph/(dataph.max()-dataph.min());
    datapd=datapd/(datapd.max()-datapd.min());
    dataphd=dataphd/(dataphd.max()-dataphd.min());
    fmat* pdataup=(&dataup); fmat* pdatadown=(&datadown);
    fmat* pdatawinbeg=(&datawinbeg); fmat* pdatawinend=(&datawinend);
    fmat* pdata=(&datavz); fmat* pdict=(&datap);
    fmat* pdataph=(&dataph); fmat* pdatapd=(&datapd);
    fmat* pdataphd=(&dataphd);
    thread* pcal;
    pcal=new thread[ncpu];
    dnx=nx/ncpu;
    nx1=0;
    for(k=0;k<ncpu-1;k++){
        nx2=nx1+dnx;
        pcal[k]=thread(single_trace_dewave_withdatawin_pthread,pdataup, pdatadown,\
            pdata, pdict,pdataph, pdatapd, pdataphd,pdatawinbeg, pdatawinend,\
            nt, nx,nwp, nwph, nwpd, nwphd, zzh,nx1,nx2);
        nx1=nx2;
    }
    pcal[ncpu-1]=thread(single_trace_dewave_withdatawin_pthread,pdataup, pdatadown,\
        pdata, pdict,pdataph, pdatapd, pdataphd,pdatawinbeg, pdatawinend,\
        nt, nx,nwp, nwph, nwpd, nwphd, zzh,nx1,nx);
    for(k=0;k<ncpu;k++){
        pcal[k].join();
    }
    delete [] pcal;

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

void single_trace_dewave_withdatawin_pthread(fmat* dataup, fmat* datadown,\
    fmat* datavz, fmat* datap,fmat*dataph, fmat* datapd, fmat*dataphd,\
    fmat * datawinbeg, fmat * datawinend,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    float zzh,int nx1,int nx2)
{
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd),kwtbeg,kwtend,nwt,dnwt;
    fmat matq0(nw,1);
    
    for(k=nx1;k<nx2;k++){
        kwtbeg=datawinbeg[0](k,0);
        kwtbeg=max(0,kwtbeg);
        kwtbeg=min(nt-1,kwtbeg);
        kwtend=datawinend[0](k,0);
        kwtend=max(0,kwtend);
        kwtend=min(nt-1,kwtend);
        nwt=kwtend-kwtbeg;
        {
        fmat mat1(nwt,nw),mat0(nwt,nw),data2(nt,nx,fill::zeros),\
            matq(nw,1),mat2(nw,1),matd(nwt,1),matD(nw,nw);
        fmat digmat(nw,nw,fill::zeros);
        digmat.diag()+=1;
        kwt=kwtbeg;
        //cout<<nwt<<","<<kwtbeg<<"tracl = "<<k<<endl;
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

            matq=inv(matD+digmat)*mat1.t()*matd;
            matq0=matq;
        }

        nwt=nt-5;
        dnwt=nwt-nw-5;
        kwt=0;

        //for(kwt=0;kwt<=(nt-nwt);kwt+=dnwt)
        {
            fmat mat1(nwt,nw),mat0(nwt,nw),data2(nt,nx,fill::zeros),\
            matq(nw,1),mat2(nw,1),matd(nwt,1),matD(nw,nw);
            
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
            
            {
                matq=matq0;
            }
            matd=mat1*matq;
            for(i=nw;i<nwt;i++){
                dataup[0](i+kwt,k)=datavz[0](i+kwt,k)-matd(i,0);
                datadown[0](i+kwt,k)=datavz[0](i+kwt,k)+matd(i,0);
                //datadown[0](i+kwt,k)=datavz[0](i+kwt,k)-dataup[0](i+kwt,k);
            }
        }
    }
}
void single_trace_dewave_withdatawin(struct demultiple2d & par, int ncpu=1)
{
    single_trace_dewave_withdatawin(par.data_remove_wave, par.data_antiremove_wave,\
    par.data2d, par.datadict,par.datah, par.datad, par.datahd,\
    par.datawinbeg, par.datawinend,\
    par.n1,par.n2,par.nwp,par.nwph,par.nwpd,par.nwphd,\
    par.lmd,ncpu);
}

/******* Coding predicts multiple waves, 2D *******/
void multiple_code2d_pthread(cx_fmat* u2cx, cx_fmat* u1cx, fmat* green,\
 float df, int fn1, int fn2);
fmat multiple_code2d(fmat u1, fmat& green,\
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
    delete [] pcal;

    fmat u2;
    u2.copy_size(u1);
    for(k=0;k<n2;k++){
        //u1cx.row(k)=fft(u1.col(k).t(),n1);
        onecol.col(0)=u2cx.row(k).st();
        u2.col(k)=real(ifft(onecol.col(0),n1));
    }
    return u2;
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

/******* Coding predicts multiple waves, 3D *******/
void multiple_code3d_onepoint_onefrequence(cx_fcube* u2, cx_fcube* u1, fmat* green,\
 int isx, int jsy, int kf, float df,bool *finish_work)
{
    int n1(u1[0].n_rows),n2(u1[0].n_cols);
    int i,j;
    float w,pi(3.1415926),t;
    cx_float a;
    
    w=-2.0*pi*df*kf;
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            t=green[0](i,j);
            a.real(0.0);
            a.imag(w*t);
            a=exp(a);
            u2[0](i,j,kf)+=a*u1[0](isx,jsy,kf);
            }
        }
    finish_work[0]=true;
}

void multiple_code3d_onepoint_allfrequence(cx_fcube* u2, cx_fcube* u1,\
 fmat* seabase_depth,fmat *coordx_data,fmat *coordy_data,float docode,\
 int isx, int jsy, int system_source_ix, int system_source_jy, float df,\
 int fn1, int fn2, int minspacewin,int maxspacewin,float water_velocity,\
 float code_pattern, int ncpu,float wavelet_delay,bool *end_of_thread)
{
    if(docode<0.5){
        end_of_thread[0]=true;
    }else{
    int n1(u1[0].n_rows),n2(u1[0].n_cols),n3(u1[0].n_slices);
    int i,j,k,i1,j1,i2,j2,sx,sy,iloopbeg,iloopend,jloopbeg,jloopend;
    float depth, half_offset, d11,d12,d21,d22, l11,l12,l21,l22,fi,fj;
    float l_min(0.0001),l_sum(0.0),l_trace;
    fmat green(n1,n2);
    float w,pi(3.1415926),t;
    code_pattern=1.0-code_pattern;
    code_pattern=std::min(code_pattern,float(1.0));
    code_pattern=std::max(code_pattern,float(0.0));

    if(isx>=system_source_ix){
        iloopbeg=system_source_ix-minspacewin;
        iloopend=isx+minspacewin;
        int wideadd(floor(abs(isx-system_source_ix)*code_pattern));
        wideadd=std::max(wideadd,int(isx-system_source_ix-maxspacewin));
        iloopbeg=iloopbeg+wideadd;
    }else{
        iloopend=system_source_ix+minspacewin;
        iloopbeg=isx-minspacewin;
        int wideadd(floor(abs(system_source_ix-isx)*code_pattern));
        wideadd=std::max(wideadd,int(system_source_ix-isx-maxspacewin));
        iloopend=iloopend-wideadd;
    }
    if(jsy>=system_source_jy){
        jloopbeg=system_source_jy-minspacewin;
        jloopend=jsy+minspacewin;
        int wideadd(floor(abs(jsy-system_source_jy)*code_pattern));
        wideadd=std::max(wideadd,int(jsy-system_source_jy-maxspacewin));
        jloopbeg=jloopbeg+wideadd;
    }else{
        jloopend=system_source_jy+minspacewin;
        jloopbeg=jsy-minspacewin;
        int wideadd(floor(abs(system_source_jy-jsy)*code_pattern));
        wideadd=std::max(wideadd,int(system_source_jy-jsy-maxspacewin));
        jloopend=jloopend-wideadd;
    }
    iloopbeg=max(iloopbeg,0);iloopend=min(iloopend,n1);
    jloopbeg=max(jloopbeg,0);jloopend=min(jloopend,n2);
    //iloopbeg=0;iloopend=n1;
    //jloopbeg=0;jloopend=n2;
    sy=jsy;sx=isx;
    float source_x=coordx_data[0](isx,jsy);
    float source_y=coordy_data[0](isx,jsy);
    for(i=iloopbeg;i<iloopend;i++){      
        i1=floor((sx+i)/2);
        i2=1+i1;
        i1=max(iloopbeg,0);
        i2=min(i2,iloopend-1);
        for(j=jloopbeg;j<jloopend;j++){
            fi=(float(source_x)+coordx_data[0](i,j))/2.0;
            fj=(float(source_y)+coordy_data[0](i,j))/2.0;
            j1=floor((sy+j)/2);
            j2=1+j1;
            j1=max(j1,jloopbeg);
            j2=min(j2,jloopend-1);
            l11=(fi-coordx_data[0](i1,j1))*(fi-coordx_data[0](i1,j1))\
                +(fj-coordy_data[0](i1,j1))*(fj-coordy_data[0](i1,j1))+l_min;
            l12=(fi-coordx_data[0](i1,j2))*(fi-coordx_data[0](i1,j2))\
                +(fj-coordy_data[0](i1,j2))*(fj-coordy_data[0](i1,j2))+l_min;
            l21=(fi-coordx_data[0](i2,j1))*(fi-coordx_data[0](i2,j1))\
                +(fj-coordy_data[0](i2,j1))*(fj-coordy_data[0](i2,j1))+l_min;
            l22=(fi-coordx_data[0](i2,j2))*(fi-coordx_data[0](i2,j2))\
                +(fj-coordy_data[0](i2,j2))*(fj-coordy_data[0](i2,j2))+l_min;
            l11=1/l11;l12=1/l12;l21=1/l21;l22=1/l22;
            l_sum=l11+l12+l21+l22;
            l11/=l_sum;l12/=l_sum;l21/=l_sum;l22/=l_sum;

            d11=seabase_depth[0](i1,j1);
            d12=seabase_depth[0](i1,j2);
            d21=seabase_depth[0](i2,j1);
            d22=seabase_depth[0](i2,j2);
            depth=l11*d11+l21*d21+l12*d12+l22*d22;
            half_offset=sqrt((fi-source_x)*(fi-source_x)+(fj-source_y)*(fj-source_y));
            l_trace=2*sqrt(half_offset*half_offset+depth*depth);
            green(i,j)=l_trace/water_velocity+wavelet_delay;
        }
    }

    cx_float a;
    float blackmanfilter_ibeg,blackmanfilter_iend;
    float blackmanfilter_ibegwide(std::min(minspacewin,isx-iloopbeg)),\
        blackmanfilter_iendwide(std::min(minspacewin,iloopend-isx));
    float blackmanfilter_jbeg,blackmanfilter_jend;
    float blackmanfilter_jbegwide(std::min(minspacewin,jsy-jloopbeg)),\
        blackmanfilter_jendwide(std::min(minspacewin,jloopend-jsy));
    for(k=fn1;k<fn2;k++){
        w=-2.0*pi*df*k;
        float li(0.0),widei(0.0);
        for(i=iloopbeg;i<iloopend;i++){
            if(i-iloopbeg<blackmanfilter_ibegwide){
                blackmanfilter_ibeg=Blackman(i-iloopbeg,blackmanfilter_ibegwide);
            }
            else{
                blackmanfilter_ibeg=1;
            }
            if(iloopend-i<blackmanfilter_iendwide){
                blackmanfilter_iend=Blackman(iloopend-i,blackmanfilter_iendwide);
            }
            else{
                blackmanfilter_iend=1;
            }
            blackmanfilter_ibeg*=blackmanfilter_iend;
            for(j=jloopbeg;j<jloopend;j++){
                float lj(0.0),widej(0.0);
                float l_blackman;
                if(j-jloopbeg<blackmanfilter_jbegwide){
                    blackmanfilter_jbeg=Blackman(j-jloopbeg,blackmanfilter_jbegwide);
                }
                else{
                    blackmanfilter_jbeg=1;
                }
                if(jloopend-j<blackmanfilter_jendwide){
                    blackmanfilter_jend=Blackman(jloopend-j,blackmanfilter_jendwide);
                }
                else{
                    blackmanfilter_jend=1;
                }
                blackmanfilter_jbeg*=blackmanfilter_jend;
                blackmanfilter_jbeg*=blackmanfilter_ibeg;

                t=green(i,j);
                a.real(0.0);
                a.imag(w*t);
                a=exp(a);
                a.real(real(a)*blackmanfilter_jbeg);
                a.imag(imag(a)*blackmanfilter_jbeg);
                u2[0](isx,jsy,k)+=a*u1[0](i,j,k);
        }}}
    end_of_thread[0]=true;
    }
}

void multiple_code3d(cx_fcube& u2, cx_fcube& u1, fmat& seabase_depth,\
 fmat& coordx_data, fmat& coordy_data, fmat& docode,\
 int system_source_ix, int system_source_jy,float water_velocity, \
 float df, int fn1, int fn2,int minspacewin,int maxspacewin,\
 float code_pattern, int ncpu, float wavelet_delay=0.0)
{
    int n1(u1.n_rows),n2(u1.n_cols),n3(u1.n_slices);
    int ix,jy;
    
    cx_fcube *pu2(&u2);
    cx_fcube *pu1(&u1);
    fmat *pcoordx_data(&coordx_data);
    fmat *pcoordy_data(&coordy_data);
    fmat *pseabase(&seabase_depth);
    thread *pcal;
    int kcpu,js;
    bool *end_of_thread;
    end_of_thread=new bool[ncpu];
    pcal=new thread[ncpu];

    for(kcpu=0;kcpu<ncpu;kcpu++){
        jy=0;ix=kcpu;
        end_of_thread[kcpu]=false;
        pcal[kcpu]=thread(multiple_code3d_onepoint_allfrequence,\
            pu2,pu1,pseabase,pcoordx_data,pcoordy_data,docode(ix,jy),ix,jy,\
            system_source_ix,system_source_jy,df,fn1,fn2,minspacewin,\
            maxspacewin,water_velocity, code_pattern,ncpu,wavelet_delay,\
            &(end_of_thread[kcpu]));
    }
    js=ncpu;
    while(js<n1*n2){
        for(kcpu=0;kcpu<ncpu;kcpu++){
        if(!end_of_thread[kcpu])continue;
        else if(end_of_thread[kcpu]){
            if(pcal[kcpu].joinable()&&js<n1*n2){
                pcal[kcpu].join();
                end_of_thread[kcpu]=false;
                jy=int(js/n1);ix=js-jy*n1;
                if(js%n1==0)
                    cout<<jy<<endl;
                pcal[kcpu]=thread(multiple_code3d_onepoint_allfrequence,\
                    pu2,pu1,pseabase,pcoordx_data,pcoordy_data,docode(ix,jy),ix,jy,\
                    system_source_ix,system_source_jy,df,fn1,fn2,minspacewin,\
                    maxspacewin,water_velocity, code_pattern,ncpu,wavelet_delay,\
                    &(end_of_thread[kcpu]));
                js++;
            }
        }}
    }
    for(kcpu=0;kcpu<ncpu;kcpu++){
        if(pcal[kcpu].joinable()){
        pcal[kcpu].join();}
    }
    delete[] pcal;
    delete[] end_of_thread;
/*
    for(jy=0;jy<n2;jy++){
        cout<<"now is run line: "<<jy<<endl;
        for(ix=0;ix<n1;ix++){
            multiple_code3d_onepoint_allfrequence(\
                pu2,pu1,pseabase,pcoordx_data,pcoordy_data,\
                ix,jy,coordx_data(ix,jy), coordy_data(ix,jy),df,fn1,fn2,\
                water_velocity,ncpu,wavelet_delay);
        }
    }
    */
}

/******* Common Functions *******/
void cxfcubeLinearInterpolation3dByCol(cx_fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),i;
    cx_fcube data3dInter;
    int nInter(n2*2-1);
    data3dInter.zeros(n1,nInter,n3);

    data3dInter.col(0)=data3d.col(0);
    for(i=1;i<n2;i++){
        int kinter=2*i;
        data3dInter.col(kinter)=data3d.col(i);
        data3dInter.col(kinter-1)=(data3dInter.col(kinter)\
            +data3dInter.col(kinter-2));
        data3dInter.col(kinter-1)=data3dInter.col(kinter-1)/2.0;
    }
    data3d=data3dInter;
}
void cxfcubeLinearInterpolation3dByRow(cx_fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),i;
    cx_fcube data3dInter;
    int nInter(n1*2-1);
    data3dInter.zeros(nInter,n2,n3);

    data3dInter.row(0)=data3d.row(0);
    for(i=1;i<n1;i++){
        int kinter=2*i;
        data3dInter.row(kinter)=data3d.row(i);
        data3dInter.row(kinter-1)=(data3dInter.row(kinter)\
            +data3dInter.row(kinter-2));
        data3dInter.row(kinter-1)=data3dInter.row(kinter-1)/2.0;
    }
    data3d=data3dInter;
}
void cxfcubeAntiLinearInterpolation3dByRow(cx_fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),kinter;
    cx_fcube data3dInter;
    int nAntiInter((n1+1)/2);
    data3dInter.zeros(nAntiInter,n2,n3);

    for(kinter=0;kinter<nAntiInter;kinter++){
        int i=2*kinter;
        data3dInter.row(kinter)=data3d.row(i);
    }
    data3d=data3dInter;
}
void cxfcubeAntiLinearInterpolation3dByCol(cx_fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),kinter;
    cx_fcube data3dInter;
    int nAntiInter((n2+1)/2);
    data3dInter.zeros(n1,nAntiInter,n3);

    for(kinter=0;kinter<nAntiInter;kinter++){
        int i=2*kinter;
        data3dInter.col(kinter)=data3d.col(i);
    }
    data3d=data3dInter;
}
void fcubeLinearInterpolation3dByRow(fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),i;
    fcube data3dInter;
    int nInter(n1*2-1);
    data3dInter.zeros(nInter,n2,n3);

    data3dInter.row(0)=data3d.row(0);
    for(i=1;i<n1;i++){
        int kinter=2*i;
        data3dInter.row(kinter)=data3d.row(i);
        data3dInter.row(kinter-1)=(data3dInter.row(kinter)\
            +data3dInter.row(kinter-2));
        data3dInter.row(kinter-1)=data3dInter.row(kinter-1)/2.0;
    }
    data3d=data3dInter;
}
void fcubeLinearInterpolation3dByCol(fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),i;
    fcube data3dInter;
    int nInter(n2*2-1);
    data3dInter.zeros(n1,nInter,n3);

    data3dInter.col(0)=data3d.col(0);
    for(i=1;i<n2;i++){
        int kinter=2*i;
        data3dInter.col(kinter)=data3d.col(i);
        data3dInter.col(kinter-1)=(data3dInter.col(kinter)\
            +data3dInter.col(kinter-2));
        data3dInter.col(kinter-1)=data3dInter.col(kinter-1)/2.0;
    }
    data3d=data3dInter;
}
void fcubeAntiLinearInterpolation3dByCol(fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),kinter;
    fcube data3dInter;
    int nAntiInter((n2+1)/2);
    data3dInter.zeros(n1,nAntiInter,n3);

    for(kinter=0;kinter<nAntiInter;kinter++){
        int i=2*kinter;
        data3dInter.col(kinter)=data3d.col(i);
    }
    data3d=data3dInter;
}
void fcubeAntiLinearInterpolation3dByRow(fcube& data3d)
{
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices),kinter;
    fcube data3dInter;
    int nAntiInter((n1+1)/2);
    data3dInter.zeros(nAntiInter,n2,n3);
    for(kinter=0;kinter<nAntiInter;kinter++){
        int i=2*kinter;
        data3dInter.row(kinter)=data3d.row(i);
    }
    data3d=data3dInter;
}
void fmatLinearInterpolation2dByRow(fmat& data2d)
{
    int n1(data2d.n_rows),n2(data2d.n_cols),i;
    fmat data2dInter;
    int nInter(n1*2-1);
    data2dInter.zeros(nInter,n2);
    data2dInter.row(0)=data2d.row(0);
    for(i=1;i<n1;i++){
        int kinter=2*i;
        data2dInter.row(kinter)=data2d.row(i);
        data2dInter.row(kinter-1)=(data2dInter.row(kinter)\
            +data2dInter.row(kinter-2));
        data2dInter.row(kinter-1)=data2dInter.row(kinter-1)/2.0;
    }
    data2d=data2dInter;
}
void fmatLinearInterpolation2dByCol(fmat& data2d)
{
    int n1(data2d.n_rows),n2(data2d.n_cols),i;
    fmat data2dInter;
    int nInter(n2*2-1);
    data2dInter.zeros(n1,nInter);
    data2dInter.col(0)=data2d.col(0);
    for(i=1;i<n2;i++){
        int kinter=2*i;
        data2dInter.col(kinter)=data2d.col(i);
        data2dInter.col(kinter-1)=(data2dInter.col(kinter)\
            +data2dInter.col(kinter-2));
        data2dInter.col(kinter-1)=data2dInter.col(kinter-1)/2.0;
    }
    data2d=data2dInter;
}

fmat get_datawin_tp3d(int nt,int npx,int npy,float dt,\
 float dpx,float dpy,float px0,float py0,float zr_obndepth,\
 float v_water=1500.0,int nx=20001,int ny=20001,float dx=0.1,float dy=0.1)
{
    float v(v_water),zr(zr_obndepth);
    int i,j,k,kx,ky,kpx,kpy;
    fmat datatime0(npx,npy);
    fmat datatime0x;
    fmat datatime0y;

    float xr,yr,t,t0,px,py;
    datatime0.fill(0);
    for(kx=0;kx<nx;kx++){
        xr=(kx*dx);
        for(ky=0;ky<ny;ky++){
            yr=(ky*dy);
            t=(1/v)*sqrt(zr*zr+xr*xr+yr*yr);
            px=(1/v)*(xr/sqrt(zr*zr+xr*xr+yr*yr));
            py=(1/v)*(yr/sqrt(zr*zr+xr*xr+yr*yr));
            t0=t-px*xr-py*yr;
            kpx=round((px-px0)/dpx);
            kpy=round((py-py0)/dpy);
            if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0){
                datatime0(kpx,kpy)=t0;}
            kpx=round((-px-px0)/dpx);
            kpy=round((py-py0)/dpy);
            if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0){
                datatime0(kpx,kpy)=t0;}
            kpx=round((px-px0)/dpx);
            kpy=round((-py-py0)/dpy);
            if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0){
                datatime0(kpx,kpy)=t0;}
            kpx=round((-px-px0)/dpx);
            kpy=round((-py-py0)/dpy);
            if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0){
                datatime0(kpx,kpy)=t0;}
    }}

    datatime0x=datatime0;
    datatime0y=datatime0;
    for(kpy=0;kpy<npy;kpy++){
        t0=0;
    for(kpx=0;kpx<npx;kpx++){
        if(datatime0y(kpx,kpy)>0)
        {t0=datatime0y(kpx,kpy);}
        else
        {datatime0y(kpx,kpy)=t0;}
    }
    for(kpx=npx-1;kpx>=0;kpx--){
        if(datatime0y(kpx,kpy)>0)
        {t0=datatime0y(kpx,kpy);}
        else
        {datatime0y(kpx,kpy)=t0;}
    }
    }
    for(kpx=0;kpx<npx;kpx++){
        t0=0;
    for(kpy=0;kpy<npy;kpy++){
        if(datatime0x(kpx,kpy)>0)
        {t0=datatime0x(kpx,kpy);}
        else
        {datatime0x(kpx,kpy)=t0;}
    }
    for(kpy=npy-1;kpy>=0;kpy--){
        if(datatime0x(kpx,kpy)>0)
        {t0=datatime0x(kpx,kpy);}
        else
        {datatime0x(kpx,kpy)=t0;}
    }
    }
    for(kpx=0;kpx<npx;kpx++){
    for(kpy=0;kpy<npy;kpy++){
        datatime0(kpx,kpy)=max(datatime0y(kpx,kpy),datatime0x(kpx,kpy));
        datatime0(kpx,kpy)=round(datatime0(kpx,kpy)/dt);
    }}
    return datatime0;
}
void getSourceIndes(int& sxindex,int& syindex,\
    fmat& coordx, fmat& coordy)
{
    int i,j;
    float offset,d(1e30);
    for(i=0;i<coordx.n_rows;i++){
    for(j=0;j<coordx.n_cols;j++){
        offset=coordx(i,j)*coordx(i,j)+coordy(i,j)*coordy(i,j);
        if(offset<d){
            d=offset;
            sxindex=i;
            syindex=j;
    }}}
}

/******* Non-functional function *******/
bool*** newboolmat(int x1, int x2, int x3)
{
    bool ***p;
    int i,j,k;
    p=new bool**[x1];   
    for(j=0;j<x1;j++){
        p[j]=new bool*[x2];
        for(i=0;i<x2;i++){
            p[j][i]=new bool[x3];
            for(k=0;k<x3;k++){
                p[j][i][k]=true;
            }
        }}
    return p;
}
fcube fcubesmooth(fcube mat2, int n_smooth)
{
    int i,j,k,n;
    int x1(mat2.n_rows), x2(mat2.n_cols), x3(mat2.n_slices);
    fcube mat1(x1,x2,x3);
    mat1=mat2;

    for(n=0;n<n_smooth;n++){
        for(i=0;i<x1;i++){
            for(j=0;j<x2;j++){
                for(k=0;k<x3;k++){
                if(i>=1 && i<x1-1 && j>=1 && j<x2-1){
                    mat1(i,j,k)=(0.4/0.6)*(((mat2(i-1,j-1,k)*1.0/12+mat2(i-1,j,k)*1.0/6+mat2(i-1,j+1,k)*1.0/12\
                        +mat2(i,j-1,k)*1.0/6+mat2(i,j+1,k)*1.0/6+mat2(i+1,j-1,k)*1.0/12\
                        +mat2(i+1,j,k)*1.0/6+mat2(i+1,j+1,k)*1.0/12+mat2(i,j,k)*1.0/3))*(3.0/4.0))\
                        +(0.2/0.6)*(((mat2(i-1,j-1,k+1)*1.0/12+mat2(i-1,j,k+1)*1.0/6+mat2(i-1,j+1,k+1)*1.0/12\
                        +mat2(i,j-1,k+1)*1.0/6+mat2(i,j+1,k+1)*1.0/6+mat2(i+1,j-1,k+1)*1.0/12\
                        +mat2(i+1,j,k+1)*1.0/6+mat2(i+1,j+1,k+1)*1.0/12+mat2(i,j,k+1)*1.0/3))*(3.0/4.0))\
                        +(0.2/0.6)*(((mat2(i-1,j-1,k-1)*1.0/12+mat2(i-1,j,k-1)*1.0/6+mat2(i-1,j+1,k-1)*1.0/12\
                        +mat2(i,j-1,k-1)*1.0/6+mat2(i,j+1,k-1)*1.0/6+mat2(i+1,j-1,k-1)*1.0/12\
                        +mat2(i+1,j,k-1)*1.0/6+mat2(i+1,j+1,k-1)*1.0/12+mat2(i,j,k-1)*1.0/3))*(3.0/4.0));
                }
                if(i==0)
                {mat1(i,j,k)=mat1(i+1,j,k);}
                if(j==0)
                {mat1(i,j,k)=mat1(i,j+1,k);}
                if(k==0)
                {mat1(i,j,k)=mat1(i,j,k+1);}

                if(k==x3-1)
                {mat1(i,j,k)=mat1(i,j,k-1);}
                if(i==x1-1)
                {mat1(i,j,k)=mat1(i-1,j,k);}
                if(j==x2-1)
                {mat1(i,j,k)=mat1(i,j-1,k);}
            }}}
        mat2=mat1;
    }
    return mat1;
}
void get_taup3d_CCA(fcube &filterdata, fcube &data1, fcube &data2,\
    float wx=20, float wy=20, float wz=20,\
    float Butterworth_factor=1, float Butterworth_n=2, int n_smooth=9)
{
    int i,j,k,i1,j1,nx,ny,nz,k1,k2;
    ny=data1.n_cols;
    nx=data1.n_rows;
    nz=data1.n_slices;
    float agcpow1,agcpow2,Butterworth;
    fcube datacov1(nx,ny,nz,fill::zeros);
    fcube datacov2(nx,ny,nz,fill::zeros);
    fcube filter(nx,ny,nz,fill::zeros);

    for(i=wx;i<nx-wx;i++){
    for(j=wy;j<ny-wy;j++){
    for(k=wz;k<nz-wz;k++){
        agcpow1=0;
        agcpow2=0;
        for(i1=-wx;i1<=wx;i1++){
        for(j1=-wy;j1<=wy;j1++){
        for(k1=-wz;k1<=wz;k1++){
            agcpow1+=(data1(i+i1,j+j1,k+k1)*data1(i-i1,j+j1,k+k1)*\
                data1(i+i1,j-j1,k+k1)*data1(i-i1,j-j1,k+k1))*\
                (data1(i+i1,j+j1,k-k1)*data1(i-i1,j+j1,k-k1)*\
                data1(i+i1,j-j1,k-k1)*data1(i-i1,j-j1,k-k1));
            agcpow2+=(data2(i+i1,j+j1,k+k1)*data2(i-i1,j+j1,k+k1)*\
                data2(i+i1,j-j1,k+k1)*data2(i-i1,j-j1,k+k1))*\
                (data2(i+i1,j+j1,k-k1)*data2(i-i1,j+j1,k-k1)*\
                data2(i+i1,j-j1,k-k1)*data2(i-i1,j-j1,k-k1));
        }}}
        datacov1(i,j,k)=abs(agcpow1);
        datacov2(i,j,k)=abs(agcpow2);
        agcpow1=0;
        agcpow2=0;
    }}}
    datacov1=fcubesmooth(datacov1,n_smooth);
    datacov1=datacov1/datacov1.max();
    datacov2=fcubesmooth(datacov2,n_smooth);
    datacov2=datacov2/datacov2.max();

    for(i=wx;i<nx-wx;i++){
    for(j=wy;j<ny-wy;j++){
    for(k=wz;k<nz-wz;k++){
        filter(i,j,k)=Butterworth_factor*datacov1(i,j,k)/datacov2(i,j,k);
        Butterworth=1/sqrt(1+pow(filter(i,j,k),Butterworth_n));
        filterdata(i,j,k)=data1(i,j,k)-data1(i,j,k)*(Butterworth);
    }}}
}
void get_taup2d_CCA(fmat &filterdata, fmat &datacov_out, fmat &data1, fmat &data2,\
 int ncpu, float wx=20, float wy=20,\
 float Butterworth_factor=1, float Butterworth_n=2,int Butterworth_smooth=99)
{
    int i,j,i1,j1,nx,ny,k1,k2;
    ny=data1.n_cols;
    nx=data1.n_rows;
    float agcpow1,agcpow2,Butterworth;
    fmat datacov1(nx,ny,fill::zeros);
    fmat datacov2(nx,ny,fill::zeros);
    fmat datacov1_unit(nx,ny,fill::zeros);
    fmat datacov2_unit(nx,ny,fill::zeros);
    fmat filter(nx,ny,fill::zeros);
    k1=(nx-wx);
omp_set_num_threads(ncpu);
#pragma omp parallel for
    for(i=wx;i<k1;i++){
    for(j=wy;j<ny-wy;j++){
        agcpow1=0;
        agcpow2=0;
        for(i1=-wx;i1<=wx;i1++){
        for(j1=-wy;j1<=wy;j1++){
            agcpow1+=(data1(i+i1,j+j1)*data1(i-i1,j+j1)*\
                data1(i+i1,j-j1)*data1(i-i1,j-j1));
            agcpow2+=(data2(i+i1,j+j1)*data2(i-i1,j+j1)*\
                data2(i+i1,j-j1)*data2(i-i1,j-j1));
        }}
        datacov1(i,j)=abs(agcpow1);
        datacov2(i,j)=abs(agcpow2);
        agcpow1=0;
        agcpow2=0;
    }}
/*
    for(i=wx;i<nx-wx;i++){
    for(j=wy;j<ny-wy;j++){
        agcpow1=0;
        agcpow2=0;
        for(i1=-wx;i1<=wx;i1++){
        for(j1=-wy;j1<=wy;j1++){
            agcpow1+=datacov1(i+i1,j+j1);
            agcpow2+=datacov2(i+i1,j+j1);
        }}
        datacov1_unit(i,j)=datacov1(i,j)/(agcpow1);
        datacov2_unit(i,j)=datacov2(i,j)/(agcpow2);
    }}
    datacov2=datacov2_unit;
    datacov1=datacov1_unit;
*/
    datacov1=fmatsmooth(datacov1,nx,ny,Butterworth_smooth);
    datacov1=datacov1/datacov1.max();
    datacov2=fmatsmooth(datacov2,nx,ny,Butterworth_smooth);
    datacov2=datacov2/datacov2.max();
    for(i=wx;i<nx-wx;i++){
    for(j=wy;j<ny-wy;j++){
        filter(i,j)=Butterworth_factor*datacov1(i,j)/datacov2(i,j);
        Butterworth=1/sqrt(1+pow(filter(i,j),Butterworth_n));  
        filterdata(i,j)=data1(i,j)-data1(i,j)*(Butterworth);
        datacov_out(i,j)=Butterworth;
    }}
}
//////////////////////////////////////////////////////////////////////
fmat invcg(fmat x,fmat hess,fmat d,\
    float residual_ratio=0.1,int iterations_num=45)
{
    int ip,jp,in,jn,k,iter(0);
    int np1(x.n_rows),np2(x.n_cols);
    cout<<np1<<"||"<<np2<<endl;
    
    fmat gradient_rk,gradient_rk_1,\
        gradient_cg_pk,gradient_cg_pk_1,\
        datatp_k,datatp_k_1,\
        recoverdatatx_uk,A_gradient_cg_pk,A_datatp_k;
    fmat sum_num(1,1),residual_pow(1,1),residual_k(1,1);
    fmat beta_k(1,1),alpha_k(1,1);
    datatp_k.copy_size(x);
    datatp_k_1.copy_size(datatp_k);
    gradient_rk.copy_size(datatp_k);
    gradient_rk_1.copy_size(datatp_k);
    gradient_cg_pk.copy_size(datatp_k);
    gradient_cg_pk_1.copy_size(datatp_k); 
    A_gradient_cg_pk.copy_size(datatp_k);
    A_datatp_k.copy_size(datatp_k);

    iter=0;
    datatp_k=x;
    sum_num=sum(abs(datatp_k));
    //datatp_k.fill(0.0);
    datatp_k=((datatp_k)/1.0);
    //sum_num=sum(sum(sum(abs(datatp_k))));
    residual_pow=sum_num(0,0);
    residual_pow*=residual_ratio;

    //cxfmatget_Ap_small(A_datatp_k,hess_A,datatp_k,kcpu,par);
    A_datatp_k=hess*datatp_k;
    gradient_rk=d-A_datatp_k;
    gradient_cg_pk=gradient_rk;

//cxfmatget_Ap_small(A_gradient_cg_pk,hess_A,gradient_cg_pk,kcpu,par);
    A_gradient_cg_pk=hess*gradient_cg_pk;
    alpha_k=(gradient_rk.t(),gradient_rk);
    alpha_k=alpha_k/(gradient_cg_pk.t(),A_gradient_cg_pk);
    //get new solution
    datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
    gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;

    beta_k=(gradient_rk_1.t(),gradient_rk_1);
    beta_k=beta_k/(gradient_rk.t(),gradient_rk);
    gradient_cg_pk_1=gradient_rk_1+beta_k(0,0)*gradient_cg_pk;
    //updata
    datatp_k=datatp_k_1;
    gradient_rk=gradient_rk_1;
    gradient_cg_pk=gradient_cg_pk_1;
    //cal residual_pow
    sum_num=sum(abs(gradient_rk));
    residual_k=sum_num;

    while(iter<iterations_num && residual_k(0,0)>residual_pow(0,0)){
        iter++;

//cxfmatget_Ap_small(A_gradient_cg_pk,hess_A,gradient_cg_pk,kcpu,par);
        A_gradient_cg_pk=hess*gradient_cg_pk;

        alpha_k=(gradient_rk.t()*gradient_rk);
        alpha_k=alpha_k/(gradient_cg_pk.t()*A_gradient_cg_pk);
        //get new solution
        datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
        gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;

        beta_k=(gradient_rk_1.t()*gradient_rk_1);
        beta_k=beta_k/(gradient_rk.t()*gradient_rk);
        gradient_cg_pk_1=gradient_rk_1+(beta_k(0,0))*gradient_cg_pk;
        //updata
        datatp_k=datatp_k_1;
        gradient_rk=gradient_rk_1;
        gradient_cg_pk=gradient_cg_pk_1;
        //cal residual_pow
        sum_num=sum(abs(gradient_rk));
        residual_k=sum_num;
    }
    cout<<residual_pow(0,0)<<"||"<<iter<<"||"\
        <<residual_k(0,0)/residual_pow(0,0)<<"||"<<sum_num(0,0)<<endl;
    return datatp_k;
    //sum_num=sum(sum(sum(abs(datatp_k))));
    
/*
    if((residual_k(0,0)/residual_pow(0,0))>(10.0)\
        ||isnan(residual_k(0,0))){
            par[0].datafP.slice(kf).fill(0.0);
            convergence[0]=false;
        }
        else{
            par[0].datafP.slice(kf)=datatp_k;
            convergence[0]=true;
        }
*/
}



#endif
