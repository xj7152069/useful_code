/************update in 2021.01.07**************/
/*
    
***********************************************/

#ifndef SLANTSTACK3D_HPP
#define SLANTSTACK3D_HPP
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

// 2-D Slant-stack in the time - space domain.
int Slantstack_CG_2D(float **trace, int ntrace, int ns, float dt,\
 float *coor, float x0, float **tauppanel, int npsr,float psrmin, float dpsr,\
 float **recoverdata, bool dorecover, int ncpu, \
 float factor_L2, int iterations_num, float residual_ratio);
/* discription.
This function will do the local LRT over local traces .
float **trace;
the INPUT local seismic-gather, i.e. trace[ntrace][ns].
int ntrace;
the trace-number of the INPUT gather.
int ns;
the trace length.
float dt;
time-sampling interval, (unit of dt should be second) ,
float *coor;
coor[ntrace] stores the offset of each trace .
float x0;
the beam-center coordinate .
float **tauppanel; 
the OUTPUT tau-p spectrum, i.e. tauppanel[npsr][ns].
int npsr;
the ray-parameter number of tau-p panel.
float psrmin;
the minimum ray-parameter, (unit is s/m) .
float dpsr;
the interval of ray-parameter， (unit is s/m) .
float **recoverdata;
the OUTPUT recover local seismic-gather, i.e. recoverdata[ntrace][ns].
bool dorecover; default true
Whether to OUTPUT recover local seismic-gather
float factor_L2; default 1.0
L2, Tikhonov regularization parameter
float fmax; default 150
The maximum frequency of the seismic signal
int ncpu; default 1
The number of threads computed by multithreading
int iterations_num; default 48
the Maximum number of iterations of conjugate-gradient method
float residual_ratio; default 0.5
Conjugate gradient method iteration error, The smaller the more accurate
*/
///////////////////////////////////////////////////////////////////
/*slantstack3d变换 传递的参数：
n: number of
d: Time or spatial sampling interval
p: Ray parameters
*/
struct slantstack3d
{
    int ntrace,nline,nt,np_trace,np_line,numthread;
    float d_trace,d_line,dt,dp_trace,dp_line,\
        ptrace_coord0,pline_coord0,line_coord0,trace_coord0;
    fcube weighted_fcube,datatx,datatp,recoverdatatx;
    fmat pline_coord,ptrace_coord,nline_coordy,ntrace_coordx;
    float factor_l1,factor_l2;
};
//Sets a General parameters of initialization parameters
void slantstack3d_parset(int ntrace, int nline, int nt,\
    struct slantstack3d & par)
{
    par.ntrace=ntrace,par.nt=nt,par.nline=nline;
    par.d_trace=10,par.d_line=10,par.dt=0.001;
    par.numthread=1;
    par.factor_l1=(1);
    par.factor_l2=(1);
} 

void slantstack3d_parupdate(struct slantstack3d & par)
{ 
    par.datatx.zeros(par.ntrace,par.nline,par.nt);
    par.recoverdatatx.zeros(par.ntrace,par.nline,par.nt);
    par.ntrace_coordx.zeros(par.ntrace,par.nline);
    par.nline_coordy.zeros(par.ntrace,par.nline);
    par.datatp.zeros(par.np_trace,par.np_line,par.nt);
    par.weighted_fcube.zeros(par.np_trace,par.np_line,par.nt);
    par.ptrace_coord.zeros(par.np_trace,1);
    par.pline_coord.zeros(par.np_line,1);
    
    int k;
    for(k=0;k<par.np_trace;k++)
    {
        par.ptrace_coord(k,0)=k*par.dp_trace+par.ptrace_coord0;
    }
    for(k=0;k<par.np_line;k++)
    {
        par.pline_coord(k,0)=k*par.dp_line+par.pline_coord0;
    }
}

void slantstack3d_cleardata(struct slantstack3d & par)
{
    par.datatp.fill(0.0);
    par.recoverdatatx.fill(0.0);
}

/////////////////////////slantstack3d_L_LT///////////////////////////

void slantstack3d_stack_LT_pthread(fcube *datatp,fcube *datatx,\
 fmat *ptrace_coord,fmat *pline_coord,fmat *ntrace_coord,fmat *nline_coord,\
 float dt, int kt,bool *end_of_thread)
{

    int ip,jp,in,jn,k,ntaoceil,ntaofloor;
    float tao0,dtao1,dtao2,dtao12,ntao12,wceil,wfloor,wsum,wmin(0.000000001);
    int n1(datatx[0].n_rows),n2(datatx[0].n_cols),n3(datatx[0].n_slices),\
        np1(datatp[0].n_rows),np2(datatp[0].n_cols);
        
k=kt;
{
    tao0=k*dt;
    for(ip=0;ip<np1;ip++){
        for(jp=0;jp<np2;jp++){
            for(in=0;in<n1;in++){
                for(jn=0;jn<n2;jn++){
                    dtao1=ntrace_coord[0](in,jn)*ptrace_coord[0](ip,0);
                    dtao2=nline_coord[0](in,jn)*pline_coord[0](jp,0);
                    dtao12=dtao1+dtao2;
                    ntao12=(tao0+dtao12)/dt;
                    if(ntao12<(n3-1) && ntao12>0){
                        ntaoceil=ceil(ntao12);
                        ntaofloor=floor(ntao12);
                        wceil=1.0/(abs(ntaoceil-ntao12)+wmin);
                        wfloor=1.0/(abs(ntaofloor-ntao12)+wmin);
                        wsum=wceil+wfloor;
                        wceil=wceil/wsum;
                        wfloor=wfloor/wsum;

                        datatp[0](ip,jp,k)=datatp[0](ip,jp,k)\
                            +datatx[0](in,jn,ntaoceil)*wceil\
                            +datatx[0](in,jn,ntaofloor)*wfloor;
                    }
    }}}}}
    end_of_thread[0]=true;
}

void slantstack3d_stack_LT_operator(fcube &datatp,fcube &datatx,\
 fmat &ptrace_coord,fmat &pline_coord,fmat &ntrace_coord,fmat &nline_coord,\
 float dt, int ncpu){

    int i,j,k,nt(datatx.n_slices);
    datatp.fill(0.0);
    fcube *pdatatp(&datatp);
    fcube *pdatatx(&datatx);
    fmat *pptrace_coord(&ptrace_coord);
    fmat *ppline_coord(&pline_coord);
    fmat *pntrace_coord(&ntrace_coord);
    fmat *pnline_coord(&nline_coord);
    bool *end_of_thread;
    end_of_thread=new bool[ncpu];

    thread *pcal;
    pcal=new thread[ncpu];
    for(k=0;k<min(ncpu,nt);k++){
        end_of_thread[k]=false;
        pcal[k]=thread(slantstack3d_stack_LT_pthread,pdatatp,pdatatx,\
            pptrace_coord,ppline_coord,pntrace_coord,pnline_coord,\
                dt, k,&(end_of_thread[k]));
    }
    i=ncpu;
    while(i<nt){
        for(k=0;k<ncpu;k++){
            if(pcal[k].joinable() && i<nt && end_of_thread[k]){
                pcal[k].join();
                end_of_thread[k]=false;
                pcal[k]=thread(slantstack3d_stack_LT_pthread,pdatatp,pdatatx,\
                    pptrace_coord,ppline_coord,pntrace_coord,pnline_coord,\
                    dt, i,&(end_of_thread[k]));
                i++;
                //cout<<i<<endl;
        }}
    }
    for(k=0;k<min(ncpu,nt);k++){
        if(pcal[k].joinable()){
        pcal[k].join();}
    }
    delete [] pcal;
    delete [] end_of_thread;
}

void slantstack3d_stack_LT_operator(fcube &datatp,fcube &datatx,\
 fmat &ptrace_coord,fmat &pline_coord,fmat &ntrace_coord,fmat &nline_coord,\
 float dt){

    int ip,jp,in,jn,k,ntaoceil,ntaofloor;
    float tao0,dtao1,dtao2,dtao12,ntao12,wceil,wfloor,wsum,wmin(0.000000001);
    int n1(datatx.n_rows),n2(datatx.n_cols),n3(datatx.n_slices),\
        np1(datatp.n_rows),np2(datatp.n_cols);
    datatp.fill(0.0);

for(k=0;k<n3;k++){
    tao0=k*dt;
    for(ip=0;ip<np1;ip++){
        for(jp=0;jp<np2;jp++){
            for(in=0;in<n1;in++){
                for(jn=0;jn<n2;jn++){
                    dtao1=ntrace_coord(in,jn)*ptrace_coord(ip,0);
                    dtao2=nline_coord(in,jn)*pline_coord(jp,0);
                    dtao12=dtao1+dtao2;
                    ntao12=(tao0+dtao12)/dt;
                    if(ntao12<(n3-1) && ntao12>0){
                        ntaoceil=ceil(ntao12);
                        ntaofloor=floor(ntao12);
                        wceil=1.0/(abs(ntaoceil-ntao12)+wmin);
                        wfloor=1.0/(abs(ntaofloor-ntao12)+wmin);
                        wsum=wceil+wfloor;
                        wceil=wceil/wsum;
                        wfloor=wfloor/wsum;

                        datatp(ip,jp,k)=datatp(ip,jp,k)+datatx(in,jn,ntaoceil)*wceil\
                            +datatx(in,jn,ntaofloor)*wfloor;
                    }
    }}}}}
}
/////////////////////////////////////////////////////////////////////////
void slantstack3d_recover_L_pthread(fcube *recoverdatatx,fcube *datatp,\
 fmat *ptrace_coord,fmat *pline_coord,fmat *ntrace_coord,fmat *nline_coord,\
 float dt,int kt,bool *end_of_thread)
{
    int ip,jp,in,jn,k,ntaoceil,ntaofloor;
    float tao0,dtao1,dtao2,dtao12,ntao12,wceil,wfloor,wsum,wmin(0.000000001);
    int n1(recoverdatatx[0].n_rows),n2(recoverdatatx[0].n_cols),n3(recoverdatatx[0].n_slices),\
        np1(datatp[0].n_rows),np2(datatp[0].n_cols);

k=kt;
{
    tao0=k*dt;
    for(in=0;in<n1;in++){
        for(jn=0;jn<n2;jn++){
            for(ip=0;ip<np1;ip++){
                for(jp=0;jp<np2;jp++){
                    dtao1=ntrace_coord[0](in,jn)*ptrace_coord[0](ip,0);
                    dtao2=nline_coord[0](in,jn)*pline_coord[0](jp,0);
                    dtao12=dtao1+dtao2;
                    ntao12=(tao0-dtao12)/dt;
                    if(ntao12<(n3-1) && ntao12>0){
                        ntaoceil=ceil(ntao12);
                        ntaofloor=floor(ntao12);
                        wceil=1.0/(abs(ntaoceil-ntao12)+wmin);
                        wfloor=1.0/(abs(ntaofloor-ntao12)+wmin);
                        wsum=wceil+wfloor;
                        wceil=wceil/wsum;
                        wfloor=wfloor/wsum;

                        recoverdatatx[0](in,jn,k)=recoverdatatx[0](in,jn,k)\
                            +datatp[0](ip,jp,ntaoceil)*wceil\
                            +datatp[0](ip,jp,ntaofloor)*wfloor;
                    }
    }}}}}
    end_of_thread[0]=true;
}
void slantstack3d_recover_L_operator(fcube &recoverdatatx,fcube &datatp,\
 fmat &ptrace_coord,fmat &pline_coord,fmat &ntrace_coord,fmat &nline_coord,\
 float dt, int ncpu){

    int i,j,k,nt(recoverdatatx.n_slices);
    recoverdatatx.fill(0.0);
    fcube *pdatatp(&datatp);
    fcube *pdatatx(&recoverdatatx);
    fmat *pptrace_coord(&ptrace_coord);
    fmat *ppline_coord(&pline_coord);
    fmat *pntrace_coord(&ntrace_coord);
    fmat *pnline_coord(&nline_coord);
    bool *end_of_thread;
    end_of_thread=new bool[ncpu];

    thread *pcal;
    pcal=new thread[ncpu];
    for(k=0;k<min(ncpu,nt);k++){
        end_of_thread[k]=false;
        pcal[k]=thread(slantstack3d_recover_L_pthread,pdatatx,pdatatp,\
            pptrace_coord,ppline_coord,pntrace_coord,pnline_coord,\
                dt, k,&(end_of_thread[k]));
    }
    i=ncpu;
    while(i<nt){
        for(k=0;k<ncpu;k++){
            if(pcal[k].joinable() && i<nt && end_of_thread[k]){
                pcal[k].join();
                end_of_thread[k]=false;
                pcal[k]=thread(slantstack3d_recover_L_pthread,pdatatx,pdatatp,\
                    pptrace_coord,ppline_coord,pntrace_coord,pnline_coord,\
                    dt, i,&(end_of_thread[k]));
                i++;
                //cout<<i<<endl;
        }}
    }
    for(k=0;k<min(ncpu,nt);k++){
        pcal[k].join();
    }
    delete [] pcal;
    delete [] end_of_thread;
}

void slantstack3d_recover_L_operator(fcube &recoverdatatx,fcube &datatp,\
 fmat &ptrace_coord,fmat &pline_coord,fmat &ntrace_coord,fmat &nline_coord,\
 float dt){

    int ip,jp,in,jn,k,ntaoceil,ntaofloor;
    float tao0,dtao1,dtao2,dtao12,ntao12,wceil,wfloor,wsum,wmin(0.000000001);
    int n1(recoverdatatx.n_rows),n2(recoverdatatx.n_cols),n3(recoverdatatx.n_slices),\
        np1(datatp.n_rows),np2(datatp.n_cols);
    recoverdatatx.fill(0.0);

for(k=0;k<n3;k++){
    tao0=k*dt;
    for(in=0;in<n1;in++){
        for(jn=0;jn<n2;jn++){
            for(ip=0;ip<np1;ip++){
                for(jp=0;jp<np2;jp++){
                    dtao1=ntrace_coord(in,jn)*ptrace_coord(ip,0);
                    dtao2=nline_coord(in,jn)*pline_coord(jp,0);
                    dtao12=dtao1+dtao2;
                    ntao12=(tao0-dtao12)/dt;
                    if(ntao12<(n3-1) && ntao12>0){
                        ntaoceil=ceil(ntao12);
                        ntaofloor=floor(ntao12);
                        wceil=1.0/(abs(ntaoceil-ntao12)+wmin);
                        wfloor=1.0/(abs(ntaofloor-ntao12)+wmin);
                        wsum=wceil+wfloor;
                        wceil=wceil/wsum;
                        wfloor=wfloor/wsum;

                        recoverdatatx(in,jn,k)=recoverdatatx(in,jn,k)\
                            +datatp(ip,jp,ntaoceil)*wceil\
                            +datatp(ip,jp,ntaofloor)*wfloor;
                    }
    }}}}}
}
//////////////////////////////////////////////////////////////
void slantstack3d_stack(struct slantstack3d &par){
    slantstack3d_stack_LT_operator(par.datatp,par.datatx,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
    par.dt);
}
void slantstack3d_stack_multithread(struct slantstack3d &par){
    slantstack3d_stack_LT_operator(par.datatp,par.datatx,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
    par.dt,par.numthread);
}
void slantstack3d_recover(struct slantstack3d &par){
    slantstack3d_recover_L_operator(par.recoverdatatx,par.datatp,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
    par.dt);
}
void slantstack3d_recover_multithread(struct slantstack3d &par){
    slantstack3d_recover_L_operator(par.recoverdatatx,par.datatp,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
    par.dt,par.numthread);
}
/////////////////////////slantstack_CG3d///////////////////////////
inline float fcube_inner_product(fcube data1, fcube data2){
    int ni(data1.n_rows),nj(data1.n_cols),nk(data1.n_slices),i,j,k;
    float inner_num(0);
    for(i=0;i<ni;i++){
        for(j=0;j<nj;j++){
            for(k=0;k<nk;k++){
                inner_num+=(data1(i,j,k)*data2(i,j,k));
    }}}
    return inner_num;
} 

//////////////////////////////////////////////////////////////////////
void fcubemul(fcube& outdata, fcube& mat1, fcube& mat2)
{
    int i,j,k,nx,ny,nz;
    ny=outdata.n_cols;
    nx=outdata.n_rows;
    nz=outdata.n_slices;

    for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
    for(k=0;k<nz;k++){
        outdata(i,j,k)=mat1(i,j,k)*mat2(i,j,k);
    }}}
}
void get_weighted_fcube(fcube &filter, fcube &data1,\
 int wx=5, int wy=5, int wz=10)
{
    int i,j,k,i1,j1,nx,ny,nz,k1,k2;
    int win;
    ny=data1.n_cols;
    nx=data1.n_rows;
    nz=data1.n_slices;
    win=(nx-1)/2;
    wx=min(win,wx);
    win=(ny-1)/2;
    wy=min(win,wy);
    win=(nz-1)/2;
    wz=min(win,wz);

    float agcpow1;
    for(i=wx;i<nx-wx;i++){
    for(j=wy;j<ny-wy;j++){
    for(k=wz;k<nz-wz;k++){
        agcpow1=0;
        for(i1=-wx;i1<=wx;i1++){
        for(j1=-wy;j1<=wy;j1++){
        for(k1=-wz;k1<=wz;k1++){
            //agcpow1+=(data1(i+i1,j+j1,k+k1)*data1(i-i1,j+j1,k+k1)*\
                data1(i+i1,j-j1,k+k1)*data1(i-i1,j-j1,k+k1))*\
                (data1(i+i1,j+j1,k-k1)*data1(i-i1,j+j1,k-k1)*\
                data1(i+i1,j-j1,k-k1)*data1(i-i1,j-j1,k-k1));
            agcpow1+=data1(i+i1,j+j1,k+k1)*data1(i+i1,j+j1,k+k1);
        }}}
        filter(i,j,k)=sqrt(agcpow1);
        agcpow1=0;
    }}}
    filter=(filter)/filter.max();
    for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
    for(k=0;k<nz;k++){
        filter(i,j,k)=1.0/(1.0+exp(10*(0.5-filter(i,j,k))));
    }}}
    filter=(filter)/filter.max();
    filter+=0.00001;
    filter=1.0/filter;
    filter=(filter)/filter.max();
}

void slantstack3d_stack_CG_invL_operator(struct slantstack3d &par,\
 int iterations_num=9, float residual_ratio=1)
{
    get_weighted_fcube(par.weighted_fcube,par.datatp);
    float dataxs(1.0),datamax;
    par.weighted_fcube*=par.factor_l1;
    par.weighted_fcube+=par.factor_l2;
    datawrite3d_bycol_transpose(par.weighted_fcube,par.nt,par.np_trace,"weight.bin");

    int ip,jp,in,jn,k,iter(0);
    int n1(par.datatx.n_rows),n2(par.datatx.n_cols),n3(par.datatx.n_slices),\
        np1(par.datatp.n_rows),np2(par.datatp.n_cols);
    fcube gradient_rk,gradient_rk_1,\
        gradient_cg_pk,gradient_cg_pk_1,\
        datatp_k,datatp_k_1,datatp_weighted,\
        recoverdatatx_uk,A_gradient_cg_pk,A_datatp_k;
    fmat sum_num(1,1);
    float beta_k,alpha_k,residual_pow,residual_k;
    datatp_weighted.copy_size(par.datatp);
    datatp_k.copy_size(par.datatp);
    datatp_k_1.copy_size(par.datatp);
    gradient_rk.copy_size(par.datatp);
    gradient_rk_1.copy_size(par.datatp);
    gradient_cg_pk.copy_size(par.datatp);
    gradient_cg_pk_1.copy_size(par.datatp); 
    recoverdatatx_uk.copy_size(par.datatx);
    A_gradient_cg_pk.copy_size(par.datatp);
    A_datatp_k.copy_size(par.datatp);

    iter=0;datamax=0;
    for(k=0;k<par.ntrace;k++){
        datamax+=par.datatx.row(k).max();
    }
    datamax/=par.ntrace;
    datatp_k=datamax*(par.datatp/par.datatp.max());
    slantstack3d_recover_L_operator(recoverdatatx_uk,datatp_k,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt,par.numthread);
    slantstack3d_stack_LT_operator(A_datatp_k,recoverdatatx_uk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt,par.numthread);
    //regularization
    fcubemul(datatp_weighted,datatp_k,par.weighted_fcube);
    A_datatp_k=dataxs*A_datatp_k+datatp_weighted;

    gradient_rk=par.datatp-A_datatp_k;
    gradient_cg_pk=gradient_rk;

    //cal residual_pow
    sum_num=sum(sum(sum(abs(gradient_rk))));
    residual_k=sum_num(0,0);
    residual_pow=residual_k;
    residual_pow*=residual_ratio;

    slantstack3d_recover_L_operator(recoverdatatx_uk,gradient_cg_pk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt,par.numthread);
    slantstack3d_stack_LT_operator(A_gradient_cg_pk,recoverdatatx_uk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt,par.numthread);
    //regularization
    fcubemul(datatp_weighted,gradient_cg_pk,par.weighted_fcube);
    A_gradient_cg_pk=dataxs*A_gradient_cg_pk+datatp_weighted;
        
    alpha_k=fcube_inner_product(gradient_rk,gradient_rk);
    alpha_k=alpha_k/fcube_inner_product(gradient_cg_pk,A_gradient_cg_pk);
    //get new solution
    datatp_k_1=datatp_k+alpha_k*gradient_cg_pk;
    gradient_rk_1=gradient_rk-alpha_k*A_gradient_cg_pk;

    beta_k=fcube_inner_product(gradient_rk_1,gradient_rk_1);
    beta_k=beta_k/fcube_inner_product(gradient_rk,gradient_rk);
    gradient_cg_pk_1=gradient_rk_1+beta_k*gradient_cg_pk;
    //updata
    datatp_k=datatp_k_1;
    gradient_rk=gradient_rk_1;
    gradient_cg_pk=gradient_cg_pk_1;

    while(iter<iterations_num && residual_k>residual_pow){
        iter++;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(gradient_rk))));
        residual_k=sum_num(0,0);

        slantstack3d_recover_L_operator(recoverdatatx_uk,gradient_cg_pk,\
            par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
            par.dt,par.numthread);
        slantstack3d_stack_LT_operator(A_gradient_cg_pk,recoverdatatx_uk,\
            par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
            par.dt,par.numthread);
        //regularization
        fcubemul(datatp_weighted,gradient_cg_pk,par.weighted_fcube);
        A_gradient_cg_pk=dataxs*A_gradient_cg_pk+datatp_weighted;
        
        alpha_k=fcube_inner_product(gradient_rk,gradient_rk);
        alpha_k=alpha_k/fcube_inner_product(gradient_cg_pk,A_gradient_cg_pk);
        //get new solution
        datatp_k_1=datatp_k+alpha_k*gradient_cg_pk;
        gradient_rk_1=gradient_rk-alpha_k*A_gradient_cg_pk;

        beta_k=fcube_inner_product(gradient_rk_1,gradient_rk_1);
        beta_k=beta_k/fcube_inner_product(gradient_rk,gradient_rk);
        gradient_cg_pk_1=gradient_rk_1+beta_k*gradient_cg_pk;
        //updata
        datatp_k=datatp_k_1;
        gradient_rk=gradient_rk_1;
        gradient_cg_pk=gradient_cg_pk_1;
    }
    par.datatp=datatp_k;
    cout<<"iterations times:"<<iter<<" || "<<\
        "err level:"<<residual_k/residual_pow<<endl;

}

// 2-D Slant-stack in the time - space domain.
int Slantstack_CG_2D(float **trace, int ntrace, int ns, float dt,\
 float *coor, float x0, float **tauppanel, int npsr,float psrmin, float dpsr,\
 float **recoverdata, bool dorecover=true, int ncpu=1,  \
 float factor_L2=1.0,float factor_L1=1.0, int iterations_num=48, float residual_ratio=0.5)
{
/* discription.
This function will do the local LRT over local traces .
float **trace;
the INPUT local seismic-gather, i.e. trace[ntrace][ns].
int ntrace;
the trace-number of the INPUT gather.
int ns;
the trace length.
float dt;
time-sampling interval, (unit of dt should be second) ,
float *coor;
coor[ntrace] stores the offset of each trace .
float x0;
the beam-center coordinate .
float **tauppanel; 
the OUTPUT tau-p spectrum, i.e. tauppanel[npsr][ns].
int npsr;
the ray-parameter number of tau-p panel.
float psrmin;
the minimum ray-parameter, (unit is s/m) .
float dpsr;
the interval of ray-parameter， (unit is s/m) .
float **recoverdata;
the OUTPUT recover local seismic-gather, i.e. recoverdata[ntrace][ns].
bool dorecover; default true
Whether to OUTPUT recover local seismic-gather
float L2; default 1.0
L2, Tikhonov regularization parameter
float fmax; default 150
The maximum frequency of the seismic signal
int ncpu; default 1
The number of threads computed by multithreading
int iterations_num=48, 
the Maximum number of iterations of conjugate-gradient method; default 48
float residual_ratio=0.5
Conjugate gradient method iteration error, The smaller the more accurate
*/
///////////////////////////////////////////////////////////////////
    struct slantstack3d par;

    int i,j,k;
    int nz2(ns),nz(ns),nx(ntrace),ny(1),nf(0);
//////////////////////////radon par-set////////////////////////////
    slantstack3d_parset(nx,ny,nz,par);
    par.dp_trace=dpsr;
    par.dp_line=dpsr;
    par.dt=dt;

//The default px of central channel is zero
    par.np_trace=npsr;
    par.ptrace_coord0=psrmin;
    par.np_line=1;
    par.pline_coord0=0.0;   

    //ncpu
    par.numthread=ncpu;

//regularization parameter
    par.factor_l2=nx*ny*factor_L2;  //L2, Tikhonov 
    par.factor_l1=factor_L1;  //L1 

//Parameters updated
    slantstack3d_parupdate(par);
//Seismic trace coordinates
    for(k=0;k<par.ntrace;k++){
        par.ntrace_coordx(k,0)=coor[k]-x0;
        par.nline_coordy(k,0)=0;
    }
//////////////////////////////////////////////////////////////////
//data input, (row,col,slice) of par.data is (nx,ny,nz)
//nx: Spatial sampling nummber
//ny: trace nummber
//nz: time sampling nummber
for(i=0;i<ntrace;i++){
    for(j=0;j<ns;j++){
        par.datatx(i,0,j)=trace[i][j];
    }
}
//////////////////////////////////////////////////
    slantstack3d_stack_multithread(par);
    slantstack3d_stack_CG_invL_operator(par,iterations_num,residual_ratio);
////////////////////output recover data/////////////////////
//output recover data, (row,col,slice) of par.realrebuildtx is (nx,ny,nz)
//npx: X Ray parameters sampling nummber
//npy: Y Ray parameters sampling nummber
//nz: time sampling nummber
    if(dorecover){
        slantstack3d_recover_multithread(par);
        for(i=0;i<ntrace;i++){
            for(j=0;j<ns;j++){
                recoverdata[i][j]=par.recoverdatatx(i,0,j);
            }
        }
    }

///////////////////output tau-p///////////////////
//output tau-p, (row,col,slice) of par.realdataTP is (npx,npy,nz)
//npx: X Ray parameters sampling nummber
//npy: Y Ray parameters sampling nummber
//nz: time sampling nummber
    for(i=0;i<npsr;i++){
        for(j=0;j<ns;j++){
            tauppanel[i][j]=par.datatp(i,0,j);
        }
    }
    return 0;
}
int Slantstack_CG_3D(fcube &tauppanel,fcube &recoverdata,fcube &recovererr,\
 fcube& trace, fmat coordx,fmat coordy, int ns, int ntrace,int nline,float dt,\
 int npx,float pxmin, float dpx,int npy,float pymin, float dpy,\
 int ncpu=1, float factor_L2=0.1,float factor_L1=0.1,int iterations_num=45,\
 float residual_ratio=0.1, bool dorecover=false,bool docginv=true)
{
    struct slantstack3d par;

    int i,j,k;
    int nz2(ns),nz(ns),nx(ntrace),ny(nline),nf(0);
//////////////////////////radon par-set////////////////////////////
    slantstack3d_parset(nx,ny,nz,par);
    par.dp_trace=dpx;
    par.dp_line=dpy;
    par.dt=dt;

//The default px of central channel is zero
    par.np_trace=npx;
    par.ptrace_coord0=pxmin;
    par.np_line=npy;
    par.pline_coord0=pymin;  

//ncpu
    par.numthread=ncpu;

//regularization parameter
    par.factor_l2=factor_L2;  //L2, Tikhonov 
    par.factor_l1=factor_L1;  //L1 
//Parameters updated
    slantstack3d_parupdate(par);
//Seismic trace coordinates
    for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
        par.ntrace_coordx(i,j)=coordx(i,j);
        par.nline_coordy(i,j)=coordy(i,j);
    }}
//////////////////////////////////////////////////
    par.datatx=trace;
    slantstack3d_stack_multithread(par);
    if(docginv){
        slantstack3d_stack_CG_invL_operator\
        (par,iterations_num,residual_ratio);
    }
    tauppanel=par.datatp;

    if(dorecover){
        slantstack3d_recover_multithread(par);
        recoverdata=par.recoverdatatx;
        recovererr=recoverdata-trace;
    }
    
    return 0;
}

#endif
