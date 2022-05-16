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
    fcube datatx,datatp,recoverdatatx;
    fmat weighted_fmat,pline_coord,ptrace_coord,nline_coord,ntrace_coord;
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
    par.ntrace_coord.zeros(par.ntrace,1);
    par.nline_coord.zeros(par.nline,1);
    par.datatp.zeros(par.np_trace,par.np_line,par.nt);
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
    for(k=0;k<par.ntrace;k++)
    {
        par.ntrace_coord(k,0)=k*par.d_trace+par.trace_coord0;
    }
    for(k=0;k<par.nline;k++)
    {
        par.nline_coord(k,0)=k*par.d_line+par.line_coord0;
    }

}

void slantstack3d_cleardata(struct slantstack3d & par)
{
    par.datatp.fill(0.0);
    par.recoverdatatx.fill(0.0);
}

/////////////////////////slantstack3d_L_LT///////////////////////////
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
                dtao1=ntrace_coord(in,0)*ptrace_coord(ip,0);
                for(jn=0;jn<n2;jn++){
                    dtao2=nline_coord(jn,0)*pline_coord(jp,0);
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
                dtao1=ntrace_coord(in,0)*ptrace_coord(ip,0);
                for(jp=0;jp<np2;jp++){
                    dtao2=nline_coord(jn,0)*pline_coord(jp,0);
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

void slantstack3d_stack(struct slantstack3d &par){
    slantstack3d_stack_LT_operator(par.datatp,par.datatx,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
    par.dt);
}
void slantstack3d_recover(struct slantstack3d &par){
    slantstack3d_recover_L_operator(par.recoverdatatx,par.datatp,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
    par.dt);
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

void slantstack3d_stack_CG_invL_operator(struct slantstack3d &par,\
 int iterations_num=9, float residual_ratio=1)
{

    int ip,jp,in,jn,k,iter(0);
    int n1(par.datatx.n_rows),n2(par.datatx.n_cols),n3(par.datatx.n_slices),\
        np1(par.datatp.n_rows),np2(par.datatp.n_cols);
    fcube gradient_rk,gradient_rk_1,\
        gradient_cg_pk,gradient_cg_pk_1,\
        datatp_k,datatp_k_1,\
        recoverdatatx_uk,A_gradient_cg_pk,A_datatp_k;
    fmat sum_num(1,1);
    float beta_k,alpha_k,residual_pow,residual_k;
    datatp_k.copy_size(par.datatp);
    datatp_k_1.copy_size(par.datatp);
    gradient_rk.copy_size(par.datatp);
    gradient_rk_1.copy_size(par.datatp);
    gradient_cg_pk.copy_size(par.datatp);
    gradient_cg_pk_1.copy_size(par.datatp); 
    recoverdatatx_uk.copy_size(par.datatx);
    A_gradient_cg_pk.copy_size(par.datatp);
    A_datatp_k.copy_size(par.datatp);
    sum_num=sum(sum(sum(abs(par.datatp))));
    residual_pow=sum_num(0,0);
    residual_pow*=residual_ratio;

    iter=0;

    datatp_k=par.datatp/par.nline/par.ntrace;
    slantstack3d_recover_L_operator(recoverdatatx_uk,datatp_k,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
    par.dt);
    slantstack3d_stack_LT_operator(A_datatp_k,recoverdatatx_uk,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
    par.dt);
    //regularization
    A_datatp_k=A_datatp_k+datatp_k*par.factor_l2;

    gradient_rk=par.datatp-A_datatp_k;
    gradient_cg_pk=gradient_rk;

    //cal residual_pow
    sum_num=sum(sum(sum(abs(gradient_rk))));
    residual_k=sum_num(0,0);

    slantstack3d_recover_L_operator(recoverdatatx_uk,gradient_cg_pk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
        par.dt);
    slantstack3d_stack_LT_operator(A_gradient_cg_pk,recoverdatatx_uk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
        par.dt);
    //regularization
    A_gradient_cg_pk=A_gradient_cg_pk+gradient_cg_pk*par.factor_l2;
        
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
        cout<<iter<<"||"<<residual_k/residual_pow<<endl;
        iter++;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(gradient_rk))));
        residual_k=sum_num(0,0);

        slantstack3d_recover_L_operator(recoverdatatx_uk,gradient_cg_pk,\
            par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
            par.dt);
        slantstack3d_stack_LT_operator(A_gradient_cg_pk,recoverdatatx_uk,\
            par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
            par.dt);
        //regularization
        A_gradient_cg_pk=A_gradient_cg_pk+gradient_cg_pk*par.factor_l2;
        
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
    cout<<iter<<"||"<<residual_k/residual_pow<<endl;

}

void slantstack3d_stack_HSCG_invL_operator(struct slantstack3d &par,\
 int iterations_num=9, int residual_ratio=0.01)
{

    int ip,jp,in,jn,k,iter(0);
    int n1(par.datatx.n_rows),n2(par.datatx.n_cols),n3(par.datatx.n_slices),\
        np1(par.datatp.n_rows),np2(par.datatp.n_cols);
    fcube gradient_gk,gradient_gk_1,\
        gradient_cg_dk,gradient_cg_dk_1,\
        datatp_k,datatp_k_1,\
        recoverdatatx_uk;
    fmat sum_num(1,1);
    float beta_k,alpha_k,residual_pow,residual_k;
    datatp_k.copy_size(par.datatp);
    datatp_k_1.copy_size(par.datatp);
    gradient_gk.copy_size(par.datatp);
    gradient_gk_1.copy_size(par.datatp);
    gradient_cg_dk.copy_size(par.datatp);
    gradient_cg_dk_1.copy_size(par.datatp); 
    recoverdatatx_uk.copy_size(par.datatx);
    datatp_k_1=par.datatp;

    iter=0;alpha_k=0.00005;

    slantstack3d_recover_L_operator(recoverdatatx_uk,datatp_k_1,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
    par.dt);
    recoverdatatx_uk=recoverdatatx_uk-par.datatx;
    //cal residual_pow
    sum_num=sum(sum(sum(abs(recoverdatatx_uk))));
    residual_pow=sum_num(0,0);
    residual_k=residual_pow;
    residual_pow*=residual_ratio;
    //
    slantstack3d_stack_LT_operator(gradient_gk,recoverdatatx_uk,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
    par.dt);
    //updata radient
    beta_k=1.0;
    gradient_gk_1=gradient_gk;
    //updata radient_cg
    gradient_cg_dk=-gradient_gk;
    gradient_cg_dk_1=gradient_cg_dk;
    //updata solution
    datatp_k=datatp_k_1+alpha_k*gradient_cg_dk;
    datatp_k_1=datatp_k;

    while(iter<iterations_num || residual_k<residual_pow){
        cout<<iter<<"||"<<residual_k<<endl;
        iter++;
        slantstack3d_recover_L_operator(recoverdatatx_uk,datatp_k_1,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
        par.dt);
        recoverdatatx_uk=recoverdatatx_uk-par.datatx;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(recoverdatatx_uk))));
        residual_k=sum_num(0,0);
        //
        slantstack3d_stack_LT_operator(gradient_gk,recoverdatatx_uk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coord,par.nline_coord,\
        par.dt);
        //updata radient
        beta_k=fcube_inner_product(gradient_gk,\
            gradient_gk-gradient_gk_1);
        beta_k=beta_k/fcube_inner_product(gradient_cg_dk_1,\
            gradient_gk-gradient_gk_1);
        gradient_gk_1=gradient_gk;
        //updata radient_cg
        gradient_cg_dk=-gradient_gk+beta_k*gradient_cg_dk_1;
        gradient_cg_dk_1=gradient_cg_dk;
        //updata solution
        datatp_k=datatp_k_1+alpha_k*gradient_cg_dk;
        datatp_k_1=datatp_k;
    }
    par.datatp=datatp_k;

}

#endif
