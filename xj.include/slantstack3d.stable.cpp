/************update in 2022.06.12**************/
/*
    
***********************************************/
#include <iostream>
#include <math.h>
#include "./slantstack3d.hpp"
#include "./armadillo.h"
using namespace std;
using namespace arma;
template<typename T1> void matdelete(T1 **mat, int x1)
{
    int i;
    for(i=0;i<x1;i++)
    {
        delete []mat[i];
        mat[i]=NULL;
    }
    delete []mat;
    mat=NULL;
}

void medianfilter(float **mat, int n1, int n2, int w1=25, int w2=1, float lmd=3.0)
{
    w1=min(w1,(n1/2)-1);
    w2=min(w2,(n2/2)-1);
    fmat data(n1,n2),datapow(n1,n2),datapowfilter(n1,n2);
    int i,j,k1,k2,k3,num((1+2*w1)*(1+2*w2));
    float *win,avera;
    fmat datawin(num,1);

    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            data(i,j)=mat[i][j];
        }
    }
    datapow=abs(data);
    datapowfilter=datapow;
    for(i=w1;i<n1-w1;i++){
    for(j=w2;j<n2-w2;j++){
        k3=0;
    for(k1=-w1;k1<=w1;k1++){
    for(k2=-w2;k2<=w2;k2++){
        datawin(k3,0)=datapow(i+k1,j+k2);
        k3++;
    }}
    avera=sum(sum(datawin))/num;
    if(datapow(i,j)>=lmd*avera)
        datapowfilter(i,j)=avera;
    }}
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            mat[i][j]=mat[i][j]*datapowfilter(i,j)/(datapow(i,j)+0.00001);
        }
    }
}
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
    for(k=0;k<par.np_trace;k++){
        par.ptrace_coord(k,0)=k*par.dp_trace+par.ptrace_coord0;
    }
    for(k=0;k<par.np_line;k++){
        par.pline_coord(k,0)=k*par.dp_line+par.pline_coord0;
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
            datatp(ip,jp,k)=datatp(ip,jp,k)\
                +datatx(in,jn,ntaoceil)*wceil\
                +datatx(in,jn,ntaofloor)*wfloor;
        }
    }}}}}
}
/////////////////////////////////////////////////////////////////////////
void slantstack3d_recover_L_operator(fcube &recoverdatatx,fcube &datatp,\
 fmat &ptrace_coord,fmat &pline_coord,fmat &ntrace_coord,fmat &nline_coord,\
 float dt){

    int ip,jp,in,jn,k,ntaoceil,ntaofloor;
    float tao0,dtao1,dtao2,dtao12,ntao12,wceil,wfloor,wsum,wmin(0.000000001);
    int n1(recoverdatatx.n_rows),n2(recoverdatatx.n_cols),\
        n3(recoverdatatx.n_slices),np1(datatp.n_rows),np2(datatp.n_cols);
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
void slantstack3d_recover(struct slantstack3d &par){
    slantstack3d_recover_L_operator(par.recoverdatatx,par.datatp,\
    par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
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
fmat fmatsmooth(fmat mat2, int x1, int x2, int k)
{
    int i,j,n;
    fmat mat1(x1,x2);
    mat1=mat2;

    for(n=0;n<k;n++){
    for(i=0;i<x1;i++){
    for(j=0;j<x2;j++){
        if(i>=1 && i<x1-1 && j>=1 && j<x2-1){
            mat1(i,j)=((mat2(i-1,j-1)*1.0/12+mat2(i-1,j)*1.0/6+mat2(i-1,j+1)*1.0/12\
                +mat2(i,j-1)*1.0/6+mat2(i,j+1)*1.0/6+mat2(i+1,j-1)*1.0/12\
                +mat2(i+1,j)*1.0/6+mat2(i+1,j+1)*1.0/12+mat2(i,j)*1.0/3))*(3.0/4.0);
        }
        if(i==0)
        {mat1(i,j)=mat1(i+1,j);}
        if(j==0)
        {mat1(i,j)=mat1(i,j+1);}
        if(i==x1-1)
        {mat1(i,j)=mat1(i-1,j);}
        if(j==x2-1)
        {mat1(i,j)=mat1(i,j-1);}
    }}
    mat2=mat1;
    }
    return mat1;
}
void get_weighted_fcube(fcube &filter, fcube &data1,\
 int wx=5, int wy=5, int wz=15)
{
    int i,j,k,i1,j1,nx,ny,nz,k1,k2,win;
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
            agcpow1+=(data1(i+i1,j+j1,k+k1)*data1(i-i1,j+j1,k+k1)*\
                data1(i+i1,j-j1,k+k1)*data1(i-i1,j-j1,k+k1))*\
                (data1(i+i1,j+j1,k-k1)*data1(i-i1,j+j1,k-k1)*\
                data1(i+i1,j-j1,k-k1)*data1(i-i1,j-j1,k-k1));
            //agcpow1+=data1(i+i1,j+j1,k+k1)*data1(i+i1,j+j1,k+k1);
        }}}
        filter(i,j,k)=sqrt(agcpow1);
        agcpow1=0;
    }}}
    filter=(filter)/filter.max();
    for(j=0;j<ny;j++){
        filter.col(j)=fmatsmooth(filter.col(j), nx, nz, 9);
    }
    filter=(filter)/filter.max();
    filter+=0.00001;
    filter=1.0/filter;
    filter=(filter)/filter.max();
}

void slantstack3d_stack_CG_invL_operator(struct slantstack3d &par,\
 int iterations_num=25, float residual_ratio=0.1)
{
    get_weighted_fcube(par.weighted_fcube,par.datatp);
    float dataxs(1.0),datamax;
    par.weighted_fcube*=par.factor_l1;
    par.weighted_fcube+=par.factor_l2;
    //datawrite3d_bycol_transpose(par.weighted_fcube,par.nt,par.np_trace,"weight.bin");

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
    //sum_num=sum(sum(sum(abs(par.datatp))));
    //residual_pow=sum_num(0,0);
    //residual_pow*=residual_ratio;

    iter=0;datamax=0;
    for(k=0;k<par.ntrace;k++){
        datamax+=par.datatx.row(k).max();
    }
    datamax/=par.ntrace;
    datatp_k=datamax*(par.datatp/par.datatp.max());
    slantstack3d_recover_L_operator(recoverdatatx_uk,datatp_k,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt);
    slantstack3d_stack_LT_operator(A_datatp_k,recoverdatatx_uk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt);
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
        par.dt);
    slantstack3d_stack_LT_operator(A_gradient_cg_pk,recoverdatatx_uk,\
        par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
        par.dt);
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
        //cout<<"iterations times:"<<iter<<" || "<<\
        "err level:"<<residual_k/residual_pow<<endl;
        iter++;
        //cal residual_pow
        sum_num=sum(sum(sum(abs(gradient_rk))));
        residual_k=sum_num(0,0);

        slantstack3d_recover_L_operator(recoverdatatx_uk,gradient_cg_pk,\
            par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
            par.dt);
        slantstack3d_stack_LT_operator(A_gradient_cg_pk,recoverdatatx_uk,\
            par.ptrace_coord,par.pline_coord,par.ntrace_coordx,par.nline_coordy,\
            par.dt);
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
void Slantstack_CG_2D(float **recoverdata, float **tauppanel, float **trace,\
 float *coor, int ntrace, float x0, int nt, int ntau, float dt,\
 int npsr,float psrmin, float dpsr,\
 float factor_L2=1.0,float factor_L1=1.0, \
 int iterations_num=48, float residual_ratio=0.5, bool dorecover=false)
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
    int nz(nt),nx(ntrace),ny(1),nf(0);
    ntau=min(ntau,nt);
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

//regularization parameter
    par.factor_l2=factor_L2;  //L2, Tikhonov 
    par.factor_l1=factor_L1;  //L1, Sparsity

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
//ns: time sampling nummber
    for(i=0;i<ntrace;i++){
        for(j=0;j<nt;j++){
            par.datatx(i,0,j)=trace[i][j];
    }}
//////////////////////////////////////////////////
    slantstack3d_stack(par);
    slantstack3d_stack_CG_invL_operator(par,iterations_num,residual_ratio);
////////////////////output recover data/////////////////////
//output recover data, (row,col,slice) of par.realrebuildtx is (nx,ny,ns)
    if(dorecover){
        slantstack3d_recover(par);
        for(i=0;i<ntrace;i++){
            for(j=0;j<nt;j++){
                recoverdata[i][j]=par.recoverdatatx(i,0,j);
        }}
    }
///////////////////output tau-p///////////////////
//output tau-p, (row,col,slice) of par.realdataTP is (npx,npy,nz)
//npx: X Ray parameters sampling nummber
//npy: Y Ray parameters sampling nummber
    for(i=0;i<npsr;i++){
        for(j=0;j<ntau;j++){
            tauppanel[i][j]=par.datatp(i,0,j);
    }}
}
void apply_GlobalLSLRT_to_CMPgather(int ntrCMP, float *offset, int ns, float dt_s,\
 float **gather, int nphr, float phrmin, float dphr, int ntau,float **taup2d_lslrt,\
 float factor_L2=1.0,float factor_L1=1.0, \
 int iterations_num=45, float residual_ratio=0.1)
{
//void Slantstack_CG_2D(float **recoverdata, float **tauppanel, float **trace,\
 float *coor, int ntrace, float x0, int ns, float dt,\
 int npsr,float psrmin, float dpsr,\
 float factor_L2=1.0,float factor_L1=1.0, \
 int iterations_num=48, float residual_ratio=0.5)
 float **recover;
 float x0(0.0);
    Slantstack_CG_2D(recover,taup2d_lslrt,gather,\
            offset,ntrCMP, x0,ns,ntau, dt_s,\
            nphr,phrmin, dphr,\
            factor_L2, factor_L1, \
            iterations_num, residual_ratio,false);
}
float** newfmat(int x1, int x2)
{
    float **p;
    int j;
    p=new float*[x1];      
    for(j=0;j<x1;j++)  
        {  
        p[j]=new float[x2];
        }
    return p;
}

void apply_LocalLSLRT_to_CMPgather(int ntrCMP, float *offset_site, int ns, float dt_s,\
 float **gather, float offsetWidth, int nphr, float phrmin, float dphr, int ntau,\
 float **taup2d_lslrt, int mintrace=25, float factor_L2=1.0,float factor_L1=1.0, \
 int iterations_num=45, float residual_ratio=0.1)
{
    float offset,offset1,offset2,doffset(offsetWidth),minnx(mintrace);
    int i,j,k,ntr_local,trace1,trace2;
    for(i=0;i<nphr;i++){
        for(j=0;j<ntau;j++){
            taup2d_lslrt[i][j]=0.0;
        }
    }
    cout<<"ntrCMP="<<ntrCMP<<endl;
    cout<<"ns="<<ns<<endl;
    cout<<"dt_s="<<dt_s<<endl;
    cout<<"nphr="<<nphr<<endl;
    cout<<"phrmin="<<phrmin<<endl;
    cout<<"dphr="<<dphr<<endl;
    cout<<"ntau="<<ntau<<endl;
    cout<<"offset0="<<offset_site[0]<<endl;
    cout<<"offset2="<<offset_site[2]<<endl;

    float **recover,*coord,**gather_local,**taup_local;
    gather_local=newfmat(ntrCMP, ns);
    taup_local=newfmat(nphr, ntau);
    coord=new float[ntrCMP];
    float x0(0.0);
    offset1=offset_site[0];
    trace1=0;
    for(k=1;k<ntrCMP;k++){
        offset2=offset_site[k];
        trace2=k;
        if(abs(offset2-offset1)>doffset && (trace2-trace1)>mintrace){
            for(i=trace1;i<trace2;i++){
                coord[i-trace1]=offset_site[i];
            for(j=0;j<ns;j++){
                gather_local[i-trace1][j]=gather[i][j];
            }}
            Slantstack_CG_2D(recover,taup_local,gather_local,\
                    coord,(trace2-trace1), x0,ns,ntau, dt_s,\
                    nphr,phrmin, dphr,\
                    factor_L2, factor_L1, \
                    iterations_num, residual_ratio,false);
            for(i=0;i<nphr;i++){
            for(j=0;j<ntau;j++){
                taup2d_lslrt[i][j]+=taup_local[i][j];
            }}
            trace1=trace2;
            offset1=offset2;
        }
        else if(k==(ntrCMP-1)){
            for(i=trace1;i<trace2;i++){
                coord[i-trace1]=offset_site[i];
            for(j=0;j<ns;j++){
                gather_local[i-trace1][j]=gather[i][j];
            }}
            Slantstack_CG_2D(recover,taup_local,gather_local,\
                    coord,(trace2-trace1), x0,ns,ntau, dt_s,\
                    nphr,phrmin, dphr,\
                    factor_L2, factor_L1, \
                    iterations_num, residual_ratio,false);
            for(i=0;i<nphr;i++){
            for(j=0;j<ntau;j++){
                taup2d_lslrt[i][j]+=taup_local[i][j];
            }}
            trace1=trace2;
            offset1=offset2;
        }
    }
    matdelete(gather_local,ntrCMP);
    matdelete(taup_local,nphr);
    delete [] coord;
}
