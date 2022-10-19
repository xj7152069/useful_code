/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef ELASTIC3D_ARMA_HPP
#define ELASTIC3D_ARMA_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <thread>
#include <future>
#include "../xjc.h"
#include <omp.h>

using namespace std;
using namespace arma;
/////////////////////////////////////////////////////////////////////////////////////////////
class elastic3D_ARMA
{ 
private:

public:
//pd_: partial derivative to x/y/z
//mpar: model parameter
    fcube mpar_1_dec_ro,mpar_lmd_add_2miu,mpar_lmd,mpar_miu;
    fcube mpar_1_dec_ro_pdx_jc1,mpar_lmd_add_2miu_pdx_jc1,mpar_lmd_pdx_jc1,mpar_miu_pdx_jc1;
    fcube mpar_1_dec_ro_pdy_jc1,mpar_lmd_add_2miu_pdy_jc1,mpar_lmd_pdy_jc1,mpar_miu_pdy_jc1;
    fcube mpar_1_dec_ro_pdz_jc1,mpar_lmd_add_2miu_pdz_jc1,mpar_lmd_pdz_jc1,mpar_miu_pdz_jc1;
    fcube mpar_ro,mpar_vp,mpar_vs;
    fcube vx,vy,vz,vx_t2,vy_t2,vz_t2,\
        vx_pdx_t1,vx_pdy_t1,vx_pdz_t1,\
        vy_pdx_t1,vy_pdy_t1,vy_pdz_t1,\
        vz_pdx_t1,vz_pdy_t1,vz_pdz_t1,\
        vx_pdx_t2,vx_pdy_t2,vx_pdz_t2,\
        vy_pdx_t2,vy_pdy_t2,vy_pdz_t2,\
        vz_pdx_t2,vz_pdy_t2,vz_pdz_t2;
    fcube txx,tyy,tzz,txx_t2,tyy_t2,tzz_t2,\
        txx_pdx_t1,tyy_pdx_t1,tzz_pdx_t1,\
        txx_pdy_t1,tyy_pdy_t1,tzz_pdy_t1,\
        txx_pdz_t1,tyy_pdz_t1,tzz_pdz_t1,\
        txx_pdx_t2,tyy_pdx_t2,tzz_pdx_t2,\
        txx_pdy_t2,tyy_pdy_t2,tzz_pdy_t2,\
        txx_pdz_t2,tyy_pdz_t2,tzz_pdz_t2;
    fcube txy,txz,tyz,txy_t2,txz_t2,tyz_t2,\
        txy_pdx_t1,txy_pdy_t1,txy_pdx_t2,txy_pdy_t2,\
        txz_pdx_t1,txz_pdz_t1,txz_pdx_t2,txz_pdz_t2,\
        tyz_pdy_t1,tyz_pdz_t1,tyz_pdy_t2,tyz_pdz_t2;

    float dx,dy,dz,dt,PML_wide,R,isPMLSurface;
    float C_X,C_Y,C_Z;
    int nx,ny,nz,nzSampleOfFreeSurface,ompThreadNum;
    
    elastic3D_ARMA();
    elastic3D_ARMA(const int x, const int y, const int z);
    ~elastic3D_ARMA();

    void cleardata();
    void updatepar();
    void calx_p(fcube &ut2 , fcube &ut1, const fcube &u, const fcube &m, const int pjc);
    void caly_p(fcube &ut2 , fcube &ut1, const fcube &u, const fcube &m, const int pjc);
    void calz_p(fcube &ut2 , fcube &ut1, const fcube &u, const fcube &m, const int pjc);
};
///////////////////////////////////////////////////////////////////////////////////////////
elastic3D_ARMA::elastic3D_ARMA()
{
    nx=0;ny=0;nz=0;
    dx=5.0;dy=5.0;dz=5.0;dt=0.0003;
    PML_wide=30.0;isPMLSurface=1.0;R=12.0;
    nzSampleOfFreeSurface=60;
    ompThreadNum=1;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

elastic3D_ARMA::elastic3D_ARMA(const int x, const int y, const int z)
{
    nx=x;ny=y;nz=z;
    dx=5.0;dy=5.0;dz=5.0;dt=0.0003;
    PML_wide=30.0;isPMLSurface=1.0;R=12.0;
    nzSampleOfFreeSurface=60;
    ompThreadNum=1;
    C_Y=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dy/dy/dy;
    C_X=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dx/dx/dx;
    C_Z=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dz/dz/dz;
    
    mpar_1_dec_ro.zeros(nx,ny,nz),mpar_lmd_add_2miu.zeros(nx,ny,nz),\
    mpar_lmd.zeros(nx,ny,nz),mpar_miu.zeros(nx,ny,nz);
    mpar_ro.zeros(nx,ny,nz),mpar_vp.zeros(nx,ny,nz),mpar_vs.zeros(nx,ny,nz);
    mpar_1_dec_ro_pdx_jc1.zeros(nx,ny,nz),mpar_lmd_add_2miu_pdx_jc1.zeros(nx,ny,nz),\
    mpar_lmd_pdx_jc1.zeros(nx,ny,nz),mpar_miu_pdx_jc1.zeros(nx,ny,nz);
    mpar_1_dec_ro_pdy_jc1.zeros(nx,ny,nz),mpar_lmd_add_2miu_pdy_jc1.zeros(nx,ny,nz),\
    mpar_lmd_pdy_jc1.zeros(nx,ny,nz),mpar_miu_pdy_jc1.zeros(nx,ny,nz);
    mpar_1_dec_ro_pdz_jc1.zeros(nx,ny,nz),mpar_lmd_add_2miu_pdz_jc1.zeros(nx,ny,nz),\
    mpar_lmd_pdz_jc1.zeros(nx,ny,nz),mpar_miu_pdz_jc1.zeros(nx,ny,nz);
    vx.zeros(nx,ny,nz),vy.zeros(nx,ny,nz),vz.zeros(nx,ny,nz),\
    vx_pdx_t1.zeros(nx,ny,nz),vx_pdy_t1.zeros(nx,ny,nz),vx_pdz_t1.zeros(nx,ny,nz),\
    vy_pdx_t1.zeros(nx,ny,nz),vy_pdy_t1.zeros(nx,ny,nz),vy_pdz_t1.zeros(nx,ny,nz),\
    vz_pdx_t1.zeros(nx,ny,nz),vz_pdy_t1.zeros(nx,ny,nz),vz_pdz_t1.zeros(nx,ny,nz),\
    vx_pdx_t2.zeros(nx,ny,nz),vx_pdy_t2.zeros(nx,ny,nz),vx_pdz_t2.zeros(nx,ny,nz),\
    vy_pdx_t2.zeros(nx,ny,nz),vy_pdy_t2.zeros(nx,ny,nz),vy_pdz_t2.zeros(nx,ny,nz),\
    vz_pdx_t2.zeros(nx,ny,nz),vz_pdy_t2.zeros(nx,ny,nz),vz_pdz_t2.zeros(nx,ny,nz);
    txx.zeros(nx,ny,nz),tyy.zeros(nx,ny,nz),tzz.zeros(nx,ny,nz),\
    txx_pdx_t1.zeros(nx,ny,nz),tyy_pdx_t1.zeros(nx,ny,nz),tzz_pdx_t1.zeros(nx,ny,nz),\
    txx_pdy_t1.zeros(nx,ny,nz),tyy_pdy_t1.zeros(nx,ny,nz),tzz_pdy_t1.zeros(nx,ny,nz),\
    txx_pdz_t1.zeros(nx,ny,nz),tyy_pdz_t1.zeros(nx,ny,nz),tzz_pdz_t1.zeros(nx,ny,nz),\
    txx_pdx_t2.zeros(nx,ny,nz),tyy_pdx_t2.zeros(nx,ny,nz),tzz_pdx_t2.zeros(nx,ny,nz),\
    txx_pdy_t2.zeros(nx,ny,nz),tyy_pdy_t2.zeros(nx,ny,nz),tzz_pdy_t2.zeros(nx,ny,nz),\
    txx_pdz_t2.zeros(nx,ny,nz),tyy_pdz_t2.zeros(nx,ny,nz),tzz_pdz_t2.zeros(nx,ny,nz);
    txy.zeros(nx,ny,nz),txz.zeros(nx,ny,nz),tyz.zeros(nx,ny,nz),\
    txy_pdx_t1.zeros(nx,ny,nz),txy_pdy_t1.zeros(nx,ny,nz),txy_pdx_t2.zeros(nx,ny,nz),txy_pdy_t2.zeros(nx,ny,nz),\
    txz_pdx_t1.zeros(nx,ny,nz),txz_pdz_t1.zeros(nx,ny,nz),txz_pdx_t2.zeros(nx,ny,nz),txz_pdz_t2.zeros(nx,ny,nz),\
    tyz_pdy_t1.zeros(nx,ny,nz),tyz_pdz_t1.zeros(nx,ny,nz),tyz_pdy_t2.zeros(nx,ny,nz),tyz_pdz_t2.zeros(nx,ny,nz);
}

elastic3D_ARMA::~elastic3D_ARMA()
{
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void elastic3D_ARMA::cleardata()
{
    vx.zeros(nx,ny,nz),vy.zeros(nx,ny,nz),vz.zeros(nx,ny,nz),\
    vx_pdx_t1.zeros(nx,ny,nz),vx_pdy_t1.zeros(nx,ny,nz),vx_pdz_t1.zeros(nx,ny,nz),\
    vy_pdx_t1.zeros(nx,ny,nz),vy_pdy_t1.zeros(nx,ny,nz),vy_pdz_t1.zeros(nx,ny,nz),\
    vz_pdx_t1.zeros(nx,ny,nz),vz_pdy_t1.zeros(nx,ny,nz),vz_pdz_t1.zeros(nx,ny,nz),\
    vx_pdx_t2.zeros(nx,ny,nz),vx_pdy_t2.zeros(nx,ny,nz),vx_pdz_t2.zeros(nx,ny,nz),\
    vy_pdx_t2.zeros(nx,ny,nz),vy_pdy_t2.zeros(nx,ny,nz),vy_pdz_t2.zeros(nx,ny,nz),\
    vz_pdx_t2.zeros(nx,ny,nz),vz_pdy_t2.zeros(nx,ny,nz),vz_pdz_t2.zeros(nx,ny,nz);
    txx.zeros(nx,ny,nz),tyy.zeros(nx,ny,nz),tzz.zeros(nx,ny,nz),\
    txx_pdx_t1.zeros(nx,ny,nz),tyy_pdx_t1.zeros(nx,ny,nz),tzz_pdx_t1.zeros(nx,ny,nz),\
    txx_pdy_t1.zeros(nx,ny,nz),tyy_pdy_t1.zeros(nx,ny,nz),tzz_pdy_t1.zeros(nx,ny,nz),\
    txx_pdz_t1.zeros(nx,ny,nz),tyy_pdz_t1.zeros(nx,ny,nz),tzz_pdz_t1.zeros(nx,ny,nz),\
    txx_pdx_t2.zeros(nx,ny,nz),tyy_pdx_t2.zeros(nx,ny,nz),tzz_pdx_t2.zeros(nx,ny,nz),\
    txx_pdy_t2.zeros(nx,ny,nz),tyy_pdy_t2.zeros(nx,ny,nz),tzz_pdy_t2.zeros(nx,ny,nz),\
    txx_pdz_t2.zeros(nx,ny,nz),tyy_pdz_t2.zeros(nx,ny,nz),tzz_pdz_t2.zeros(nx,ny,nz);
    txy.zeros(nx,ny,nz),txz.zeros(nx,ny,nz),tyz.zeros(nx,ny,nz),\
    txy_pdx_t1.zeros(nx,ny,nz),txy_pdy_t1.zeros(nx,ny,nz),txy_pdx_t2.zeros(nx,ny,nz),txy_pdy_t2.zeros(nx,ny,nz),\
    txz_pdx_t1.zeros(nx,ny,nz),txz_pdz_t1.zeros(nx,ny,nz),txz_pdx_t2.zeros(nx,ny,nz),txz_pdz_t2.zeros(nx,ny,nz),\
    tyz_pdy_t1.zeros(nx,ny,nz),tyz_pdz_t1.zeros(nx,ny,nz),tyz_pdy_t2.zeros(nx,ny,nz),tyz_pdz_t2.zeros(nx,ny,nz);
    cout<<"All data has cleaar!"<<endl;
}
 
void elastic3D_ARMA::updatepar()
{
    int i,j,k;
    C_Y=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dy/dy/dy;
    C_X=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dx/dx/dx;
    C_Z=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dz/dz/dz;    
    for(i=0;i<nx;i++){ 
    for(j=0;j<ny;j++){
    for(k=0;k<nz;k++){
        mpar_miu(i,j,k)=mpar_vs(i,j,k)*mpar_vs(i,j,k)*mpar_ro(i,j,k);
        mpar_lmd(i,j,k)=mpar_ro(i,j,k)*(mpar_vp(i,j,k)*mpar_vp(i,j,k)\
            -2.0*mpar_vs(i,j,k)*mpar_vs(i,j,k));
        mpar_lmd_add_2miu(i,j,k)=mpar_lmd(i,j,k)+2.0*mpar_miu(i,j,k);
        mpar_1_dec_ro(i,j,k)=1.0/mpar_ro(i,j,k);
        //mo1[i][j]=(lmd[i][j]+2*miu[i][j])/(2*lmd[i][j]+2*miu[i][j])/ro[i][j];
    }}}
    for(i=0;i<nx-1;i++){ 
    for(j=0;j<ny-1;j++){
    for(k=0;k<nz-1;k++){
        mpar_1_dec_ro_pdx_jc1(i,j,k)=(mpar_1_dec_ro(i+1,j,k)+mpar_1_dec_ro(i,j,k))*0.5;
        mpar_1_dec_ro_pdy_jc1(i,j,k)=(mpar_1_dec_ro(i,j+1,k)+mpar_1_dec_ro(i,j,k))*0.5;
        mpar_1_dec_ro_pdz_jc1(i,j,k)=(mpar_1_dec_ro(i,j,k+1)+mpar_1_dec_ro(i,j,k))*0.5;
        mpar_lmd_add_2miu_pdx_jc1(i,j,k)=(mpar_lmd_add_2miu(i+1,j,k)+mpar_lmd_add_2miu(i,j,k))*0.5;
        mpar_lmd_add_2miu_pdy_jc1(i,j,k)=(mpar_lmd_add_2miu(i,j+1,k)+mpar_lmd_add_2miu(i,j,k))*0.5;
        mpar_lmd_add_2miu_pdz_jc1(i,j,k)=(mpar_lmd_add_2miu(i,j,k+1)+mpar_lmd_add_2miu(i,j,k))*0.5;
        mpar_lmd_pdx_jc1(i,j,k)=(mpar_lmd(i+1,j,k)+mpar_lmd(i,j,k))*0.5;
        mpar_lmd_pdy_jc1(i,j,k)=(mpar_lmd(i,j+1,k)+mpar_lmd(i,j,k))*0.5;
        mpar_lmd_pdz_jc1(i,j,k)=(mpar_lmd(i,j,k+1)+mpar_lmd(i,j,k))*0.5;
        mpar_miu_pdx_jc1(i,j,k)=(mpar_miu(i+1,j,k)+mpar_miu(i,j,k))*0.5;
        mpar_miu_pdy_jc1(i,j,k)=(mpar_miu(i,j+1,k)+mpar_miu(i,j,k))*0.5;
        mpar_miu_pdz_jc1(i,j,k)=(mpar_miu(i,j,k+1)+mpar_miu(i,j,k))*0.5;
        //mo1[i][j]=(lmd[i][j]+2*miu[i][j])/(2*lmd[i][j]+2*miu[i][j])/ro[i][j];
    }}}
}

void calx_3d(fcube *ut2 ,fcube *ut1,fcube *u, fcube *m,int pjc,\
    class elastic3D_ARMA *obj)
{
    obj->calx_p(ut2[0],ut1[0],u[0],m[0],pjc);
}
void elastic3D_ARMA::calx_p(fcube &ut2 , fcube &ut1, \
    const fcube &u, const fcube &m, const int pjc)
{
    /*
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数
*/
    int Z,k;
    Z=this->nz;
omp_set_num_threads(this->ompThreadNum);
#pragma omp parallel for
    for(k=8;k<Z-8;k++){
        float C1      = 1.2340911;
        float C2      = -1.0664985e-01;
        float C3      = 2.3036367e-02;
        float C4      = -5.3423856e-03;
        float C5      = 1.0772712e-03;
        float C6      = -1.6641888e-04;
        float C7      = 1.7021711e-005;
        float C8      = -8.5234642e-007;//差分系数
        int i,j,jc(pjc);
        float du1(0),du2(0),du(0);
        float DX,DT,xshd,DT1,DX1;
        int X,Y;
        DX=this->dx,DT=this->dt,xshd=this->PML_wide;
        X=this->nx,Y=this->ny,DT1=1.0/DT,DX1=1.0/DX;
    for(j=8;j<Y-8;j++){
        for(i=8;i<xshd+8;i++){  
            du=(
                C1*(u(i+0+jc,j,k)-u(i-1+jc,j,k))+\
                C2*(u(i+1+jc,j,k)-u(i-2+jc,j,k))+\
                C3*(u(i+2+jc,j,k)-u(i-3+jc,j,k))+\
                C4*(u(i+3+jc,j,k)-u(i-4+jc,j,k))+\
                C5*(u(i+4+jc,j,k)-u(i-5+jc,j,k))+\
                C6*(u(i+5+jc,j,k)-u(i-6+jc,j,k))+\
                C7*(u(i+6+jc,j,k)-u(i-7+jc,j,k))+\
                C8*(u(i+7+jc,j,k)-u(i-8+jc,j,k))\
                );
            du=m(i,j,k)*du*DX1;
            du1=C_X*3000.0*DX*(xshd+8-i)*DX*(xshd+8-i);
            //du1=C_X*this->vp(i,j,k)*DX*(xshd+8-i)*DX*(xshd+8-i);
            ut2(i,j,k)=((du+ut1(i,j,k)*(DT1-du1*0.5))/(DT1+du1*0.5));
        }
        for(i=xshd+8;i<X-xshd-8;i++)
        {  
            du=(
                C1*(u(i+0+jc,j,k)-u(i-1+jc,j,k))+\
                C2*(u(i+1+jc,j,k)-u(i-2+jc,j,k))+\
                C3*(u(i+2+jc,j,k)-u(i-3+jc,j,k))+\
                C4*(u(i+3+jc,j,k)-u(i-4+jc,j,k))+\
                C5*(u(i+4+jc,j,k)-u(i-5+jc,j,k))+\
                C6*(u(i+5+jc,j,k)-u(i-6+jc,j,k))+\
                C7*(u(i+6+jc,j,k)-u(i-7+jc,j,k))+\
                C8*(u(i+7+jc,j,k)-u(i-8+jc,j,k))\
                );
            du=m(i,j,k)*du*DX1;
            ut2(i,j,k)=((du)*DT+ut1(i,j,k));
        }
        for(i=X-xshd-8;i<X-8;i++)
        {  
            du=(
                C1*(u(i+0+jc,j,k)-u(i-1+jc,j,k))+\
                C2*(u(i+1+jc,j,k)-u(i-2+jc,j,k))+\
                C3*(u(i+2+jc,j,k)-u(i-3+jc,j,k))+\
                C4*(u(i+3+jc,j,k)-u(i-4+jc,j,k))+\
                C5*(u(i+4+jc,j,k)-u(i-5+jc,j,k))+\
                C6*(u(i+5+jc,j,k)-u(i-6+jc,j,k))+\
                C7*(u(i+6+jc,j,k)-u(i-7+jc,j,k))+\
                C8*(u(i+7+jc,j,k)-u(i-8+jc,j,k))\
                );
            du=0.5*(m(i+jc,j,k)+m(i,j,k))*du*DX1;
            du1=C_X*3000.0*DX*(i-X+xshd+8)*DX*(i-X+xshd+8);
            //du1=C_X*this->vp(i,j,k)*DX*(i-X+xshd+8)*DX*(i-X+xshd+8);
            ut2(i,j,k)=((du+ut1(i,j,k)*(DT1-du1*0.5))/(DT1+du1*0.5));
        }
    }}
    ut1.swap(ut2);
}

void caly_3d(fcube *ut2 ,fcube *ut1,fcube *u, fcube *m,int pjc,\
    class elastic3D_ARMA *obj)
{
    obj->caly_p(ut2[0],ut1[0],u[0],m[0],pjc);
}
void elastic3D_ARMA::caly_p(fcube &ut2 , fcube &ut1, \
    const fcube &u, const fcube &m, const int pjc)
{
    /*
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数
*/
    int X,i;
    X=this->nx;
omp_set_num_threads(this->ompThreadNum);
#pragma omp parallel for
    for(i=8;i<X-8;i++){
        float DY,DT,xshd,DY1,DT1;
        int Y,Z;
        int j,k,jc(pjc);
        float du1(0),du2(0),du(0);
        DY=this->dy,DT=this->dt,xshd=this->PML_wide;
        Y=this->ny,Z=this->nz,DY1=1.0/DY,DT1=1.0/DT;
        
        float C1      = 1.2340911;
        float C2      = -1.0664985e-01;
        float C3      = 2.3036367e-02;
        float C4      = -5.3423856e-03;
        float C5      = 1.0772712e-03;
        float C6      = -1.6641888e-04;
        float C7      = 1.7021711e-005;
        float C8      = -8.5234642e-007;//差分系数
    for(k=8;k<Z-8;k++){
        for(j=8;j<xshd+8;j++){
            du=(
                C1*(u(i,j+0+jc,k)-u(i,j-1+jc,k))+\
                C2*(u(i,j+1+jc,k)-u(i,j-2+jc,k))+\
                C3*(u(i,j+2+jc,k)-u(i,j-3+jc,k))+\
                C4*(u(i,j+3+jc,k)-u(i,j-4+jc,k))+\
                C5*(u(i,j+4+jc,k)-u(i,j-5+jc,k))+\
                C6*(u(i,j+5+jc,k)-u(i,j-6+jc,k))+\
                C7*(u(i,j+6+jc,k)-u(i,j-7+jc,k))+\
                C8*(u(i,j+7+jc,k)-u(i,j-8+jc,k))\
                );
            du=m(i,j,k)*du*DY1;
            du1=C_Y*3000.0*DY*(xshd+8-j)*DY*(xshd+8-j);
            //du1=C_X*this->vp(i,j,k)*DX*(xshd+8-i)*DX*(xshd+8-i);
            ut2(i,j,k)=((du+ut1(i,j,k)*(DT1-du1*0.5))/(DT1+du1*0.5));
        }
        for(j=xshd+8;j<Y-xshd-8;j++)
        {
            du=(
                C1*(u(i,j+0+jc,k)-u(i,j-1+jc,k))+\
                C2*(u(i,j+1+jc,k)-u(i,j-2+jc,k))+\
                C3*(u(i,j+2+jc,k)-u(i,j-3+jc,k))+\
                C4*(u(i,j+3+jc,k)-u(i,j-4+jc,k))+\
                C5*(u(i,j+4+jc,k)-u(i,j-5+jc,k))+\
                C6*(u(i,j+5+jc,k)-u(i,j-6+jc,k))+\
                C7*(u(i,j+6+jc,k)-u(i,j-7+jc,k))+\
                C8*(u(i,j+7+jc,k)-u(i,j-8+jc,k))\
                );
            du=m(i,j,k)*du*DY1;
            ut2(i,j,k)=(du*DT+ut1(i,j,k));
        }
        for(j=Y-xshd-8;j<Y-8;j++)
        {
            du=(
                C1*(u(i,j+0+jc,k)-u(i,j-1+jc,k))+\
                C2*(u(i,j+1+jc,k)-u(i,j-2+jc,k))+\
                C3*(u(i,j+2+jc,k)-u(i,j-3+jc,k))+\
                C4*(u(i,j+3+jc,k)-u(i,j-4+jc,k))+\
                C5*(u(i,j+4+jc,k)-u(i,j-5+jc,k))+\
                C6*(u(i,j+5+jc,k)-u(i,j-6+jc,k))+\
                C7*(u(i,j+6+jc,k)-u(i,j-7+jc,k))+\
                C8*(u(i,j+7+jc,k)-u(i,j-8+jc,k))\
                );
            du=m(i,j,k)*du*DY1;
            du1=C_Y*3000.0*DY*(j-Y+xshd+8)*DY*(j-Y+xshd+8);
            //du1=C_Y*this->vp(i,j,k)*DY*(j-Y+xshd+8)*DY*(j-Y+xshd+8);
            ut2(i,j,k)=((du+ut1(i,j,k)*(DT1-du1*0.5))/(DT1+du1*0.5));
        }
    }}
    ut1.swap(ut2);
}

void calz_3d(fcube *ut2 ,fcube *ut1,fcube *u, fcube *m,int pjc,\
    class elastic3D_ARMA *obj)
{
    obj->calz_p(ut2[0],ut1[0],u[0],m[0],pjc);
}
void elastic3D_ARMA::calz_p(fcube &ut2 , fcube &ut1,\
    const fcube &u, const fcube &m, const int pjc)
{
    /*
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数
*/
    int X,i;
    X=this->nx;
omp_set_num_threads(this->ompThreadNum);
#pragma omp parallel for
    for(i=8;i<X-8;i++){
        float DZ,DT,xshd,suface_PML,DZ1,DT1;
        int Y,Z,nSurface;
        int j,k,jc(pjc);
        float du1(0),du2(0),du(0);
        DZ=this->dz,DT=this->dt,xshd=this->PML_wide;
        DZ1=1.0/DZ,DT1=1.0/DT;
        Y=this->ny,Z=this->nz,suface_PML=this->isPMLSurface;
        nSurface=this->nzSampleOfFreeSurface;
        float C1      = 1.2340911;
        float C2      = -1.0664985e-01;
        float C3      = 2.3036367e-02;
        float C4      = -5.3423856e-03;
        float C5      = 1.0772712e-03;
        float C6      = -1.6641888e-04;
        float C7      = 1.7021711e-005;
        float C8      = -8.5234642e-007;//差分系数
    for(j=8;j<Y-8;j++){
        for(k=8;k<xshd+8;k++){  
            du=(
                C1*(u(i,j,k+0+jc)-u(i,j,k-1+jc))+\
                C2*(u(i,j,k+1+jc)-u(i,j,k-2+jc))+\
                C3*(u(i,j,k+2+jc)-u(i,j,k-3+jc))+\
                C4*(u(i,j,k+3+jc)-u(i,j,k-4+jc))+\
                C5*(u(i,j,k+4+jc)-u(i,j,k-5+jc))+\
                C6*(u(i,j,k+5+jc)-u(i,j,k-6+jc))+\
                C7*(u(i,j,k+6+jc)-u(i,j,k-7+jc))+\
                C8*(u(i,j,k+7+jc)-u(i,j,k-8+jc))\
                );
            //take care of the FreeSurface_PML_boundary!!!
            //the freesurface should be: ro=1000.0, vp=0.0
            du=m(i,j,k)*du*DZ1;
            du1=C_Z*3000.0*DZ*suface_PML*(xshd+8-k)*DZ\
                *suface_PML*(xshd+8-k);
            //du1=C_Y*this->vp[j][i]*DY*suface_PML*(xshd+8-j)*DY\
                *suface_PML*(xshd+8-j);
            ut2(i,j,k)=((du+ut1(i,j,k)*(DT1-du1*0.5))/(DT1+du1*0.5));
        }
        for(k=xshd+8;k<Z-xshd-8;k++)
        {  
            du=(
                C1*(u(i,j,k+0+jc)-u(i,j,k-1+jc))+\
                C2*(u(i,j,k+1+jc)-u(i,j,k-2+jc))+\
                C3*(u(i,j,k+2+jc)-u(i,j,k-3+jc))+\
                C4*(u(i,j,k+3+jc)-u(i,j,k-4+jc))+\
                C5*(u(i,j,k+4+jc)-u(i,j,k-5+jc))+\
                C6*(u(i,j,k+5+jc)-u(i,j,k-6+jc))+\
                C7*(u(i,j,k+6+jc)-u(i,j,k-7+jc))+\
                C8*(u(i,j,k+7+jc)-u(i,j,k-8+jc))\
                );
            //take care of the FreeSurface_PML_boundary!!!
            //the freesurface should be: ro=1000.0, vp=0.0
            du=m(i,j,k)*du*DZ1;
            ut2(i,j,k)=((du)*DT+ut1(i,j,k));
        }
        for(k=Z-xshd-8;k<Z-8;k++)
        {  
            du=(
                C1*(u(i,j,k+0+jc)-u(i,j,k-1+jc))+\
                C2*(u(i,j,k+1+jc)-u(i,j,k-2+jc))+\
                C3*(u(i,j,k+2+jc)-u(i,j,k-3+jc))+\
                C4*(u(i,j,k+3+jc)-u(i,j,k-4+jc))+\
                C5*(u(i,j,k+4+jc)-u(i,j,k-5+jc))+\
                C6*(u(i,j,k+5+jc)-u(i,j,k-6+jc))+\
                C7*(u(i,j,k+6+jc)-u(i,j,k-7+jc))+\
                C8*(u(i,j,k+7+jc)-u(i,j,k-8+jc))\
                );
            du=m(i,j,k)*du*DZ1;
            du1=C_Z*3000.0*DZ*(k-Z+xshd+8)*DZ*(k-Z+xshd+8);
            //du1=C_Z*this->vp(i,j,k)*DZ*(k-Z+xshd+8)*DX*(k-Z+xshd+8);
            ut2(i,j,k)=((du+ut1(i,j,k)*(DT1-du1*0.5))/(DT1+du1*0.5));
        }
    }}
    ut1.swap(ut2);
}

void TimeSliceCal_elastic3D_ARMA_MultiThread(class elastic3D_ARMA & obj)
{
    int k,jc0(0),jc1(1);
    thread *useThread;
    useThread=new thread[24];

    useThread[0]=thread(calx_3d,&obj.vx_pdx_t2,&obj.vx_pdx_t1,&obj.txx,&obj.mpar_1_dec_ro_pdx_jc1,jc1,&obj);
    useThread[1]=thread(caly_3d,&obj.vx_pdy_t2,&obj.vx_pdy_t1,&obj.txy,&obj.mpar_1_dec_ro,jc0,&obj);
    useThread[2]=thread(calz_3d,&obj.vx_pdz_t2,&obj.vx_pdz_t1,&obj.txz,&obj.mpar_1_dec_ro,jc0,&obj);
    useThread[3]=thread(calx_3d,&obj.vy_pdx_t2,&obj.vy_pdx_t1,&obj.txy,&obj.mpar_1_dec_ro,jc0,&obj);
    useThread[4]=thread(caly_3d,&obj.vy_pdy_t2,&obj.vy_pdy_t1,&obj.tyy,&obj.mpar_1_dec_ro_pdy_jc1,jc1,&obj);
    useThread[5]=thread(calz_3d,&obj.vy_pdz_t2,&obj.vy_pdz_t1,&obj.tyz,&obj.mpar_1_dec_ro,jc0,&obj);
    useThread[6]=thread(calx_3d,&obj.vz_pdx_t2,&obj.vz_pdx_t1,&obj.txz,&obj.mpar_1_dec_ro,jc0,&obj);
    useThread[7]=thread(caly_3d,&obj.vz_pdy_t2,&obj.vz_pdy_t1,&obj.tyz,&obj.mpar_1_dec_ro,jc0,&obj);
    useThread[8]=thread(calz_3d,&obj.vz_pdz_t2,&obj.vz_pdz_t1,&obj.tzz,&obj.mpar_1_dec_ro_pdz_jc1,jc1,&obj);
    for(k=0;k<9;k++){
        if(useThread[k].joinable())
            useThread[k].join();
    }
    obj.vx=obj.vx_pdx_t1+obj.vx_pdy_t1+obj.vx_pdz_t1;
    obj.vy=obj.vy_pdx_t1+obj.vy_pdy_t1+obj.vy_pdz_t1;
    obj.vz=obj.vz_pdx_t1+obj.vz_pdy_t1+obj.vz_pdz_t1;

    useThread[9]=thread(calx_3d,&obj.txx_pdx_t2,&obj.txx_pdx_t1,&obj.vx,&obj.mpar_lmd_add_2miu,jc0,&obj);
    useThread[10]=thread(caly_3d,&obj.txx_pdy_t2,&obj.txx_pdy_t1,&obj.vy,&obj.mpar_lmd,jc0,&obj);
    useThread[11]=thread(calz_3d,&obj.txx_pdz_t2,&obj.txx_pdz_t1,&obj.vz,&obj.mpar_lmd,jc0,&obj);
    useThread[12]=thread(caly_3d,&obj.tyy_pdy_t2,&obj.tyy_pdy_t1,&obj.vy,&obj.mpar_lmd_add_2miu,jc0,&obj);
    useThread[13]=thread(calx_3d,&obj.tyy_pdx_t2,&obj.tyy_pdx_t1,&obj.vx,&obj.mpar_lmd,jc0,&obj);
    useThread[14]=thread(calz_3d,&obj.tyy_pdz_t2,&obj.tyy_pdz_t1,&obj.vz,&obj.mpar_lmd,jc0,&obj);
    useThread[15]=thread(calz_3d,&obj.tzz_pdz_t2,&obj.tzz_pdz_t1,&obj.vz,&obj.mpar_lmd_add_2miu,jc0,&obj);
    useThread[16]=thread(calx_3d,&obj.tzz_pdx_t2,&obj.tzz_pdx_t1,&obj.vx,&obj.mpar_lmd,jc0,&obj);
    useThread[17]=thread(caly_3d,&obj.tzz_pdy_t2,&obj.tzz_pdy_t1,&obj.vy,&obj.mpar_lmd,jc0,&obj);
    for(k=9;k<18;k++){
        if(useThread[k].joinable())
            useThread[k].join();
    }
    obj.txx=obj.txx_pdx_t1+obj.txx_pdy_t1+obj.txx_pdz_t1;
    obj.tyy=obj.tyy_pdx_t1+obj.tyy_pdy_t1+obj.tyy_pdz_t1;
    obj.tzz=obj.tzz_pdx_t1+obj.tzz_pdy_t1+obj.tzz_pdz_t1;

    useThread[18]=thread(caly_3d,&obj.txy_pdy_t2,&obj.txy_pdy_t1,&obj.vx,&obj.mpar_miu_pdy_jc1,jc1,&obj);
    useThread[19]=thread(calx_3d,&obj.txy_pdx_t2,&obj.txy_pdx_t1,&obj.vy,&obj.mpar_miu_pdx_jc1,jc1,&obj);
    useThread[20]=thread(calz_3d,&obj.txz_pdz_t2,&obj.txz_pdz_t1,&obj.vx,&obj.mpar_miu_pdz_jc1,jc1,&obj);
    useThread[21]=thread(calx_3d,&obj.txz_pdx_t2,&obj.txz_pdx_t1,&obj.vz,&obj.mpar_miu_pdx_jc1,jc1,&obj);
    useThread[22]=thread(calz_3d,&obj.tyz_pdz_t2,&obj.tyz_pdz_t1,&obj.vy,&obj.mpar_miu_pdz_jc1,jc1,&obj);
    useThread[23]=thread(caly_3d,&obj.tyz_pdy_t2,&obj.tyz_pdy_t1,&obj.vz,&obj.mpar_miu_pdy_jc1,jc1,&obj);
    for(k=18;k<24;k++){
        if(useThread[k].joinable())
            useThread[k].join();
    }
    obj.txy=obj.txy_pdy_t1+obj.txy_pdx_t1;
    obj.txz=obj.txz_pdz_t1+obj.txz_pdx_t1;
    obj.tyz=obj.tyz_pdz_t1+obj.tyz_pdy_t1;
    delete [] useThread;
}

#endif

