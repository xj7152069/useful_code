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
    fcube mpar_ro,mpar_vp,mpar_vs;
    fcube vx_t1,vy_t1,vz_t1,vx_t2,vy_t2,vz_t2,\
        vx_pdx,vx_pdy,vx_pdz,vy_pdx,vy_pdy,vy_pdz,vz_pdx,vz_pdy,vz_pdz;
    fcube txx_t1,tyy_t1,tzz_t1,txx_t2,tyy_t2,tzz_t2,\
        txx_pdx,tyy_pdy,tzz_pdz;
    fcube txy_t1,txz_t1,tyz_t1,txy_t2,txz_t2,tyz_t2,\
        txy_pdx,txy_pdy,txz_pdx,txz_pdz,tyz_pdy,tyz_pdz;
    fcube data3d1,data3d2,data3d3,data3d4,data3d5,data3d6;  

    float dx,dy,dz,dt,PML_wide,R,isFreeSurface;
    float C_X,C_Y,C_Z;
    int nx,ny,nz,nzSampleOfFreeSurface,ompThreadNum;
    
    elastic3D_ARMA();
    elastic3D_ARMA(const int x, const int y, const int z);
    ~elastic3D_ARMA();

    void cleardata();
    void updatepar();
    void prepareForMultiThread();
    void calx_p(fcube &ut2 , const fcube &ut1, const fcube &u, const fcube &m, const int pjc);
    void caly_p(fcube &ut2 , const fcube &ut1, const fcube &u, const fcube &m, const int pjc);
    void calz_p(fcube &ut2 , const fcube &ut1, const fcube &u, const fcube &m, const int pjc);
};
///////////////////////////////////////////////////////////////////////////////////////////
void elastic3D_ARMA::prepareForMultiThread(){
    this->data3d1.set_size(nx,ny,nz);
    this->data3d2.set_size(nx,ny,nz);
    this->data3d3.set_size(nx,ny,nz);
    this->data3d4.set_size(nx,ny,nz);
    this->data3d5.set_size(nx,ny,nz);
    this->data3d6.set_size(nx,ny,nz);
}
elastic3D_ARMA::elastic3D_ARMA()
{
    nx=0;ny=0;nz=0;
    dx=5.0;dy=5.0;dz=5.0;dt=0.0003;
    PML_wide=40.0;isFreeSurface=1.0;R=25.0;
    nzSampleOfFreeSurface=50;
    ompThreadNum=1;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

elastic3D_ARMA::elastic3D_ARMA(const int x, const int y, const int z)
{
    nx=x;ny=y;nz=z;
    dx=5.0;dy=5.0;dz=5.0;dt=0.0005;
    PML_wide=40.0;isFreeSurface=1.0;R=25.0;
    nzSampleOfFreeSurface=50;
    ompThreadNum=1;
    C_Y=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dy/dy/dy;
    C_X=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dx/dx/dx;
    C_Z=(R)*3.0/2.0/(PML_wide)/(PML_wide)/(PML_wide)/dz/dz/dz;
    
    mpar_1_dec_ro.zeros(nx,ny,nz),mpar_lmd_add_2miu.zeros(nx,ny,nz),\
    mpar_lmd.zeros(nx,ny,nz),mpar_miu.zeros(nx,ny,nz);
    mpar_ro.zeros(nx,ny,nz),mpar_vp.zeros(nx,ny,nz),mpar_vs.zeros(nx,ny,nz);
    vx_t1.zeros(nx,ny,nz),vy_t1.zeros(nx,ny,nz),vz_t1.zeros(nx,ny,nz),\
    vx_t2.zeros(nx,ny,nz),vy_t2.zeros(nx,ny,nz),vz_t2.zeros(nx,ny,nz),\
    vx_pdx.zeros(nx,ny,nz),vx_pdy.zeros(nx,ny,nz),vx_pdz.zeros(nx,ny,nz),\
    vy_pdx.zeros(nx,ny,nz),vy_pdy.zeros(nx,ny,nz),vy_pdz.zeros(nx,ny,nz),\
    vz_pdx.zeros(nx,ny,nz),vz_pdy.zeros(nx,ny,nz),vz_pdz.zeros(nx,ny,nz);
    txx_t1.zeros(nx,ny,nz),tyy_t1.zeros(nx,ny,nz),tzz_t1.zeros(nx,ny,nz),\
    txx_t2.zeros(nx,ny,nz),tyy_t2.zeros(nx,ny,nz),tzz_t2.zeros(nx,ny,nz),\
    txx_pdx.zeros(nx,ny,nz),tyy_pdy.zeros(nx,ny,nz),tzz_pdz.zeros(nx,ny,nz);
    txy_t1.zeros(nx,ny,nz),txz_t1.zeros(nx,ny,nz),tyz_t1.zeros(nx,ny,nz),\
    txy_t2.zeros(nx,ny,nz),txz_t2.zeros(nx,ny,nz),tyz_t2.zeros(nx,ny,nz),\
    txy_pdx.zeros(nx,ny,nz),txy_pdy.zeros(nx,ny,nz),txz_pdx.zeros(nx,ny,nz),\
    txz_pdz.zeros(nx,ny,nz),tyz_pdy.zeros(nx,ny,nz),tyz_pdz.zeros(nx,ny,nz);
}

elastic3D_ARMA::~elastic3D_ARMA()
{
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void elastic3D_ARMA::cleardata()
{
    mpar_1_dec_ro.zeros(nx,ny,nz),mpar_lmd_add_2miu.zeros(nx,ny,nz),\
    mpar_lmd.zeros(nx,ny,nz),mpar_miu.zeros(nx,ny,nz);
    mpar_ro.zeros(nx,ny,nz),mpar_vp.zeros(nx,ny,nz),mpar_vs.zeros(nx,ny,nz);
    vx_t1.zeros(nx,ny,nz),vy_t1.zeros(nx,ny,nz),vz_t1.zeros(nx,ny,nz),\
    vx_t2.zeros(nx,ny,nz),vy_t2.zeros(nx,ny,nz),vz_t2.zeros(nx,ny,nz),\
    vx_pdx.zeros(nx,ny,nz),vx_pdy.zeros(nx,ny,nz),vx_pdz.zeros(nx,ny,nz),\
    vy_pdx.zeros(nx,ny,nz),vy_pdy.zeros(nx,ny,nz),vy_pdz.zeros(nx,ny,nz),\
    vz_pdx.zeros(nx,ny,nz),vz_pdy.zeros(nx,ny,nz),vz_pdz.zeros(nx,ny,nz);
    txx_t1.zeros(nx,ny,nz),tyy_t1.zeros(nx,ny,nz),tzz_t1.zeros(nx,ny,nz),\
    txx_t2.zeros(nx,ny,nz),tyy_t2.zeros(nx,ny,nz),tzz_t2.zeros(nx,ny,nz),\
    txx_pdx.zeros(nx,ny,nz),tyy_pdy.zeros(nx,ny,nz),tzz_pdz.zeros(nx,ny,nz);
    txy_t1.zeros(nx,ny,nz),txz_t1.zeros(nx,ny,nz),tyz_t1.zeros(nx,ny,nz),\
    txy_t2.zeros(nx,ny,nz),txz_t2.zeros(nx,ny,nz),tyz_t2.zeros(nx,ny,nz),\
    txy_pdx.zeros(nx,ny,nz),txy_pdy.zeros(nx,ny,nz),txz_pdx.zeros(nx,ny,nz),\
    txz_pdz.zeros(nx,ny,nz),tyz_pdy.zeros(nx,ny,nz),tyz_pdz.zeros(nx,ny,nz);
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
            mpar_lmd_add_2miu(i,j,k)=mpar_ro(i,j,k)*mpar_vp(i,j,k)*mpar_vp(i,j,k);
            mpar_1_dec_ro(i,j,k)=1.0/mpar_ro(i,j,k);
            //mo1[i][j]=(lmd[i][j]+2*miu[i][j])/(2*lmd[i][j]+2*miu[i][j])/ro[i][j];
    }}}
}

void calx_3d(fcube *ut2 ,fcube *ut1,fcube *u, fcube *m,int pjc,\
    class elastic3D_ARMA *obj)
{
    obj->calx_p(ut2[0],ut1[0],u[0],m[0],pjc);
}
void elastic3D_ARMA::calx_p(fcube &ut2 , const fcube &ut1, \
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
    Z=this->nz,
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
        float DX,DY,DZ,DT,xshd,suface_PML;
        int X,Y;
        DX=this->dx,DY=this->dy,DZ=this->dz,DT=this->dt,xshd=this->PML_wide;
        X=this->nx,Y=this->ny,suface_PML=this->isFreeSurface;
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
            du=0.5*(m(i+jc,j,k)+m(i,j,k))*du/DX;
            du1=C_X*3000.0*DX*(xshd+8-i)*DX*(xshd+8-i);
            //du1=C_X*this->vp(i,j,k)*DX*(xshd+8-i)*DX*(xshd+8-i);
            ut2(i,j,k)+=((du+ut1(i,j,k)*(1.0/DT-du1/2.0))/(1.0/DT+du1/2.0));
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
            du=0.5*(m(i+jc,j,k)+m(i,j,k))*du/DX;
            du1=0;
            ut2(i,j,k)+=((du-du1)*DT+ut1(i,j,k));
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
            du=0.5*(m(i+jc,j,k)+m(i,j,k))*du/DX;
            du1=C_X*3000.0*DX*(i-X+xshd+8)*DX*(i-X+xshd+8);
            //du1=C_X*this->vp(i,j,k)*DX*(i-X+xshd+8)*DX*(i-X+xshd+8);
            ut2(i,j,k)+=((du+ut1(i,j,k)*(1.0/DT-du1/2.0))/(1.0/DT+du1/2.0));
        }
    }}
}

void caly_3d(fcube *ut2 ,fcube *ut1,fcube *u, fcube *m,int pjc,\
    class elastic3D_ARMA *obj)
{
    obj->caly_p(ut2[0],ut1[0],u[0],m[0],pjc);
}
void elastic3D_ARMA::caly_p(fcube &ut2 , const fcube &ut1, \
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
        float DX,DY,DZ,DT,xshd,suface_PML;
        int Y,Z;
        int j,k,jc(pjc);
        float du1(0),du2(0),du(0);
        DX=this->dx,DY=this->dy,DZ=this->dz,DT=this->dt,xshd=this->PML_wide;
        Y=this->ny,Z=this->nz,suface_PML=this->isFreeSurface;
        
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
            du=0.5*(m(i,j+jc,k)+m(i,j,k))*du/DY;
            du1=C_Y*3000.0*DY*(xshd+8-j)*DY*(xshd+8-j);
            //du1=C_X*this->vp(i,j,k)*DX*(xshd+8-i)*DX*(xshd+8-i);
            ut2(i,j,k)+=((du+ut1(i,j,k)*(1.0/DT-du1/2.0))/(1.0/DT+du1/2.0));
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
            du=0.5*(m(i,j+jc,k)+m(i,j,k))*du/DY;
            du1=0;
            ut2(i,j,k)+=((du-du1)*DT+ut1(i,j,k));
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
            du=0.5*(m(i,j+jc,k)+m(i,j,k))*du/DY;
            du1=C_Y*3000.0*DY*(j-Y+xshd+8)*DY*(j-Y+xshd+8);
            //du1=C_Y*this->vp(i,j,k)*DY*(j-Y+xshd+8)*DY*(j-Y+xshd+8);
            ut2(i,j,k)+=((du+ut1(i,j,k)*(1.0/DT-du1/2.0))/(1.0/DT+du1/2.0));
        }
    }}
}

void calz_3d(fcube *ut2 ,fcube *ut1,fcube *u, fcube *m,int pjc,\
    class elastic3D_ARMA *obj)
{
    obj->calz_p(ut2[0],ut1[0],u[0],m[0],pjc);
}
void elastic3D_ARMA::calz_p(fcube &ut2 , const fcube &ut1,\
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
        float DX,DY,DZ,DT,xshd,suface_PML;
        int Y,Z,nSurface;
        int j,k,jc(pjc);
        float du1(0),du2(0),du(0);
        DX=this->dx,DY=this->dy,DZ=this->dz,DT=this->dt,xshd=this->PML_wide;
        Y=this->ny,Z=this->nz,suface_PML=this->isFreeSurface;
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
            //take care of the PML_boundary!!!
            if((suface_PML)<0.1 && j<nSurface)
                du=0.5*(m(i,j,k)+m(i,j,k))*du/DZ;
            else
                du=0.5*(m(i,j,k+jc)+m(i,j,k))*du/DZ;
            du1=C_Z*3000.0*DZ*suface_PML*(xshd+8-k)*DZ\
                *suface_PML*(xshd+8-k);
            //du1=C_Y*this->vp[j][i]*DY*suface_PML*(xshd+8-j)*DY\
                *suface_PML*(xshd+8-j);
            ut2(i,j,k)+=((du+ut1(i,j,k)*(1.0/DT-du1/2.0))/(1.0/DT+du1/2.0));
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
            //take care of the PML_boundary!!!
            if((suface_PML)<0.1  && j<nSurface)
                du=0.5*(m(i,j,k)+m(i,j,k))*du/DZ;
            else
                du=0.5*(m(i,j+jc,k)+m(i,j,k))*du/DZ;
            du1=0;
            ut2(i,j,k)+=((du-du1)*DT+ut1(i,j,k));
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
            du=0.5*(m(i,j,k+jc)+m(i,j,k))*du/DZ;
            du1=C_Z*3000.0*DZ*(k-Z+xshd+8)*DZ*(k-Z+xshd+8);
            //du1=C_Z*this->vp(i,j,k)*DZ*(k-Z+xshd+8)*DX*(k-Z+xshd+8);
            ut2(i,j,k)+=((du+ut1(i,j,k)*(1.0/DT-du1/2.0))/(1.0/DT+du1/2.0));
        }
    }}
}

void TimeSliceCal_elastic3D_ARMA_MultiThread(class elastic3D_ARMA & obj)
{
    int i,j,k;
    float **swap=NULL;
    thread *pcal;
    pcal=new thread[9];

    obj.data3d1.fill(0.0);obj.data3d2.fill(0.0);
    obj.data3d3.fill(0.0);obj.data3d4.fill(0.0);
    obj.data3d5.fill(0.0);obj.data3d6.fill(0.0);
    obj.vx_t2.fill(0.0),obj.vy_t2.fill(0.0),obj.vz_t2.fill(0.0);
    pcal[0]=thread(calx_3d,&obj.vx_t2,&obj.vx_t1,&obj.txx_t1,&obj.mpar_1_dec_ro,1,&obj);
    pcal[1]=thread(caly_3d,&obj.data3d1,&obj.vx_t1,&obj.txy_t1,&obj.mpar_1_dec_ro,0,&obj);
    pcal[2]=thread(calz_3d,&obj.data3d2,&obj.vx_t1,&obj.txz_t1,&obj.mpar_1_dec_ro,0,&obj);
    pcal[3]=thread(calx_3d,&obj.vy_t2,&obj.vy_t1,&obj.txy_t1,&obj.mpar_1_dec_ro,0,&obj);
    pcal[4]=thread(caly_3d,&obj.data3d3,&obj.vy_t1,&obj.tyy_t1,&obj.mpar_1_dec_ro,1,&obj);
    pcal[5]=thread(calz_3d,&obj.data3d4,&obj.vy_t1,&obj.tyz_t1,&obj.mpar_1_dec_ro,0,&obj);
    pcal[6]=thread(calx_3d,&obj.vz_t2,&obj.vz_t1,&obj.txz_t1,&obj.mpar_1_dec_ro,0,&obj);
    pcal[7]=thread(caly_3d,&obj.data3d5,&obj.vz_t1,&obj.tyz_t1,&obj.mpar_1_dec_ro,0,&obj);
    pcal[8]=thread(calz_3d,&obj.data3d6,&obj.vz_t1,&obj.tzz_t1,&obj.mpar_1_dec_ro,1,&obj);
    for(k=0;k<9;k++){
        pcal[k].join();
    }
    obj.vx_t1=obj.vx_t2+obj.data3d1+obj.data3d2;
    obj.vy_t1=obj.vy_t2+obj.data3d3+obj.data3d4;
    obj.vz_t1=obj.vz_t2+obj.data3d5+obj.data3d6;

    obj.data3d1.fill(0.0);obj.data3d2.fill(0.0);
    obj.data3d3.fill(0.0);obj.data3d4.fill(0.0);
    obj.data3d5.fill(0.0);obj.data3d6.fill(0.0);
    obj.txx_t2.fill(0.0),obj.tyy_t2.fill(0.0),obj.tzz_t2.fill(0.0);
    pcal[0]=thread(calx_3d,&obj.txx_t2,&obj.txx_t1,&obj.vx_t1,&obj.mpar_lmd_add_2miu,0,&obj);
    pcal[1]=thread(caly_3d,&obj.data3d1,&obj.txx_t1,&obj.vy_t1,&obj.mpar_lmd,0,&obj);
    pcal[2]=thread(calz_3d,&obj.data3d2,&obj.txx_t1,&obj.vz_t1,&obj.mpar_lmd,0,&obj);
    pcal[3]=thread(caly_3d,&obj.tyy_t2,&obj.tyy_t1,&obj.vy_t1,&obj.mpar_lmd_add_2miu,0,&obj);
    pcal[4]=thread(calx_3d,&obj.data3d3,&obj.tyy_t1,&obj.vx_t1,&obj.mpar_lmd,0,&obj);
    pcal[5]=thread(calz_3d,&obj.data3d4,&obj.tyy_t1,&obj.vz_t1,&obj.mpar_lmd,0,&obj);
    pcal[6]=thread(calz_3d,&obj.tzz_t2,&obj.tzz_t1,&obj.vz_t1,&obj.mpar_lmd_add_2miu,0,&obj);
    pcal[7]=thread(calx_3d,&obj.data3d5,&obj.tzz_t1,&obj.vx_t1,&obj.mpar_lmd,0,&obj);
    pcal[8]=thread(caly_3d,&obj.data3d6,&obj.tzz_t1,&obj.vy_t1,&obj.mpar_lmd,0,&obj);
    for(k=0;k<9;k++){
        pcal[k].join();
    }
    obj.txx_t1=obj.txx_t2+obj.data3d1+obj.data3d2;
    obj.tyy_t1=obj.tyy_t2+obj.data3d3+obj.data3d4;
    obj.tzz_t1=obj.tzz_t2+obj.data3d5+obj.data3d6;

    obj.data3d1.fill(0.0);obj.data3d2.fill(0.0);obj.data3d3.fill(0.0);
    obj.txz_t2.fill(0.0),obj.txy_t2.fill(0.0),obj.tyz_t2.fill(0.0);
    pcal[0]=thread(caly_3d,&obj.txy_t2,&obj.txy_t1,&obj.vx_t1,&obj.mpar_miu,1,&obj);
    pcal[1]=thread(calx_3d,&obj.data3d1,&obj.txy_t1,&obj.vy_t1,&obj.mpar_miu,1,&obj);
    pcal[2]=thread(calz_3d,&obj.txz_t2,&obj.txz_t1,&obj.vx_t1,&obj.mpar_miu,1,&obj);
    pcal[3]=thread(calx_3d,&obj.data3d2,&obj.txz_t1,&obj.vz_t1,&obj.mpar_miu,1,&obj);
    pcal[4]=thread(calz_3d,&obj.tyz_t2,&obj.tyz_t1,&obj.vy_t1,&obj.mpar_miu,1,&obj);
    pcal[5]=thread(caly_3d,&obj.data3d3,&obj.tyz_t1,&obj.vz_t1,&obj.mpar_miu,1,&obj);
    for(k=0;k<6;k++){
        pcal[k].join();
    }
    obj.txy_t1=obj.txy_t2+obj.data3d1;
    obj.txz_t1=obj.txz_t2+obj.data3d2;
    obj.tyz_t1=obj.tyz_t2+obj.data3d3;
    delete [] pcal; 
}

#endif

