/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef ELASTIC2D_HPP
#define ELASTIC_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <thread>
#include <future>
using namespace std;

#include "../xjc.h"
/////////////////////////////////////////////////////////////////////////////////////////////
float** newfmatcs(int x1, int x2, float f)
{
    float **p,**p2;
    int j;
    p=new float*[x1];      
    for(j=0;j<x1;j++)  
        {  
        p[j]=new float[x2];
        }
    p2=p;
    p=NULL;
    matcopy(p2,f,x1,x2); 
    return p2;
}

struct elastic_wave_data
{
    float **vxx=NULL, **vxy=NULL, **vyx=NULL, **vyy=NULL; //PML boundary
    float **vxx2=NULL, **vxy2=NULL, **vyx2=NULL, **vyy2=NULL; //PML boundary
    float **vxx3=NULL, **vxy3=NULL, **vyx3=NULL, **vyy3=NULL; //PML boundary
    float **txx=NULL, **txy=NULL, **tyy=NULL; //PML boundary
    float **txxx=NULL, **txyx=NULL, **tyyx=NULL; //PML boundary
    float **txxy=NULL, **txyy=NULL, **tyyy=NULL; //PML boundary
    float **txxx2=NULL, **txyx2=NULL, **tyyx2=NULL; //PML boundary
    float **txxy2=NULL, **txyy2=NULL, **tyyy2=NULL; //PML boundary
    float **vpx13=NULL,**vpx12=NULL, **vpx1=NULL;
    float **vpx23=NULL,**vpx22=NULL, **vpx2=NULL;
    float **vpy13=NULL,**vpy12=NULL, **vpy1=NULL;
    float **vpy23=NULL,**vpy22=NULL, **vpy2=NULL;
    float **vsy=NULL, **vsx=NULL;
    float **vpy=NULL, **vpx=NULL;

};

void timeslicecal_T_thread(class elastic2D &par);
class elastic2D
{ 
private:
//����ƫ΢����ɢ���ϵ��
    //float xs2[5]={1.666667/1.463612\
    12,-0.238095/1.463612,0.039683/1.463612,\
    -0.004960/1.463612,0.000317/1.463612};
    float xs2[5]={1.66666639,-0.238095194,3.96825150E-02,-4.96031437E-03,3.17460130E-04};
//һ��ƫ΢����ɢ���ϵ��
    float xs1[10]={-7.93651154E-04,9.92063992E-03,-5.95238246E-02,0.238095269,-0.833333373,\
0.833333373,-0.238095269,5.95238246E-02,-9.92063992E-03,7.93651154E-04};
    void calx(float** ut2 ,float** ut1,float **u, float **m,int jc);
    void caly(float** ut2 ,float** ut1,float **u, float **m,int jc);
    void calx0(float** ut,float** u, float **m,int jc);
    void caly0(float** ut,float** u, float **m,int jc);
    void calx2pml(float** ut ,float **u, float **m,int jc);
    void caly2pml(float** ut ,float **u, float **m,int jc);
    void calx2(float** ut2 ,float** ut1 ,float **ut0,float **u, float **m,int jc);
    void caly2(float** ut2 ,float** ut1 ,float **ut0,float **u, float **m,int jc);

public:
    elastic_wave_data data;
    float **vp=NULL,**vs=NULL,**lmd=NULL,**uz=NULL,**ux=NULL,\
        **ro1=NULL,**miu=NULL,**ro=NULL,**mo=NULL,**mo1=NULL,\
        **up=NULL,\
        **Txy=NULL,**Tyy=NULL,**Txx=NULL,vp1,vs1,roo; //velocity model
    float dx,dy,dt,R,dt2,t3,t5,C_X,C_Y;
    int nx,ny,suface,PML_wide,Zsiteofseasuface;
    int *par0=NULL,*par1=NULL;
    
    elastic2D();
    elastic2D(int x, int y);
    ~elastic2D();

    //void timeslicecal(int restart=0);
    void timeslicecal_T();
    void timeslicecal_V();
    void timeslicecal_U();
    void cleardata();
    void updatepar();
    void calx_p(float** ut2 ,float** ut1,float **u, float **m,int *pjc);
    void caly_p(float** ut2 ,float** ut1,float **u, float **m,int *pjc);
};
///////////////////////////////////////////////////////////////////////////////////////////

elastic2D::elastic2D()
{
    nx=0;ny=0;
    dx=5.0;dy=5.0;dt=0.0005;
    dt2=dt;
    PML_wide=40;suface=1;R=25;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

elastic2D::elastic2D(int z, int x)
{
    nx=x;ny=z;t3=1;
    dx=5.0;dy=5.0;dt=0.0005;dt2=dt;
    PML_wide=40;suface=1;R=25;
    vp1=3000,vs1=2000,roo=2000;
    C_Y=(R)*3/2/(PML_wide)/(PML_wide)/(PML_wide)/dy/dy/dy;
    C_X=(R)*3/2/(PML_wide)/(PML_wide)/(PML_wide)/dx/dx/dx;
    par0=new int[1],par0[0]=0;
    par1=new int[1],par1[0]=1;
    
    data.txx=newfmatcs(ny,nx,0.0),data.txxx2=newfmatcs(ny,nx,0.0),data.txxx=newfmatcs(ny,nx,0.0),\
    data.txxy2=newfmatcs(ny,nx,0.0),data.txxy=newfmatcs(ny,nx,0.0),data.txy=newfmatcs(ny,nx,0.0),\
    data.txyx2=newfmatcs(ny,nx,0.0),data.txyx=newfmatcs(ny,nx,0.0),data.txyy2=newfmatcs(ny,nx,0.0),\
    data.txyy=newfmatcs(ny,nx,0.0),data.tyy=newfmatcs(ny,nx,0.0),data.tyyx2=newfmatcs(ny,nx,0.0),\
    data.tyyx=newfmatcs(ny,nx,0.0),data.tyyy2=newfmatcs(ny,nx,0.0),data.tyyy=newfmatcs(ny,nx,0.0),\
    data.vxx2=newfmatcs(ny,nx,0.0),data.vxx=newfmatcs(ny,nx,0.0),data.vxy2=newfmatcs(ny,nx,0.0),\
    data.vxy=newfmatcs(ny,nx,0.0),data.vyx2=newfmatcs(ny,nx,0.0),data.vyx=newfmatcs(ny,nx,0.0),\
    data.vxx3=newfmatcs(ny,nx,0.0),data.vxy3=newfmatcs(ny,nx,0.0),\
    data.vyx3=newfmatcs(ny,nx,0.0),data.vyy3=newfmatcs(ny,nx,0.0),\
    data.vyy2=newfmatcs(ny,nx,0.0),data.vyy=newfmatcs(ny,nx,0.0);
    data.vpx13=newfmatcs(ny,nx,0.0),data.vpx12=newfmatcs(ny,nx,0.0),data.vpx1=newfmatcs(ny,nx,0.0);
    data.vpx23=newfmatcs(ny,nx,0.0),data.vpx22=newfmatcs(ny,nx,0.0),data.vpx2=newfmatcs(ny,nx,0.0);
    data.vpy13=newfmatcs(ny,nx,0.0),data.vpy12=newfmatcs(ny,nx,0.0),data.vpy1=newfmatcs(ny,nx,0.0);
    data.vpy23=newfmatcs(ny,nx,0.0),data.vpy22=newfmatcs(ny,nx,0.0),data.vpy2=newfmatcs(ny,nx,0.0);
    data.vsy=newfmatcs(ny,nx,0.0),data.vsx=newfmatcs(ny,nx,0.0);
    data.vpy=newfmatcs(ny,nx,0.0),data.vpx=newfmatcs(ny,nx,0.0);
//pure-pressure modeling

    vs=newfmatcs(ny,nx,vs1),vp=newfmatcs(ny,nx,vp1),ro=newfmatcs(ny,nx,roo);
    miu=newfmatcs(ny,nx,vs1*vs1*roo),lmd=newfmatcs(ny,nx,(vp1*vp1*roo-2*vs1*vs1*roo));
    mo=newfmatcs(ny,nx,vp1*vp1*roo),ro1=newfmatcs(ny,nx,1.0/roo); 
    ux=newfmatcs(ny,nx,0.0),uz=newfmatcs(ny,nx,0.0),up=newfmatcs(ny,nx,0.0);
    mo1=newfmatcs(ny,nx,(lmd[0][0]+2*miu[0][0])/(2*lmd[0][0]+2*miu[0][0])/ro[0][0]); 
    Txx=newfmatcs(ny,nx,0.0),Tyy=newfmatcs(ny,nx,0.0),Txy=newfmatcs(ny,nx,0.0);
}

elastic2D::~elastic2D()
{
    matdelete(vs,ny),matdelete(vp,ny),matdelete(lmd,ny),\
    matdelete(miu,ny),matdelete(ro,ny),matdelete(ro1,ny),matdelete(mo,ny);
    vp=NULL,vs=NULL,lmd=NULL,miu=NULL,ro=NULL,ro1=NULL,mo=NULL;

    matdelete(data.txx,ny),matdelete(data.txxx2,ny),matdelete(data.txxx,ny),\
    matdelete(data.txxy2,ny),matdelete(data.txxy,ny),matdelete(data.txy,ny),\
    matdelete(data.txyx2,ny),matdelete(data.txyx,ny),matdelete(data.txyy2,ny),\
    matdelete(data.txyy,ny),matdelete(data.tyy,ny),matdelete(data.tyyx2,ny),\
    matdelete(data.tyyx,ny),matdelete(data.tyyy2,ny),matdelete(data.tyyy,ny),\
    matdelete(data.vxx2,ny),matdelete(data.vxx,ny),matdelete(data.vxy2,ny),\
    matdelete(data.vxy,ny),matdelete(data.vyx2,ny),matdelete(data.vyx,ny),\
    matdelete(data.vyy2,ny),matdelete(data.vyy,ny);
    matdelete(data.vxx3,ny),matdelete(data.vxy3,ny),\
    matdelete(data.vyx3,ny),matdelete(data.vyy3,ny);
    matdelete(ux,ny),matdelete(uz,ny);
    matdelete(data.vpx23,ny),matdelete(data.vpx22,ny),matdelete(data.vpx1,ny);
    matdelete(data.vpx13,ny),matdelete(data.vpx12,ny),matdelete(data.vpx2,ny);
    matdelete(data.vpy23,ny),matdelete(data.vpy22,ny),matdelete(data.vpy1,ny);
    matdelete(data.vpy13,ny),matdelete(data.vpy12,ny),matdelete(data.vpy2,ny);
    matdelete(data.vsy,ny),matdelete(data.vsx,ny);
    matdelete(data.vpy,ny),matdelete(data.vpx,ny);
    matdelete(mo1,ny); 
    matdelete(Txx,ny),matdelete(Tyy,ny),matdelete(Txy,ny);
    uz=NULL,ux=NULL,mo1=NULL;
    Txy=NULL,Txx=NULL,Tyy=NULL;
    data.vxx=NULL, data.vxy=NULL, data.vyx=NULL, data.vyy=NULL; //PML boundary
    data.vxx2=NULL, data.vxy2=NULL, data.vyx2=NULL, data.vyy2=NULL; //PML boundary
    data.txx=NULL, data.txy=NULL, data.tyy=NULL; //PML boundary
    data.txxx=NULL, data.txyx=NULL, data.tyyx=NULL; //PML boundary
    data.txxy=NULL, data.txyy=NULL, data.tyyy=NULL; //PML boundary
    data.txxx2=NULL, data.txyx2=NULL, data.tyyx2=NULL; //PML boundary
    data.txxy2=NULL, data.txyy2=NULL, data.tyyy2=NULL; //PML boundaryt

    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void elastic2D::cleardata()
{
    //this->setvelocity(0.0);
    matdelete(data.txx,ny),matdelete(data.txxx2,ny),matdelete(data.txxx,ny),\
    matdelete(data.txxy2,ny),matdelete(data.txxy,ny),matdelete(data.txy,ny),\
    matdelete(data.txyx2,ny),matdelete(data.txyx,ny),matdelete(data.txyy2,ny),\
    matdelete(data.txyy,ny),matdelete(data.tyy,ny),matdelete(data.tyyx2,ny),\
    matdelete(data.tyyx,ny),matdelete(data.tyyy2,ny),matdelete(data.tyyy,ny),\
    matdelete(data.vxx2,ny),matdelete(data.vxx,ny),matdelete(data.vxy2,ny),\
    matdelete(data.vxy,ny),matdelete(data.vyx2,ny),matdelete(data.vyx,ny),\
    matdelete(data.vyy2,ny),matdelete(data.vyy,ny);
    matdelete(ux,ny),matdelete(uz,ny);
    uz=NULL,ux=NULL;
    data.vxx=NULL, data.vxy=NULL, data.vyx=NULL, data.vyy=NULL; //PML boundary
    data.vxx2=NULL, data.vxy2=NULL, data.vyx2=NULL, data.vyy2=NULL; //PML boundary
    data.txx=NULL, data.txy=NULL, data.tyy=NULL; //PML boundary
    data.txxx=NULL, data.txyx=NULL, data.tyyx=NULL; //PML boundary
    data.txxy=NULL, data.txyy=NULL, data.tyyy=NULL; //PML boundary
    data.txxx2=NULL, data.txyx2=NULL, data.tyyx2=NULL; //PML boundary
    data.txxy2=NULL, data.txyy2=NULL, data.tyyy2=NULL; //PML boundaryt
    cout<<"All Matrix data has been clear!"<<endl;

    data.txx=newfmatcs(ny,nx,0.0),data.txxx2=newfmatcs(ny,nx,0.0),data.txxx=newfmatcs(ny,nx,0.0),\
    data.txxy2=newfmatcs(ny,nx,0.0),data.txxy=newfmatcs(ny,nx,0.0),data.txy=newfmatcs(ny,nx,0.0),\
    data.txyx2=newfmatcs(ny,nx,0.0),data.txyx=newfmatcs(ny,nx,0.0),data.txyy2=newfmatcs(ny,nx,0.0),\
    data.txyy=newfmatcs(ny,nx,0.0),data.tyy=newfmatcs(ny,nx,0.0),data.tyyx2=newfmatcs(ny,nx,0.0),\
    data.tyyx=newfmatcs(ny,nx,0.0),data.tyyy2=newfmatcs(ny,nx,0.0),data.tyyy=newfmatcs(ny,nx,0.0),\
    data.vxx2=newfmatcs(ny,nx,0.0),data.vxx=newfmatcs(ny,nx,0.0),data.vxy2=newfmatcs(ny,nx,0.0),\
    data.vxy=newfmatcs(ny,nx,0.0),data.vyx2=newfmatcs(ny,nx,0.0),data.vyx=newfmatcs(ny,nx,0.0),\
    data.vyy2=newfmatcs(ny,nx,0.0),data.vyy=newfmatcs(ny,nx,0.0);
}
 
void elastic2D::updatepar()
{
    int i,j,n;
    dt2=dt;
    C_Y=(R)*3/2/(PML_wide)/(PML_wide)/(PML_wide)/dy/dy/dy;
    C_X=(R)*3/2/(PML_wide)/(PML_wide)/(PML_wide)/dx/dx/dx;
    t5=1.0;    
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            miu[i][j]=vs[i][j]*vs[i][j]*ro[i][j];
            lmd[i][j]=ro[i][j]*(vp[i][j]*vp[i][j]-2.0*vs[i][j]*vs[i][j]);
            mo[i][j]=ro[i][j]*vp[i][j]*vp[i][j];
            ro1[i][j]=1.0/ro[i][j];
            mo1[i][j]=(lmd[i][j]+2*miu[i][j])/(2*lmd[i][j]+2*miu[i][j])/ro[i][j];
        }
    }
}

void elastic2D::caly(float** ut2 ,float **ut1,float **u, float **m,int jc)
{
    int *pjc;
    pjc=new int[1];
    pjc[0]=jc;
    this->caly_p(ut2 ,ut1,u,m,pjc);
}
void caly_pp(float** ut2 ,float **ut1,float **u, float **m,int *pjc,\
    class elastic2D *par)
{
    par->caly_p(ut2 ,ut1,u,m,pjc);
}
void elastic2D::caly_p(float** ut2 ,float **ut1,float **u, float **m,int *pjc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1,jc;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    jc=pjc[0];
    
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
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
    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    {
    for(i=8;i<X-8;i++)
    {
        for(j=8;j<xshd+8;j++)
        {
            du=(
                C1*(u[j+0+jc][i]-u[j-1+jc][i])+\
                C2*(u[j+1+jc][i]-u[j-2+jc][i])+\
                C3*(u[j+2+jc][i]-u[j-3+jc][i])+\
                C4*(u[j+3+jc][i]-u[j-4+jc][i])+\
                C5*(u[j+4+jc][i]-u[j-5+jc][i])+\
                C6*(u[j+5+jc][i]-u[j-6+jc][i])+\
                C7*(u[j+6+jc][i]-u[j-7+jc][i])+\
                C8*(u[j+7+jc][i]-u[j-8+jc][i])\
                );
            //take care of the PML_boundary!!!
            if((suface_PML)==0 && j<this->Zsiteofseasuface)
                du=0.5*(m[j][i]+m[j][i])*du/DY;
            else
            du=0.5*(m[j+(jc)][i]+m[j][i])*du/DY;

            du1=C_Y*3000.0*DY*suface_PML*(xshd+8-j)*DY\
                *suface_PML*(xshd+8-j);
            //du1=C_Y*this->vp[j][i]*DY*suface_PML*(xshd+8-j)*DY\
                *suface_PML*(xshd+8-j);
            ut2[j][i]=(du+ut1[j][i]*(1/DT-du1/2))/(1/DT+du1/2);
        }
        for(j=xshd+8;j<Y-xshd-8;j++)
        {
            du=(
                C1*(u[j+0+jc][i]-u[j-1+jc][i])+\
                C2*(u[j+1+jc][i]-u[j-2+jc][i])+\
                C3*(u[j+2+jc][i]-u[j-3+jc][i])+\
                C4*(u[j+3+jc][i]-u[j-4+jc][i])+\
                C5*(u[j+4+jc][i]-u[j-5+jc][i])+\
                C6*(u[j+5+jc][i]-u[j-6+jc][i])+\
                C7*(u[j+6+jc][i]-u[j-7+jc][i])+\
                C8*(u[j+7+jc][i]-u[j-8+jc][i])\
                );
            //take care of the PML_boundary!!!
            if((suface_PML)==0  && j<this->Zsiteofseasuface)
                du=0.5*(m[j][i]+m[j][i])*du/DY;
            else
            du=0.5*(m[j+(jc)][i]+m[j][i])*du/DY;

            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=Y-xshd-8;j<Y-8;j++)
        {
            du=(
                C1*(u[j+0+jc][i]-u[j-1+jc][i])+\
                C2*(u[j+1+jc][i]-u[j-2+jc][i])+\
                C3*(u[j+2+jc][i]-u[j-3+jc][i])+\
                C4*(u[j+3+jc][i]-u[j-4+jc][i])+\
                C5*(u[j+4+jc][i]-u[j-5+jc][i])+\
                C6*(u[j+5+jc][i]-u[j-6+jc][i])+\
                C7*(u[j+6+jc][i]-u[j-7+jc][i])+\
                C8*(u[j+7+jc][i]-u[j-8+jc][i])\
                );
            du=0.5*(m[j+jc][i]+m[j][i])*du/DY;
            du1=C_Y*3000.0*DY*(j-Y+xshd+8)*DY*(j-Y+xshd+8);
            //du1=C_Y*this->vp[j][i]*DY*(j-Y+xshd+8)*DY*(j-Y+xshd+8);
            ut2[j][i]=(du+ut1[j][i]*(1/DT-du1/2))/(1/DT+du1/2);
        }
    }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::calx(float** ut2 ,float **ut1,float **u, float **m,int jc)
{
    int *pjc;
    pjc=new int[1];
    pjc[0]=jc;
    this->calx_p(ut2 ,ut1,u,m,pjc);
}
void calx_pp(float** ut2 ,float **ut1,float **u, float **m,int *pjc,\
    class elastic2D *par)
{
    par->calx_p(ut2 ,ut1,u,m,pjc);
}
void elastic2D::calx_p(float** ut2 ,float **ut1,float **u, float **m,int *pjc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1,jc(pjc[0]);
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
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
    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    {
    for(j=8;j<Y-8;j++)
    {
        for(i=8;i<xshd+8;i++)
        {  
            du=(
                C1*(u[j][i+0+jc]-u[j][i-1+jc])+\
                C2*(u[j][i+1+jc]-u[j][i-2+jc])+\
                C3*(u[j][i+2+jc]-u[j][i-3+jc])+\
                C4*(u[j][i+3+jc]-u[j][i-4+jc])+\
                C5*(u[j][i+4+jc]-u[j][i-5+jc])+\
                C6*(u[j][i+5+jc]-u[j][i-6+jc])+\
                C7*(u[j][i+6+jc]-u[j][i-7+jc])+\
                C8*(u[j][i+7+jc]-u[j][i-8+jc])\
                );
            du=0.5*(m[j][i+jc]+m[j][i])*du/DX;
            du1=C_X*3000.0*DX*(xshd+8-i)*DX*(xshd+8-i);
            //du1=C_X*this->vp[j][i]*DX*(xshd+8-i)*DX*(xshd+8-i);
            ut2[j][i]=(du+ut1[j][i]*(1/DT-du1/2))/(1/DT+du1/2);
        }
        for(i=xshd+8;i<X-xshd-8;i++)
        {  
            du=(
                C1*(u[j][i+0+jc]-u[j][i-1+jc])+\
                C2*(u[j][i+1+jc]-u[j][i-2+jc])+\
                C3*(u[j][i+2+jc]-u[j][i-3+jc])+\
                C4*(u[j][i+3+jc]-u[j][i-4+jc])+\
                C5*(u[j][i+4+jc]-u[j][i-5+jc])+\
                C6*(u[j][i+5+jc]-u[j][i-6+jc])+\
                C7*(u[j][i+6+jc]-u[j][i-7+jc])+\
                C8*(u[j][i+7+jc]-u[j][i-8+jc])\
                );
            du=0.5*(m[j][i+jc]+m[j][i])*du/DX;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=X-xshd-8;i<X-8;i++)
        {  
            du=(
                C1*(u[j][i+0+jc]-u[j][i-1+jc])+\
                C2*(u[j][i+1+jc]-u[j][i-2+jc])+\
                C3*(u[j][i+2+jc]-u[j][i-3+jc])+\
                C4*(u[j][i+3+jc]-u[j][i-4+jc])+\
                C5*(u[j][i+4+jc]-u[j][i-5+jc])+\
                C6*(u[j][i+5+jc]-u[j][i-6+jc])+\
                C7*(u[j][i+6+jc]-u[j][i-7+jc])+\
                C8*(u[j][i+7+jc]-u[j][i-8+jc])\
                );
            du=0.5*(m[j][i+jc]+m[j][i])*du/DX;
            du1=C_X*3000.0*DX*(i-X+xshd+8)*DX*(i-X+xshd+8);
            //du1=C_X*this->vp[j][i]*DX*(i-X+xshd+8)*DX*(i-X+xshd+8);
            ut2[j][i]=(du+ut1[j][i]*(1/DT-du1/2))/(1/DT+du1/2);
        }
    }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::caly2(float** ut2 ,float **ut1,float **ut0,float **u, float **m,int jc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    /*
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
    */
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    for(i=8;i<X-8;i++)
    {
        for(j=8;j<Y-8;j++)
        {
            du=(
                C1*(u[j+0+jc][i]-u[j-1+jc][i])+\
                C2*(u[j+1+jc][i]-u[j-2+jc][i])+\
                C3*(u[j+2+jc][i]-u[j-3+jc][i])+\
                C4*(u[j+3+jc][i]-u[j-4+jc][i])+\
                C5*(u[j+4+jc][i]-u[j-5+jc][i])+\
                C6*(u[j+5+jc][i]-u[j-6+jc][i])+\
                C7*(u[j+6+jc][i]-u[j-7+jc][i])+\
                C8*(u[j+7+jc][i]-u[j-8+jc][i])\
                );
            du=m[j][i]*du/DY;
            ut2[j][i]=DT*DT*(du)+2.0*ut1[j][i]-ut0[j][i];
        }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::calx2(float** ut2 ,float **ut1,float **ut0,float **u, float **m,int jc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    /*
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
    */
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    for(j=8;j<Y-8;j++)
    {
        for(i=8;i<X-8;i++)
        {  
            du=(
                C1*(u[j][i+0+jc]-u[j][i-1+jc])+\
                C2*(u[j][i+1+jc]-u[j][i-2+jc])+\
                C3*(u[j][i+2+jc]-u[j][i-3+jc])+\
                C4*(u[j][i+3+jc]-u[j][i-4+jc])+\
                C5*(u[j][i+4+jc]-u[j][i-5+jc])+\
                C6*(u[j][i+5+jc]-u[j][i-6+jc])+\
                C7*(u[j][i+6+jc]-u[j][i-7+jc])+\
                C8*(u[j][i+7+jc]-u[j][i-8+jc])\
                );
            du=m[j][i]*du/DX;
            ut2[j][i]=DT*DT*(du)+2.0*ut1[j][i]-ut0[j][i];
        }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::caly0(float** ut ,float **u, float **m,int jc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    /*
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
    */
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    for(i=8;i<X-8;i++)
    {
        for(j=8;j<Y-8;j++)
        {
            du=(
                C1*(u[j+0+jc][i]-u[j-1+jc][i])+\
                C2*(u[j+1+jc][i]-u[j-2+jc][i])+\
                C3*(u[j+2+jc][i]-u[j-3+jc][i])+\
                C4*(u[j+3+jc][i]-u[j-4+jc][i])+\
                C5*(u[j+4+jc][i]-u[j-5+jc][i])+\
                C6*(u[j+5+jc][i]-u[j-6+jc][i])+\
                C7*(u[j+6+jc][i]-u[j-7+jc][i])+\
                C8*(u[j+7+jc][i]-u[j-8+jc][i])\
                );
            du=m[j][i]*du/DY;
            ut[j][i]=du;
        }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::calx0(float** ut ,float **u, float **m,int jc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    /*
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
    */
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    for(j=8;j<Y-8;j++)
    {
        for(i=8;i<X-8;i++)
        {  
            du=(
                C1*(u[j][i+0+jc]-u[j][i-1+jc])+\
                C2*(u[j][i+1+jc]-u[j][i-2+jc])+\
                C3*(u[j][i+2+jc]-u[j][i-3+jc])+\
                C4*(u[j][i+3+jc]-u[j][i-4+jc])+\
                C5*(u[j][i+4+jc]-u[j][i-5+jc])+\
                C6*(u[j][i+5+jc]-u[j][i-6+jc])+\
                C7*(u[j][i+6+jc]-u[j][i-7+jc])+\
                C8*(u[j][i+7+jc]-u[j][i-8+jc])\
                );
            du=m[j][i]*du/DX;
            ut[j][i]=du;
        }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::caly2pml(float** ut ,float **u, float **m,int jc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    /*
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
    */
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    for(i=8;i<X-8;i++)
    {
        for(j=8;j<Y-8;j++)
        {
            du=(
                C1*(u[j+0+jc][i]-u[j-1+jc][i])+\
                C2*(u[j+1+jc][i]-u[j-2+jc][i])+\
                C3*(u[j+2+jc][i]-u[j-3+jc][i])+\
                C4*(u[j+3+jc][i]-u[j-4+jc][i])+\
                C5*(u[j+4+jc][i]-u[j-5+jc][i])+\
                C6*(u[j+5+jc][i]-u[j-6+jc][i])+\
                C7*(u[j+6+jc][i]-u[j-7+jc][i])+\
                C8*(u[j+7+jc][i]-u[j-8+jc][i])\
                );
            du=m[j][i]*du/DY;
            ut[j][i]=du;
        }
    }
    du1=0,du2=0,du=0;
}

void elastic2D::calx2pml(float** ut ,float **u, float **m,int jc)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    int i,j,i1,j1;
    float du1(0),du2(0),du(0);
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 
    /*
    float C1      = 1.2340911;
    float C2      = -1.0664985e-01;
    float C3      = 2.3036367e-02;
    float C4      = -5.3423856e-03;
    float C5      = 1.0772712e-03;
    float C6      = -1.6641888e-04;
    float C7      = 1.7021711e-005;
    float C8      = -8.5234642e-007;//差分系数
    */
    float C1=1.21124268;
    float C2=-8.97216797E-02;
    float C3=1.38427736E-02;
    float C4=-1.76565989E-03;   
    float C5=1.18679469E-04;
    float C6      = 0;
    float C7      = 0;
    float C8      = 0;//差分系数

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    for(j=8;j<Y-8;j++)
    {
        for(i=8;i<X-8;i++)
        {  
            du=(
                C1*(u[j][i+0+jc]-u[j][i-1+jc])+\
                C2*(u[j][i+1+jc]-u[j][i-2+jc])+\
                C3*(u[j][i+2+jc]-u[j][i-3+jc])+\
                C4*(u[j][i+3+jc]-u[j][i-4+jc])+\
                C5*(u[j][i+4+jc]-u[j][i-5+jc])+\
                C6*(u[j][i+5+jc]-u[j][i-6+jc])+\
                C7*(u[j][i+6+jc]-u[j][i-7+jc])+\
                C8*(u[j][i+7+jc]-u[j][i-8+jc])\
                );
            du=m[j][i]*du/DX;
            ut[j][i]=du;
        }
    }
    du1=0,du2=0,du=0;
}


void timeslicecal_T_thread(class elastic2D & par)
{
    int i,j,k;
    float **swap=NULL;
    thread *pcal;
    pcal=new thread[6];
    
    pcal[0]=thread(caly_pp,par.data.vyy2,par.data.vyy,par.Tyy,par.ro1,par.par0,&par);
    pcal[1]=thread(calx_pp,par.data.vyx2,par.data.vyx,par.Txy,par.ro1,par.par0,&par);
    pcal[2]=thread(caly_pp,par.data.vxy2,par.data.vxy,par.Txy,par.ro1,par.par1,&par);
    pcal[3]=thread(calx_pp,par.data.vxx2,par.data.vxx,par.Txx,par.ro1,par.par1,&par);
    for(k=0;k<4;k++){
        pcal[k].join();
    }
    //par.caly(par.data.vyy2,par.data.vyy,par.Tyy,par.ro1,0);
    swap=par.data.vyy,par.data.vyy=par.data.vyy2,par.data.vyy2=swap;
    //par.calx(par.data.vyx2,par.data.vyx,par.Txy,par.ro1,0);
    swap=par.data.vyx,par.data.vyx=par.data.vyx2,par.data.vyx2=swap;
    //par.caly(par.data.vxy2,par.data.vxy,par.Txy,par.ro1,1);
    swap=par.data.vxy,par.data.vxy=par.data.vxy2,par.data.vxy2=swap;
    //par.calx(par.data.vxx2,par.data.vxx,par.Txx,par.ro1,1);
    swap=par.data.vxx,par.data.vxx=par.data.vxx2,par.data.vxx2=swap;

    for(i=0;i<par.ny;i++){
        for(j=0;j<par.nx;j++){
            par.uz[i][j]=par.data.vyx[i][j]+par.data.vyy[i][j];
            par.ux[i][j]=par.data.vxx[i][j]+par.data.vxy[i][j];
        }}

    pcal[0]=thread(caly_pp,par.data.txyy2,par.data.txyy,par.ux,par.miu,par.par0,&par);
    pcal[1]=thread(calx_pp,par.data.txyx2,par.data.txyx,par.uz,par.miu,par.par1,&par);
    pcal[2]=thread(caly_pp,par.data.tyyy2,par.data.tyyy,par.uz,par.mo,par.par1,&par);
    pcal[3]=thread(calx_pp,par.data.tyyx2,par.data.tyyx,par.ux,par.lmd,par.par0,&par);
    pcal[4]=thread(caly_pp,par.data.txxy2,par.data.txxy,par.uz,par.lmd,par.par1,&par);
    pcal[5]=thread(calx_pp,par.data.txxx2,par.data.txxx,par.ux,par.mo,par.par0,&par);
    for(k=0;k<6;k++){
        pcal[k].join();
    }
    //par.caly(par.data.txyy2,par.data.txyy,par.ux,par.miu,0);
    swap=par.data.txyy,par.data.txyy=par.data.txyy2,par.data.txyy2=swap;
    //par.calx(par.data.txyx2,par.data.txyx,par.uz,par.miu,1);
    swap=par.data.txyx,par.data.txyx=par.data.txyx2,par.data.txyx2=swap;
    //par.caly(par.data.tyyy2,par.data.tyyy,par.uz,par.mo,1);
    swap=par.data.tyyy,par.data.tyyy=par.data.tyyy2,par.data.tyyy2=swap;
    //par.calx(par.data.tyyx2,par.data.tyyx,par.ux,par.lmd,0);
    swap=par.data.tyyx,par.data.tyyx=par.data.tyyx2,par.data.tyyx2=swap;
    //par.caly(par.data.txxy2,par.data.txxy,par.uz,par.lmd,1);
    swap=par.data.txxy,par.data.txxy=par.data.txxy2,par.data.txxy2=swap;
    //par.calx(par.data.txxx2,par.data.txxx,par.ux,par.mo,0);
    swap=par.data.txxx,par.data.txxx=par.data.txxx2,par.data.txxx2=swap;

    for(i=0;i<par.ny;i++){
        for(j=0;j<par.nx;j++){
            par.Txx[i][j]=par.data.txxx[i][j]+par.data.txxy[i][j];
            par.Tyy[i][j]=par.data.tyyx[i][j]+par.data.tyyy[i][j];
            par.Txy[i][j]=par.data.txyx[i][j]+par.data.txyy[i][j];
        }}
    swap=NULL;
}

void elastic2D::timeslicecal_T()
{
    int i,j;
    float **swap=NULL;
    this->caly(this->data.vyy2,this->data.vyy,this->Tyy,this->ro1,0);
    swap=this->data.vyy,this->data.vyy=this->data.vyy2,this->data.vyy2=swap;
    this->calx(this->data.vyx2,this->data.vyx,this->Txy,this->ro1,0);
    swap=this->data.vyx,this->data.vyx=this->data.vyx2,this->data.vyx2=swap;

    this->caly(this->data.vxy2,this->data.vxy,this->Txy,this->ro1,1);
    swap=this->data.vxy,this->data.vxy=this->data.vxy2,this->data.vxy2=swap;
    this->calx(this->data.vxx2,this->data.vxx,this->Txx,this->ro1,1);
    swap=this->data.vxx,this->data.vxx=this->data.vxx2,this->data.vxx2=swap;
/*
    this->calx(this->data.vpx12,this->data.vpx1,this->Txx,this->mo1,0);
    swap=this->data.vpx12,this->data.vpx12=this->data.vpx1,this->data.vpx1=swap;
    this->calx(this->data.vpx22,this->data.vpx2,this->Tyy,this->mo1,0);
    swap=this->data.vpx22,this->data.vpx22=this->data.vpx2,this->data.vpx2=swap;
    this->caly(this->data.vpy12,this->data.vpy1,this->Txx,this->mo1,0);
    swap=this->data.vpy12,this->data.vpy12=this->data.vpy1,this->data.vpy1=swap;
    this->caly(this->data.vpy22,this->data.vpy2,this->Tyy,this->mo1,0);
    swap=this->data.vpy22,this->data.vpy22=this->data.vpy2,this->data.vpy2=swap;
*/
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->uz[i][j]=this->data.vyx[i][j]+this->data.vyy[i][j];
            this->ux[i][j]=this->data.vxx[i][j]+this->data.vxy[i][j];/*
            this->data.vpx[i][j]=(this->data.vpx1[i][j]+this->data.vpx2[i][j]);
            this->data.vpy[i][j]=(this->data.vpy1[i][j]+this->data.vpy2[i][j]);
            this->data.vsx[i][j]=this->ux[i][j]-this->data.vpx[i][j];
            this->data.vsy[i][j]=this->uz[i][j]-this->data.vpy[i][j];*/
        }
    }
    this->caly(this->data.txyy2,this->data.txyy,this->ux,this->miu,0);
    swap=this->data.txyy,this->data.txyy=this->data.txyy2,this->data.txyy2=swap;
    this->calx(this->data.txyx2,this->data.txyx,this->uz,this->miu,1);
    swap=this->data.txyx,this->data.txyx=this->data.txyx2,this->data.txyx2=swap;

    this->caly(this->data.tyyy2,this->data.tyyy,this->uz,this->mo,1);
    swap=this->data.tyyy,this->data.tyyy=this->data.tyyy2,this->data.tyyy2=swap;
    this->calx(this->data.tyyx2,this->data.tyyx,this->ux,this->lmd,0);
    swap=this->data.tyyx,this->data.tyyx=this->data.tyyx2,this->data.tyyx2=swap;

    this->caly(this->data.txxy2,this->data.txxy,this->uz,this->lmd,1);
    swap=this->data.txxy,this->data.txxy=this->data.txxy2,this->data.txxy2=swap;
    this->calx(this->data.txxx2,this->data.txxx,this->ux,this->mo,0);
    swap=this->data.txxx,this->data.txxx=this->data.txxx2,this->data.txxx2=swap;

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->Txx[i][j]=this->data.txxx[i][j]+this->data.txxy[i][j];
            this->Tyy[i][j]=this->data.tyyx[i][j]+this->data.tyyy[i][j];
            this->Txy[i][j]=this->data.txyx[i][j]+this->data.txyy[i][j];
        }
    }
    swap=NULL;
}
void elastic2D::timeslicecal_V()
{this->timeslicecal_T();}
void elastic2D::timeslicecal_U()
{
    int i,j;
    float **swap=NULL;
    this->caly2(this->data.vyy3,this->data.vyy2,this->data.vyy,this->Tyy,this->ro1,0);
    swap=this->data.vyy,this->data.vyy=this->data.vyy2,\
    this->data.vyy2=this->data.vyy3,this->data.vyy3=swap;
    this->calx2(this->data.vyx3,this->data.vyx2,this->data.vyx,this->Txy,this->ro1,0);
    swap=this->data.vyx,this->data.vyx=this->data.vyx2,\
    this->data.vyx2=this->data.vyx3,this->data.vyx3=swap;

    this->caly2(this->data.vxy3,this->data.vxy2,this->data.vxy,this->Txy,this->ro1,1);
    swap=this->data.vxy,this->data.vxy=this->data.vxy2,\
    this->data.vxy2=this->data.vxy3,this->data.vxy3=swap;
    this->calx2(this->data.vxx3,this->data.vxx2,this->data.vxx,this->Txx,this->ro1,1);
    swap=this->data.vxx,this->data.vxx=this->data.vxx2,\
    this->data.vxx2=this->data.vxx3,this->data.vxx3=swap;

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->uz[i][j]=this->data.vyx2[i][j]+this->data.vyy2[i][j];
            this->ux[i][j]=this->data.vxx2[i][j]+this->data.vxy2[i][j];
        }
    }
    
    this->caly0(this->data.txyy,this->ux,this->miu,0);
    this->calx0(this->data.txyx,this->uz,this->miu,1);

    this->caly0(this->data.tyyy,this->uz,this->mo,1);
    this->calx0(this->data.tyyx,this->ux,this->lmd,0);

    this->caly0(this->data.txxy,this->uz,this->lmd,1);
    this->calx0(this->data.txxx,this->ux,this->mo,0);

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->Txx[i][j]=this->data.txxx[i][j]+this->data.txxy[i][j];
            this->Tyy[i][j]=this->data.tyyx[i][j]+this->data.tyyy[i][j];
            this->Txy[i][j]=this->data.txyx[i][j]+this->data.txyy[i][j];
        }
    }

//cal pure-pressure modeling
    this->calx2(this->data.vpx13,this->data.vpx12,this->data.vpx1,this->Txx,this->mo1,0);
    swap=this->data.vpx1,this->data.vpx1=this->data.vpx12,\
    this->data.vpx12=this->data.vpx13,this->data.vpx13=swap;
    this->calx2(this->data.vpx23,this->data.vpx22,this->data.vpx2,this->Tyy,this->mo1,0);
    swap=this->data.vpx2,this->data.vpx2=this->data.vpx22,\
    this->data.vpx22=this->data.vpx23,this->data.vpx23=swap;
    this->caly2(this->data.vpy13,this->data.vpy12,this->data.vpy1,this->Txx,this->mo1,1);
    swap=this->data.vpy1,this->data.vpy1=this->data.vpy12,\
    this->data.vpy12=this->data.vpy13,this->data.vpy13=swap;
    this->caly2(this->data.vpy23,this->data.vpy22,this->data.vpy2,this->Tyy,this->mo1,1);
    swap=this->data.vpy2,this->data.vpy2=this->data.vpy22,\
    this->data.vpy22=this->data.vpy23,this->data.vpy23=swap;
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            //this->up[i][j]=this->data.vpx22[i][j]+this->data.vpx12[i][j];
            this->up[i][j]=this->data.vpy22[i][j]+this->data.vpy12[i][j];
        }
    }

    swap=NULL;

}

float** gauss_souce_windows(int ny,int nx, int n1, int n2, int n)
{
    int i,j;
    float **w;
    w=newfmatcs(ny,nx,0.0);
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            //g=wavelet02(nt,dt,f,1.0);        
            w[i][j]=float(exp(-((i-n1)*(i-n1)+(j-n2)*(j-n2))/2.0/n)); 
        } 
    }
    return w;
}

void elastic_test(int dmovie=1)
{
    int Z(350),X(350),T(1500);
    elastic2D A(Z,X);
    ofstream outf1,outf2,outf3,outf4;
    float **uu,**w;
    uu=newfmatcs(Z,X,0.0);
    //outf1.open("u1movie.bin");
    //outf2.open("u2movie.bin");
    outf3.open("u3movie.bin");
    //outf4.open("u4movie.bin");
    w=gauss_souce_windows(Z,X,Z/2,X/2,5);
    datawrite(w,Z,X,"sourcewin.bin");
    int k,i,j;
    for(k=0;k<T;k++)
    {
        //addsouce(A.ux,Z/2,X/2,5,30,A.dt,k,1);
        for(i=0;i<A.ny;i++)
        {
            for(j=0;j<A.nx;j++)
            {
            A.ux[i][j]+=w[i][j]*wavelet02(k,A.dt,30,1);
            }
        }
        A.timeslicecal_T();

        //if(k%dmovie==0)
        {
            //matsmooth(uu,A.ux,Z,X,0);
            //datawrite(A.ux, Z, X, outf1);
            //matsmooth(uu,A.uz,Z,X,0);
            //datawrite(A.uz, Z, X, outf2);
            datawrite(A.data.vpx, Z, X, outf3);
            //datawrite(A.data.vsx, Z, X, outf4);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
   // outf1.close();
   // outf2.close();
    outf3.close();
   // outf4.close();
    cout<<"Have output test movie."<<endl;
}
////////////////////////////////////////////////////
void multiple_code(fmat& u2, fmat& u1, fmat& tu2u1,\
 float df, int fn1, int fn2)
{
    int n1(u1.n_rows),n2(u1.n_cols);
    int i,j,k;
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
    
    for(k=fn1;k<fn2;k++){
        w=-2.0*pi*df*k;
        cout<<"k="<<k<<" | "<<"w="<<w<<endl;
        for(i=0;i<n2;i++){
            for(j=0;j<n2;j++){
                a.imag(w*tu2u1(i,j));
                codemat(i,j)=exp(a);
            }
        }
        u2cx.col(k)=codemat*u1cx.col(k);
    }

    for(k=0;k<n2;k++){
        //u1cx.row(k)=fft(u1.col(k).t(),n1);
        onecol.col(0)=u2cx.row(k).st();
        u2.col(k)=real(ifft(onecol.col(0),n1));
    }
}

#endif

