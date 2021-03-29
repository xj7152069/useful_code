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
    float **txx=NULL, **txy=NULL, **tyy=NULL; //PML boundary
    float **txxx=NULL, **txyx=NULL, **tyyx=NULL; //PML boundary
    float **txxy=NULL, **txyy=NULL, **tyyy=NULL; //PML boundary
    float **txxx2=NULL, **txyx2=NULL, **tyyx2=NULL; //PML boundary
    float **txxy2=NULL, **txyy2=NULL, **tyyy2=NULL; //PML boundary
    float **vpx12=NULL, **vpx1=NULL;
    float **vpx22=NULL, **vpx2=NULL;
    float **vpy12=NULL, **vpy1=NULL;
    float **vpy22=NULL, **vpy2=NULL;
    float **vsy=NULL, **vsx=NULL;
    float **vpy=NULL, **vpx=NULL;
};

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
    void cal(float** ut2 ,float **ut1,float **u, float **m,const char & label);
    void calx(float** ut2 ,float **ut1,float **u, float **m,const char & label);
    void caly(float** ut2 ,float **ut1,float **u, float **m,const char & label);

public:
    elastic_wave_data data;
    float **vp=NULL,**vs=NULL,**lmd=NULL,**uz=NULL,**ux=NULL,\
        **ro1=NULL,**miu=NULL,**ro=NULL,**mo=NULL,**mo1=NULL,\
        **Txy=NULL,**Tyy=NULL,**Txx=NULL,vp1,vs1,roo; //velocity model
    float dx,dy,dt,R,dt2,t3;
    int nx,ny,suface,PML_wide;
    
    elastic2D();
    elastic2D(int x, int y);
    ~elastic2D();

    //void timeslicecal(int restart=0);
    void timeslicecal_u();
    void timeslicecal_T();
    void cleardata();
    void updatepar();
};
///////////////////////////////////////////////////////////////////////////////////////////

elastic2D::elastic2D()
{
    nx=0;ny=0;
    dx=5.0;dy=5.0;dt=0.0005;
    dt2=dt;
    PML_wide=30;suface=1;R=20;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

elastic2D::elastic2D(int z, int x)
{
    nx=x;ny=z;t3=1;
    dx=5.0;dy=5.0;dt=0.0005;dt2=dt;
    PML_wide=30;suface=1;R=20;
    vp1=3000,vs1=2000,roo=2000;
    
    data.txx=newfmatcs(ny,nx,0.0),data.txxx2=newfmatcs(ny,nx,0.0),data.txxx=newfmatcs(ny,nx,0.0),\
    data.txxy2=newfmatcs(ny,nx,0.0),data.txxy=newfmatcs(ny,nx,0.0),data.txy=newfmatcs(ny,nx,0.0),\
    data.txyx2=newfmatcs(ny,nx,0.0),data.txyx=newfmatcs(ny,nx,0.0),data.txyy2=newfmatcs(ny,nx,0.0),\
    data.txyy=newfmatcs(ny,nx,0.0),data.tyy=newfmatcs(ny,nx,0.0),data.tyyx2=newfmatcs(ny,nx,0.0),\
    data.tyyx=newfmatcs(ny,nx,0.0),data.tyyy2=newfmatcs(ny,nx,0.0),data.tyyy=newfmatcs(ny,nx,0.0),\
    data.vxx2=newfmatcs(ny,nx,0.0),data.vxx=newfmatcs(ny,nx,0.0),data.vxy2=newfmatcs(ny,nx,0.0),\
    data.vxy=newfmatcs(ny,nx,0.0),data.vyx2=newfmatcs(ny,nx,0.0),data.vyx=newfmatcs(ny,nx,0.0),\
    data.vyy2=newfmatcs(ny,nx,0.0),data.vyy=newfmatcs(ny,nx,0.0);
    data.vpx12=newfmatcs(ny,nx,0.0),data.vpx1=newfmatcs(ny,nx,0.0);
    data.vpx22=newfmatcs(ny,nx,0.0),data.vpx2=newfmatcs(ny,nx,0.0);
    data.vpy12=newfmatcs(ny,nx,0.0),data.vpy1=newfmatcs(ny,nx,0.0);
    data.vpy22=newfmatcs(ny,nx,0.0),data.vpy2=newfmatcs(ny,nx,0.0);
    data.vsy=newfmatcs(ny,nx,0.0),data.vsx=newfmatcs(ny,nx,0.0);
    data.vpy=newfmatcs(ny,nx,0.0),data.vpx=newfmatcs(ny,nx,0.0);

    vs=newfmatcs(ny,nx,vs1),vp=newfmatcs(ny,nx,vp1),ro=newfmatcs(ny,nx,roo);
    miu=newfmatcs(ny,nx,vs1*vs1*roo),lmd=newfmatcs(ny,nx,(vp1*vp1*roo-2*vs1*vs1*roo));
    mo=newfmatcs(ny,nx,vp1*vp1*roo),ro1=newfmatcs(ny,nx,1.0/roo); 
    ux=newfmatcs(ny,nx,0.0),uz=newfmatcs(ny,nx,0.0);
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
    matdelete(ux,ny),matdelete(uz,ny);
    matdelete(data.vpx22,ny),matdelete(data.vpx1,ny);
    matdelete(data.vpx12,ny),matdelete(data.vpx2,ny);
    matdelete(data.vpy22,ny),matdelete(data.vpy1,ny);
    matdelete(data.vpy12,ny),matdelete(data.vpy2,ny);
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
    int i,j;
    dt2=dt;
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

void elastic2D::caly(float** ut2 ,float **ut1,float **u, float **m,const char & label)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    float fdx,fdy,fddx,fddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,i1,j1,n,t4;
    float du1(0),du2(0),du(0),dux(0),duy(0),duxy(0);
    float DT2,DT3,DX2,DY2,mo2;
    float C_X, C_Y;
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    DT2=DT*DT,DT3=DT2*DT,DX2=DX*DX,DY2=DY*DY;
    t5=float(DY)/DX;
    C_Y=(this->R)*3/2/(xshd-5)/(xshd-5)/(xshd-5)/DY2/DY;
    C_X=(this->R)*3/2/(xshd-5)/(xshd-5)/(xshd-5)/DX2/DX;

    //else if(label=='y')
    {
    for(i=5;i<X-5;i++)
    {
        for(j=5;j<xshd;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=C_Y*this->vp[j][i]*DY*(xshd-j)*DY*(xshd-j)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=xshd;j<Y-xshd;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=Y-xshd;j<Y-5;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=C_Y*this->vp[j][i]*DY*(j-Y+xshd+1)*DY*(j-Y+xshd+1)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
    }
    }
    du1=0,du2=0,du=0,dux=0,duy=0,duxy=0;
    for(j=0;j<Y;j++)
    {
        for(i=0;i<X;i++)
        {
            ut1[j][i]=(ut2[j][i]);
        }
    }
}

void elastic2D::calx(float** ut2 ,float **ut1,float **u, float **m,const char & label)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    float fdx,fdy,fddx,fddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,i1,j1,n,t4;
    float du1(0),du2(0),du(0),dux(0),duy(0),duxy(0);
    float DT2,DT3,DX2,DY2,mo2;
    float C_X, C_Y;
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    DT2=DT*DT,DT3=DT2*DT,DX2=DX*DX,DY2=DY*DY;
    t5=float(DY)/DX;
    C_Y=(this->R)*3/2/(xshd-5)/(xshd-5)/(xshd-5)/DY2/DY;
    C_X=(this->R)*3/2/(xshd-5)/(xshd-5)/(xshd-5)/DX2/DX;

    //if(label=='x')
    {
    for(j=5;j<Y-5;j++)
    {
        for(i=5;i<xshd;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=C_X*this->vp[j][i]*DX*(xshd-i)*DX*(xshd-i)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=xshd;i<X-xshd;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=X-xshd;i<X-5;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=C_X*this->vp[j][i]*DX*(i-X+xshd+1)*DX*(i-X+xshd+1)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
    }
    }
    du1=0,du2=0,du=0,dux=0,duy=0,duxy=0;
    for(j=0;j<Y;j++)
    {
        for(i=0;i<X;i++)
        {
            ut1[j][i]=(ut2[j][i]);
        }
    }
}

void elastic2D::cal(float** ut2 ,float **ut1,float **u, float **m,const char & label)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    float fdx,fdy,fddx,fddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,i1,j1,n,t4;
    float du1(0),du2(0),du(0),dux(0),duy(0),duxy(0);
    float DT2,DT3,DX2,DY2,mo2;
    float C_X, C_Y;
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    DT2=DT*DT,DT3=DT2*DT,DX2=DX*DX,DY2=DY*DY;
    t5=float(DY)/DX;
    C_Y=this->R*3/2/(xshd-5)/(xshd-5)/(xshd-5)/DY2/DY;
    C_X=this->R*3/2/(xshd-5)/(xshd-5)/(xshd-5)/DX2/DX;

    if(label=='x')
    {
    for(j=5;j<Y-5;j++)
    {
        for(i=5;i<xshd;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=C_X*this->vp[j][i]*DX*(xshd-i)*DX*(xshd-i)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=xshd;i<X-xshd;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=X-xshd;i<X-5;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=C_X*this->vp[j][i]*DX*(i-X+xshd+1)*DX*(i-X+xshd+1)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
    }
    }
    else if(label=='y')
    {
    for(i=5;i<X-5;i++)
    {
        for(j=5;j<xshd;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=C_Y*this->vp[j][i]*DY*(xshd-j)*DY*(xshd-j)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=xshd;j<Y-xshd;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=Y-xshd;j<Y-5;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=C_Y*this->vp[j][i]*DY*(j-Y+xshd+1)*DY*(j-Y+xshd+1)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
    }
    }
    else
    {
        cout<<"error!"<<endl;
    }
            du1=0,du2=0,du=0,dux=0,duy=0,duxy=0;
    for(j=0;j<Y;j++)
    {
        for(i=0;i<X;i++)
        {
            ut1[j][i]=(ut2[j][i]);
        }
    }
}

void elastic2D::timeslicecal_u()
{
    int i,j;
    //cal(float** ut2 ,float **ut1,float **u, float **m,const char & label)
    this->calx(this->data.txyy2,this->data.txyy,this->ux,this->miu,'y');
    this->caly(this->data.txyx2,this->data.txyx,this->uz,this->miu,'x');

    this->calx(this->data.tyyy2,this->data.tyyy,this->uz,this->mo,'y');
    this->caly(this->data.tyyx2,this->data.tyyx,this->ux,this->lmd,'x');

    this->calx(this->data.txxy2,this->data.txxy,this->uz,this->lmd,'y');
    this->caly(this->data.txxx2,this->data.txxx,this->ux,this->mo,'x');

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->Txx[i][j]=this->data.txxx[i][j]+this->data.txxy[i][j];
            this->Tyy[i][j]=this->data.tyyx[i][j]+this->data.tyyy[i][j];
            this->Txy[i][j]=this->data.txyx[i][j]+this->data.txyy[i][j];
        }
    }
    this->calx(this->data.vyy2,this->data.vyy,this->Tyy,this->ro1,'n');
    this->caly(this->data.vyx2,this->data.vyx,this->Txy,this->ro1,'n');

    this->calx(this->data.vxy2,this->data.vxy,this->Txy,this->ro1,'n');
    this->caly(this->data.vxx2,this->data.vxx,this->Txx,this->ro1,'n');

    this->caly(this->data.vpx12,this->data.vpx1,this->Txx,this->mo1,'n');
    this->caly(this->data.vpx22,this->data.vpx2,this->Tyy,this->mo1,'n');
    this->calx(this->data.vpy12,this->data.vpy1,this->Txx,this->mo1,'n');
    this->calx(this->data.vpy22,this->data.vpy2,this->Tyy,this->mo1,'n');

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->uz[i][j]=this->data.vyx[i][j]+this->data.vyy[i][j];
            this->ux[i][j]=this->data.vxx[i][j]+this->data.vxy[i][j];
            this->data.vpx[i][j]=(this->data.vpx1[i][j]+this->data.vpx2[i][j]);
            this->data.vpy[i][j]=(this->data.vpy1[i][j]+this->data.vpy2[i][j]);
            this->data.vsx[i][j]=this->ux[i][j]-this->data.vpx[i][j];
            this->data.vsy[i][j]=this->uz[i][j]-this->data.vpy[i][j];
        }
    }
    //this->t3=this->t3*(-1);
}

void elastic2D::timeslicecal_T()
{
    int i,j;
    this->caly(this->data.vyy2,this->data.vyy,this->Tyy,this->ro1,'y');
    this->calx(this->data.vyx2,this->data.vyx,this->Txy,this->ro1,'x');

    this->caly(this->data.vxy2,this->data.vxy,this->Txy,this->ro1,'y');
    this->calx(this->data.vxx2,this->data.vxx,this->Txx,this->ro1,'x');

    this->calx(this->data.vpx12,this->data.vpx1,this->Txx,this->mo1,'x');
    this->calx(this->data.vpx22,this->data.vpx2,this->Tyy,this->mo1,'x');
    this->caly(this->data.vpy12,this->data.vpy1,this->Txx,this->mo1,'y');
    this->caly(this->data.vpy22,this->data.vpy2,this->Tyy,this->mo1,'y');

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->uz[i][j]=this->data.vyx[i][j]+this->data.vyy[i][j];
            this->ux[i][j]=this->data.vxx[i][j]+this->data.vxy[i][j];
            this->data.vpx[i][j]=(this->data.vpx1[i][j]+this->data.vpx2[i][j]);
            this->data.vpy[i][j]=(this->data.vpy1[i][j]+this->data.vpy2[i][j]);
            this->data.vsx[i][j]=this->ux[i][j]-this->data.vpx[i][j];
            this->data.vsy[i][j]=this->uz[i][j]-this->data.vpy[i][j];
        }
    }
    //cal(float** ut2 ,float **ut1,float **u, float **m,const char & label)
    this->caly(this->data.txyy2,this->data.txyy,this->ux,this->miu,'y');
    this->calx(this->data.txyx2,this->data.txyx,this->uz,this->miu,'x');

    this->caly(this->data.tyyy2,this->data.tyyy,this->uz,this->mo,'y');
    this->calx(this->data.tyyx2,this->data.tyyx,this->ux,this->lmd,'x');

    this->caly(this->data.txxy2,this->data.txxy,this->uz,this->lmd,'y');
    this->calx(this->data.txxx2,this->data.txxx,this->ux,this->mo,'x');

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->Txx[i][j]=this->data.txxx[i][j]+this->data.txxy[i][j];
            this->Tyy[i][j]=this->data.tyyx[i][j]+this->data.tyyy[i][j];
            this->Txy[i][j]=this->data.txyx[i][j]+this->data.txyy[i][j];
        }
    }
    //this->t3=this->t3*(-1);
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
        A.timeslicecal_u();

        if(k%dmovie==0)
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

#endif

