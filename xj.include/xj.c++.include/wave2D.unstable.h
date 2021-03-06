/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef WAVE2D_UNSTABLE_H_H
#define WAVE2D_UNSTABLE_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "mat.h"
#include "wave2D.h"
/////////////////////////////////////////////////////////////////////////////////////////////
void wave2D_unstable_stabletest(int Z, int X, int T, float dz, float dx, float dt, float v, int sz, int sx, float sf, float pmlwide, int dmovie);
void wave2D_unstable_stabletest(int Z, int X, int T, float dz, float dx, float dt, float **v, int sz, int sx, float sf, float pmlwide, int dmovie);


class wave2D_unstable
{
private:
    float xs2[5]={1.666667,-0.238095,0.039683,-0.004960,0.000317};
    float xs1[10]={-0.0007926,0.00991800,-0.0595200,0.238080,-0.833333,\
0.833333,-0.238080,0.0595200,-0.00991800,0.0007926};
    float **sx11=NULL, **sx12=NULL, **sx13=NULL; //PML boundary
    float **sx21=NULL, **sx22=NULL, **sx23=NULL, **sx24=NULL; //PML boundary
    float **sx31=NULL,  **sx32=NULL, **sx33=NULL; //PML boundary

public:
    float **p2=NULL; //velocity model
    float **s1=NULL, **s2=NULL, **s3=NULL; //time slices, add source to "s2"
    float dx,dy,dt,PML_wide,R;
    int nx,ny,suface;
    
    wave2D_unstable();
    wave2D_unstable(int x, int y);
    ~wave2D_unstable();

    void setvelocity(float v=3000.0);
    void timeslicecal();
    void timeslicecopy();
    void cleardata();
};
///////////////////////////////////////////////////////////////////////////////////////////

wave2D_unstable::wave2D_unstable()
{
    nx=0;ny=0;
    dx=5.0;dy=5.0;dt=0.0005;
    PML_wide=30;suface=1;R=10000000;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

wave2D_unstable::wave2D_unstable(int z, int x)
{
    nx=x;ny=z;
    dx=5.0;dy=5.0;dt=0.0005;
    PML_wide=30;suface=1;R=10000000;
    int i,j;
    s1=new float*[ny];s2=new float*[ny];s3=new float*[ny];  
    sx11=new float*[ny];sx12=new float*[ny];sx13=new float*[ny];     
    sx21=new float*[ny];sx22=new float*[ny];sx23=new float*[ny];sx24=new float*[ny];     
    sx31=new float*[ny];sx32=new float*[ny];sx33=new float*[ny];

    for(j=0;j<ny;j++)  
        {  
        s1[j]=new float[nx];s2[j]=new float[nx];s3[j]=new float[nx];   
        sx11[j]=new float[nx];sx12[j]=new float[nx];sx13[j]=new float[nx];    
        sx21[j]=new float[nx];sx22[j]=new float[nx];sx23[j]=new float[nx];sx24[j]=new float[nx];   
        sx31[j]=new float[nx];sx32[j]=new float[nx];sx33[j]=new float[nx];   
        }

    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            s1[j][i]=0.0;s2[j][i]=0.0;s3[j][i]=0.0;
            sx11[j][i]=0.0;sx12[j][i]=0.0;sx13[j][i]=0.0;
            sx21[j][i]=0.0;sx22[j][i]=0.0;sx23[j][i]=0.0;sx24[j][i]=0.0;
            sx31[j][i]=0.0;sx32[j][i]=0.0;sx33[j][i]=0.0;
            }
        }

    p2=new float*[ny];      
    for(j=0;j<ny;j++)  
        {  
        p2[j]=new float[nx];
        } 

    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            p2[j][i]=3000;
            }
        }
}

wave2D_unstable::~wave2D_unstable()
{
    int i,j;

    for(int i=0;i<ny;i++)  
        delete []p2[i]; 
    delete []p2;

    for(int i=0;i<ny;i++)  
       {
        delete []s1[i];delete []s2[i];delete []s3[i];
        delete []sx11[i];delete []sx12[i];delete []sx13[i];
        delete []sx21[i];delete []sx22[i];delete []sx23[i];delete []sx24[i];
        delete []sx31[i];delete []sx32[i];delete []sx33[i];
       }
    delete []s1;delete []s2;delete []s3;
    delete []sx11;delete []sx12;delete []sx13;
    delete []sx21;delete []sx22;delete []sx23;delete []sx24;
    delete []sx31;delete []sx32;delete []sx33;
    p2=NULL;s1=NULL;s2=NULL;s3=NULL;
    sx11=NULL;sx12=NULL;sx13=NULL;
    sx21=NULL;sx22=NULL;sx23=NULL;sx24=NULL;
    sx31=NULL;sx32=NULL;sx33=NULL;
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void wave2D_unstable::setvelocity(float v)
{
    int i,j;
    for(i=0;i<ny;i++)
        {
        for(j=0;j<nx;j++)
            {
            p2[i][j]=v;
            }
        }
}

void wave2D_unstable::cleardata()
{
    //this->setvelocity(0.0);
    int i,j;
    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            s1[j][i]=0.0;s2[j][i]=0.0;s3[j][i]=0.0;
            sx11[j][i]=0.0;sx12[j][i]=0.0;sx13[j][i]=0.0;
            sx21[j][i]=0.0;sx22[j][i]=0.0;sx23[j][i]=0.0;sx24[j][i]=0.0;
            sx31[j][i]=0.0;sx32[j][i]=0.0;sx33[j][i]=0.0;
            }
        }
    cout<<"All Matrix data has been clear!"<<endl;
}

void wave2D_unstable::timeslicecal()
{
    static float DX=dx,DY=dy,DT=dt,xshd=PML_wide,*xs1_in=xs1,*xs2_in=xs2;
    static int X=nx,Y=ny,suface_PML=suface;

    static float dx,dy,ddx,ddy,snx1,sny1,snx2,sny2,t2,t5=float(Y)/X;
    static int i,j,n,t3,t4;
    static float u1(0),u2(0),u(0),ux(0),uy(0);
    static float DT2=DT*DT,DT3=DT2*DT,DX2=DX*DX,DY2=DY*DY,mo2;
    static float C_X=3/2/(xshd)/DX*log(R)/(xshd)/DX/(xshd)/DX;
    static float C_Y=3/2/(xshd)/DY*log(R)/(xshd)/DY/(xshd)/DY;

    static float **sx11_in=sx11, **sx12_in=sx12, **sx13_in=sx13; //PML boundary
    static float **sx21_in=sx21, **sx22_in=sx22, **sx23_in=sx23, **sx24_in=sx24; //PML boundary
    static float **sx31_in=sx31, **sx32_in=sx32, **sx33_in=sx33; //PML boundary
    static float **p2_in=p2, **swap; //velocity model and swap
    static float **s1_in=s1, **s2_in=s2, **s3_in=s3; //time slices, add source to "s2"

    for(i=5;i<X-5;i++)
        {

        if(i>=0.5*(X))
            {t3=1;}
        else
            {t3=-1;}

        for(j=5;j<Y-5;j++)
            {

            for(n=0;n<5;n++)
                {  
                u=u+2*xs2_in[n];
                u1=u1+xs2_in[n]*(s2_in[j-n-1][i]+s2_in[j+n+1][i]);
                u2=u2+xs2_in[n]*(s2_in[j][i-n-1]+s2_in[j][i+n+1]);
                }

            for(n=0;n<10;n++)
                {
                if(n<5)
                    {
                    ux=ux+s2_in[j][i+n-5]*xs1_in[n];
                    uy=uy+s2_in[j+n-5][i]*xs1_in[n];
                    }
                else
                    {
                    ux=ux+s2_in[j][i+n-4]*xs1_in[n];
                    uy=uy+s2_in[j+n-4][i]*xs1_in[n];
                    }
                }

            snx1=0.0;snx2=0.0;
            sny1=0.0;sny2=0.0;
            dx=0;ddx=0;dy=0;ddy=0;

            if(suface_PML==1)
                {
                if(i>=X-xshd-5 && j<=t5*i && j>=-t5*i+Y)	
                    snx1=i-(X-xshd-5);
                if(i<=xshd+5 && j>=t5*i && j<=-t5*i+Y)		
                    snx2=xshd+5-i;
                if(j>=Y-xshd-5 && j>=t5*i && j>=-t5*i+Y)		
                    sny1=j-(Y-xshd-5);
                if(j<=xshd+5 && j<=t5*i && j<=-t5*i+Y)		
                    sny2=xshd+5-j;		
                }  //处理角落
            else
                {
                if(i>=X-xshd-5 && j<=t5*i)	
                    snx1=i-(X-xshd-5);
                if(i<=xshd+5 && j<=-t5*i+Y)		
                    snx2=xshd+5-i;
                if(j>=Y-xshd-5 && j>t5*i && j>-t5*i+Y)		
                    sny1=j-(Y-xshd-5);
                if(j<=xshd+5)		
                    sny2=0;		
                }

            if(sny1 !=0)
                {
                dy=p2_in[j][i]*C_Y*sny1*sny1*DY2;
                ddy=p2_in[j][i]*C_Y*2*sny1*DY;	
                }
            if(sny2 !=0)
                {
                dy=p2_in[j][i]*C_Y*sny2*sny2*DY2;
                ddy=p2_in[j][i]*C_Y*2*sny2*DY;	
                }
            if(snx1 !=0)
                {
                dx=p2_in[j][i]*C_X*snx1*snx1*DX2;
                ddx=p2_in[j][i]*C_X*2*snx1*DX;	
                }
            if(snx2 !=0)
                {
                dx=p2_in[j][i]*C_X*snx2*snx2*DX2;
                ddx=p2_in[j][i]*C_X*2*snx2*DX;	
                }

            if(j>=0.5*(Y))
                {t4=1;}
            else
                {t4=-1;}

            mo2=p2_in[j][i]*p2_in[j][i];

            if(snx1!=0 || snx2!=0)
                {
            //根据系数求得一阶偏微分的离散算子
                for(n=0;n<10;n++)
                    {
                    if(n<5)
                        {
                        ux=ux+s2_in[j][i+n-5]*xs1_in[n];
                        uy=uy+s2_in[j+n-5][i]*xs1_in[n];
                        }
                    else
                        {
                        ux=ux+s2_in[j][i+n-4]*xs1_in[n];
                        uy=uy+s2_in[j+n-4][i]*xs1_in[n];
                        }
                    }

                sx11_in[j][i]=mo2\
                *DT2*(u2-u*s2_in[j][i])*(1.0/(DX2))\
                -dx*dx*DT2*sx12_in[j][i]+(2*sx12_in[j][i]\
                -sx13_in[j][i])+DT*(2*dx*(sx13_in[j][i]-sx12_in[j][i]));

                sx21_in[j][i]=(-mo2\
                *ddx*(1.0/(DX))*(ux*t3)-dx*dx*dx*sx22_in[j][i]\
                +3*dx*dx*(sx23_in[j][i]-sx22_in[j][i])/DT+3*dx*\
                (2*sx23_in[j][i]-sx22_in[j][i]-sx24_in[j][i])/(DT2)\
                +(3*sx22_in[j][i]-3*sx23_in[j][i]+sx24_in[j][i])\
                /(DT3))*(DT3);

                sx31_in[j][i]=DT2*mo2\
                *(1.0/(DY2))*(u1-u*s2_in[j][i])+2*sx32_in[j][i]\
                -sx33_in[j][i];
                }
            else if(sny1!=0 || sny2!=0)
                {
            //根据系数求得一阶偏微分的离散算子
                for(n=0;n<10;n++)
                    {
                    if(n<5)
                        {
                        ux=ux+s2_in[j][i+n-5]*xs1_in[n];
                        uy=uy+s2_in[j+n-5][i]*xs1_in[n];
                        }
                    else
                        {
                        ux=ux+s2_in[j][i+n-4]*xs1_in[n];
                        uy=uy+s2_in[j+n-4][i]*xs1_in[n];
                        }
                    }

                sx11_in[j][i]=mo2\
                *DT2*(u1-u*s2_in[j][i])*(1.0/(DY2))\
                -dy*dy*DT2*sx12_in[j][i]+(2*sx12_in[j][i]\
                -sx13_in[j][i])+DT*(2*dy*(sx13_in[j][i]-sx12_in[j][i]));

                sx21_in[j][i]=(-mo2\
                *ddy*(1.0/(DY))*(uy*t4)-dy*dy*dy*sx22_in[j][i]\
                +3*dy*dy*(sx23_in[j][i]-sx22_in[j][i])/DT\
                +3*dy*(2*sx23_in[j][i]-sx22_in[j][i]-sx24_in[j][i])\
                /(DT2)+(3*sx22_in[j][i]-3*sx23_in[j][i]+sx24_in[j][i])\
                /(DT3))*(DT3);

                sx31_in[j][i]=DT2*mo2\
                *(1.0/(DX2))*(u2-u*s2_in[j][i])+2*sx32_in[j][i]\
                -sx33_in[j][i];
                }

            if(snx1==0 && snx2==0 && sny1==0 && sny2==0)
            {
                s3_in[j][i]=(mo2*DT2*(u2-u*s2_in[j][i])*(1.0/(DX2))\
                +DT2*mo2*(1.0/(DY2))*(u1-u*s2_in[j][i]))+2*s2_in[j][i]-s1_in[j][i];
            }
            else
            {
                s3_in[j][i]=sx11_in[j][i]+sx21_in[j][i]+sx31_in[j][i];
            }
            u1=0,u2=0,u=0,ux=0,uy=0;
            }
        }

    swap=s1_in;s1_in=s2_in;s2_in=s3_in;s3_in=swap;
    swap=sx13_in;sx13_in=sx12_in;sx12_in=sx11_in;sx11_in=swap;
    swap=sx24_in;sx24_in=sx23_in;sx23_in=sx22_in;sx22_in=sx21_in;sx21_in=swap;
    swap=sx33_in;sx33_in=sx32_in;sx32_in=sx31_in;sx31_in=swap;

}

void wave2D_unstable::timeslicecopy()
{
    ;
}

void wave2D_unstable_stabletest(int Z, int X, int T, float dz, float dx, float dt, float v, int sz, int sx, float sf, float pmlwide, int dmovie)
{
    wave2D_unstable A(Z,X);
    ofstream outf1;
    outf1.open("testmovie.bin");
    A.dx=dx,A.dy=dz,A.dt=dt;
    A.setvelocity(v);
    A.PML_wide=pmlwide;

    int k;
    float *f;
    f=new float[T];
    wavelet01(f,T,dt,sf);
    for(k=0;k<T;k++)
    {
        A.s2[sz][sx]=A.s2[sz][sx]+f[k];
        A.timeslicecal();
        if(k%dmovie==0)
        {
        datawrite(A.s2, Z, X, outf1);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    cout<<"Have output test movie."<<endl;
}


void wave2D_unstable_stabletest(int Z, int X, int T, float dz, float dx, float dt, float **v, int sz, int sx, float sf, float pmlwide, int dmovie)
{
    wave2D_unstable A(Z,X);
    ofstream outf1;
    outf1.open("testmovie.bin");
    A.dx=dx,A.dy=dz,A.dt=dt;
    matcopy(A.p2,v,Z,X);
    A.PML_wide=pmlwide;

    int k;
    float *f;
    f=new float[T];
    wavelet01(f,T,dt,sf);
    for(k=0;k<T;k++)
    {
        A.s2[sz][sx]=A.s2[sz][sx]+f[k];
        A.timeslicecal();
        if(k%dmovie==0)
        {
        datawrite(A.s2, Z, X, outf1);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    cout<<"Have output test movie."<<endl;
}
///////////////////////////////////////////////////////////


#endif

