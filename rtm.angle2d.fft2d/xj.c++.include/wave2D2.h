/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef WAVE2D_H_H
#define WAVE2D_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "mat.h"
/////////////////////////////////////////////////////////////////////////////////////////////
void wave2Dtest(int Z, int X);
float wavelet01(int k, float DT, float hz=30.0, int delay=100);
template <typename T1> void wavelet01(T1 *w, int N, float DT, float hz=30.0, int delay=100);
float wavelet02(int k, float DT, float hz=30.0, int delay=100);
template <typename T1> void wavelet02(T1 *w, int N, float DT, float hz=30.0, int delay=100);
template <typename T1> float* hilbert1D(T1 *s, int n, float dt);
float Blackman(float n, float N);

class wave2D
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
    
    wave2D();
    wave2D(int x, int y);
    ~wave2D();

    void setvelocity(float v=3000.0);
    void timeslicecal();
    void timeslicecopy();
    void cleardata();
};
///////////////////////////////////////////////////////////////////////////////////////////

wave2D::wave2D()
{
    nx=0;
    ny=0;
    dx=10.0;dy=10.0;dt=0.001;PML_wide=25;suface=1;R=10000000;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

wave2D::wave2D(int z, int x)
{
    nx=x;
    ny=z;
    dx=10.0;dy=10.0;dt=0.001;PML_wide=25;suface=1;R=10000000;
    int i,j;
    s1=new float*[ny];  
    s2=new float*[ny];   
    s3=new float*[ny];  
    sx11=new float*[ny];    
    sx12=new float*[ny];     
    sx13=new float*[ny];     
    sx21=new float*[ny];     
    sx22=new float*[ny];     
    sx23=new float*[ny];     
    sx24=new float*[ny];     
    sx31=new float*[ny];     
    sx32=new float*[ny];     
    sx33=new float*[ny];      
    for(j=0;j<ny;j++)  
        {  
        s1[j]=new float[nx];  
        s2[j]=new float[nx];   
        s3[j]=new float[nx];   
        sx11[j]=new float[nx];   
        sx12[j]=new float[nx];   
        sx13[j]=new float[nx];    
        sx21[j]=new float[nx];   
        sx22[j]=new float[nx];   
        sx23[j]=new float[nx];   
        sx24[j]=new float[nx];   
        sx31[j]=new float[nx];   
        sx32[j]=new float[nx];    
        sx33[j]=new float[nx];   
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

wave2D::~wave2D()
{
    int i,j;

    for(int i=0;i<ny;i++)  
        delete []p2[i]; 
    delete []p2;

    for(int i=0;i<ny;i++)  
       {delete []s1[i];
        delete []s2[i];
        delete []s3[i];
        delete []sx11[i];
        delete []sx12[i];
        delete []sx13[i];
        delete []sx21[i];
        delete []sx22[i];
        delete []sx23[i];
        delete []sx24[i];
        delete []sx31[i];
        delete []sx32[i];
        delete []sx33[i];
       }
    delete []s1;
    delete []s2;
    delete []s3;
    delete []sx11;
    delete []sx12;
    delete []sx13;
    delete []sx21;
    delete []sx22;
    delete []sx23;
    delete []sx24;
    delete []sx31;
    delete []sx32;
    delete []sx33;
    s1=NULL;
    s2=NULL;
    s3=NULL;
    sx11=NULL;
    sx12=NULL;
    sx13=NULL;
    sx21=NULL;
    sx22=NULL;
    sx23=NULL;
    sx24=NULL;
    sx31=NULL;
    sx32=NULL;
    sx33=NULL;
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void wave2D::setvelocity(float v)
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

void wave2D::cleardata()
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

void wave2D::timeslicecal()
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    DX=dx;
    DY=dy;
    DT=dt;
    xshd=PML_wide;
    X=nx;
    Y=ny;
    suface_PML=suface;
 
    float dx,dy,ddx,ddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,n,i1,j1,t3,t4;
    float u1(0),u2(0),u(0),ux(0),uy(0);
    float C3=3/2/(xshd)/DX*log(R)/(xshd)/DX/(xshd)/DX;
    t5=float(Y)/X;

    for(i=5;i<X-5;i++)
        {
        for(j=5;j<Y-5;j++)
            {

            for(n=0;n<5;n++)
                {  
                u=u+2*xs2[n];
                u1=u1+xs2[n]*(s2[j-n-1][i]+s2[j+n+1][i]);
                u2=u2+xs2[n]*(s2[j][i-n-1]+s2[j][i+n+1]);
                }

            for(n=0;n<10;n++)
                {
                if(n<5)
                    {
                    ux=ux+s2[j][i+n-5]*xs1[n];
                    uy=uy+s2[j+n-5][i]*xs1[n];
                    }
                else
                    {
                    ux=ux+s2[j][i+n-4]*xs1[n];
                    uy=uy+s2[j+n-4][i]*xs1[n];
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
                if(j>=Y-xshd-5 && j>t5*i && j>-t5*i+Y)		
                    sny1=j-(Y-xshd-5);
                if(j<=xshd+5 && j<t5*i && j<-t5*i+Y)		
                    sny2=xshd+5-j;		
                }  //��������
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
                dy=p2[j][i]*C3*sny1*sny1*DY*DY;
                ddy=p2[j][i]*C3*2*sny1*DY;	
                }
            if(sny2 !=0)
                {
                dy=p2[j][i]*C3*sny2*sny2*DY*DY;
                ddy=p2[j][i]*C3*2*sny2*DY;	
                }
            if(snx1 !=0)
                {
                dx=p2[j][i]*C3*snx1*snx1*DX*DX;
                ddx=p2[j][i]*C3*2*snx1*DX;	
                }
            if(snx2 !=0)
                {
                dx=p2[j][i]*C3*snx2*snx2*DX*DX;
                ddx=p2[j][i]*C3*2*snx2*DX;	
                }

            if(i>=0.5*(X))
                t3=1;
            else
                t3=-1;

            if(j>=0.5*(Y))
                t4=1;
            else
                t4=-1;

            if(snx1!=0 || snx2!=0)
                {
                sx11[j][i]=p2[j][i]*p2[j][i]\
                *DT*DT*(u2-u*s2[j][i])*(1.0/(DX*DX))\
                -dx*dx*DT*DT*sx12[j][i]+(2*sx12[j][i]\
                -sx13[j][i])+DT*(2*dx*(sx13[j][i]-sx12[j][i]));

                sx21[j][i]=(-p2[j][i]*p2[j][i]\
                *ddx*(1.0/(DX))*(ux*t3)-dx*dx*dx*sx22[j][i]\
                +3*dx*dx*(sx23[j][i]-sx22[j][i])/DT+3*dx*\
                (2*sx23[j][i]-sx22[j][i]-sx24[j][i])/(DT*DT)\
                +(3*sx22[j][i]-3*sx23[j][i]+sx24[j][i])\
                /(DT*DT*DT))*(DT*DT*DT);

                sx31[j][i]=DT*DT*p2[j][i]*p2[j][i]\
                *(1.0/(DX*DX))*(u1-u*s2[j][i])+2*sx32[j][i]\
                -sx33[j][i];
                }
            else
                {
                sx11[j][i]=p2[j][i]*p2[j][i]\
                *DT*DT*(u1-u*s2[j][i])*(1.0/(DX*DX))\
                -dy*dy*DT*DT*sx12[j][i]+(2*sx12[j][i]\
                -sx13[j][i])+DT*(2*dy*(sx13[j][i]-sx12[j][i]));

                sx21[j][i]=(-p2[j][i]*p2[j][i]\
                *ddy*(1.0/(DX))*(uy*t4)-dy*dy*dy*sx22[j][i]\
                +3*dy*dy*(sx23[j][i]-sx22[j][i])/DT\
                +3*dy*(2*sx23[j][i]-sx22[j][i]-sx24[j][i])\
                /(DT*DT)+(3*sx22[j][i]-3*sx23[j][i]+sx24[j][i])\
                /(DT*DT*DT))*(DT*DT*DT);

                sx31[j][i]=DT*DT*p2[j][i]*p2[j][i]\
                *(1.0/(DX*DX))*(u2-u*s2[j][i])+2*sx32[j][i]\
                -sx33[j][i];
                }   

            // if(snx1!=0 || snx2!=0 || sny1!=0 || sny2!=0)
            s3[j][i]=sx11[j][i]+sx21[j][i]+sx31[j][i];
            u1=0,u2=0,u=0,ux=0,uy=0;
            }
        }
}

void wave2D::timeslicecopy()
{
    int j1,i1;
    for(j1=0;j1<ny;j1++)
    {
    for(i1=0;i1<nx;i1++)
        {
        s1[j1][i1]=s2[j1][i1];s2[j1][i1]=s3[j1][i1];
        sx13[j1][i1]=sx12[j1][i1],sx12[j1][i1]=sx11[j1][i1];
        sx33[j1][i1]=sx32[j1][i1],sx32[j1][i1]=sx31[j1][i1];
        sx24[j1][i1]=sx23[j1][i1],sx23[j1][i1]=sx22[j1][i1];
        sx22[j1][i1]=sx21[j1][i1];
        }
    }
}

void wave2Dtest(int Z, int X, int T)
{
    wave2D A(Z,X);
    ofstream outf1;
    outf1.open("testmovie.bin");

    int k;
    float f;
    for(k=0;k<T;k++)
    {
        f=wavelet01(k,A.dt);
        A.s2[Z/2][X/2]=A.s2[Z/2][X/2]+f;
        A.timeslicecal();
        A.timeslicecopy();
        datawrite(A.s3, Z, X, outf1);
    }
    outf1.close();
    cout<<"Have output test movie."<<endl;
}

///////////////////////////////////////////////////////////

float wavelet02(int k, float DT, float hz, int delay)
{
    float pi(3.1415926);
	float f,det;
    det=delay*DT;

    f=(pi)*(pi)*hz*hz*(k*DT-det)*\
    exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
    *(3.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));

	return f;
}

template<typename T1>
void wavelet02(T1 *w, int N, float DT, float hz, int delay)
{
    float pi(3.1415926);
	float f,det;
    det=delay*DT;
    int k;

    for(k=0;k<N;k++)
    {
        f=(pi)*(pi)*hz*hz*(k*DT-det)*\
        exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
        *(3.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));
        w[k]=f;
    }
}

float wavelet01(int k, float DT, float hz, int delay)
{
    float pi(3.1415926);
	float f,det;
    det=delay*DT;

    f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
    *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));

	return f;
}

template<typename T1>
void wavelet01(T1 *w, int N, float DT, float hz, int delay)
{
    float pi(3.1415926);
	float f,det;
    det=delay*DT;
    int k;

    for(k=0;k<N;k++)
    {
        f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
        *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));
        w[k]=f;
    }
}

template<typename T1>
float* hilbert1D(T1 *s, int n, float dt)
{
    float *h, pi(3.1415926);
    int i,j,z,hn;
    hn=3/dt;
    z=int(hn/2);
    h=new float[hn];
    for(i=0;i<hn;i++)
    {
        if((i-z)!=0)
        {
            h[i]=1/(pi*dt*(i-z));
        }
        else
        {
            h[i]=0.0;
        }   
    }

    float *h2, *s2, *sh, *sh2;
    s2=new float[hn+n];
    h2=new float[hn+n];
    sh=new float[hn+n];
    sh2=new float[hn+n];
    for(i=0;i<hn+n;i++)
    {
        s2[i]=0.0;
        h2[i]=0.0;
        sh[i]=0.0;
        sh2[i]=0.0;
        if(i<n)
        {
            s2[i]=s[i];
        }
        if(i<hn)
        {
            h2[i]=h[i];
        }
    }

    for(i=0;i<hn+n;i++)
    {
        for(j=0;j<=i;j++)
        {
            sh[i]+=s2[i-j]*h2[j];
        }
        if(i>=(z))
        {
            sh2[i-z]=sh[i]/(1.0/dt);
        }
    }
    return sh2;
}

float Blackman(float n, float N)
{
    float xs, pi(3.1415926);
    xs=0.42-0.5*cos(2*pi*n/N/2)+0.08*cos(4*pi*n/N/2);
    return xs;
}

#endif

