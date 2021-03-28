#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
#include "../include/xjc.h"
using namespace std;
using namespace arma;

//float wavelet01(int k, float DT, float hz=30.0);

int main()
{
    int i,j,k,nx,nz,nt;
    float dx,dz,dt,df,dkx,dkz,v,pi(3.1415926);
    
    //定义数组大小，离散化傅里叶参数
    v=3000.0,nx=512,nz=256,nt=51200,dx=10.0,dz=10.0,dt=0.001;
    df=1.0/(dt*nt),dkx=1.0/(nx*dx),dkz=1.0/(nz*dz);
    
    //生成数据，单道子波
    fmat data(nt,nx);
    data.fill(0.0);
    for(i=0;i<nt-500;i++)
    {
        data(i+500,int(nx/2))=wavelet01(i,dt);
    }
    datawrite(data,nt,nx,"data.bin");

    //将数据转换到w-kx域，并输出振幅谱查看
    cx_fmat datafft(nt,nx);
    datafft=fft2(data,nt,nx); //w-kx数据矩阵
    fmat realdatafft(nt,nx);
    realdatafft.fill(0.0);
    for(i=0;i<nt;i++)
    {
        for(j=0;j<nx;j++)
        {
            realdatafft(i,j)=sqrt(datafft(i,j).real()*datafft(i,j).real()\
                +datafft(i,j).imag()*datafft(i,j).imag());
        }
    }
    datawrite(realdatafft,nt,nx,"datafft.bin");

    //将w-kx数据转换到kz-kx域的B矩阵
    int n1,n2;
    float kz,kx,f,xs;
    cx_fmat B(nz,nx),Bifft(nz,nx);
    fmat realBifft(nz,nx);
    B.fill(0.0);
    for(i=1;i<nz/2;i++) //处理正波数
    {
        for(j=0;j<nx/2;j++)
        {
            kx=j*dkx;
            kz=i*dkz;
            f=v*sqrt(kz*kz+kx*kx);
            n1=f/df;
            if(n1<nt/2)
            {
                B(i,j)=B(i,j)+datafft(n1,j);
            }
            xs=v*kz/sqrt(kz*kz+kx*kx);
            B(i,j)=B(i,j)*xs;
        }
    }    

    for(i=1;i<nz/2;i++) //处理负波数
    {
        for(j=nx/2;j<nx;j++)
        {
            kx=-(nx-j)*dkx;
            kz=i*dkz;
            f=v*sqrt(kz*kz+kx*kx);
            n1=f/df;
            if(n1<nt/2)
            {
                B(i,j)=B(i,j)+datafft(n1,j);
            }
            xs=v*kz/sqrt(kz*kz+kx*kx);
            B(i,j)=B(i,j)*xs;
        }
    }    
    Bifft=ifft2(B,nz,nx);
    realBifft=real(Bifft); //傅里叶反变换后，结果为其实部
    datawrite(realBifft,nz,nx,"stolt.bin");

    realBifft=real(B);
    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
        {
            realBifft(i,j)=sqrt(B(i,j).real()*B(i,j).real()\
                +B(i,j).imag()*B(i,j).imag());
        }
    }
    datawrite(realBifft,nz,nx,"B.bin");

    return 0;
}





