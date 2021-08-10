/************update in 2021.01.07**************/

    
    
/*********************************************/

#ifndef RADON_BEAMFORMING_HPP
#define RADON_BEAMFORMING_HPP

#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace arma;
#include "../xjc.h"


void tx_to_fk2d(fmat realdatafk, cx_fmat datafk, fmat data);
void tx_to_fk2d(struct linerradon2d & par);
void tx_to_fx2d(cx_fmat datafx, fmat data);
void tx_to_fx2d(struct linerradon2d & par);
void smoothdig(fmat & dig,int l,int n);
fmat get_blackman_win2d(int n1, int n2, float w1, float w2);

//从复数矩阵提取振幅谱
fmat amplitude_arma(cx_fmat & m, int n1, int n2)
{
    fmat a(n1,n2);
    int i,j;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            a(i,j)=pow(m(i,j).real(),2)+pow(m(i,j).imag(),2);
            a(i,j)=sqrt(a(i,j));
        }
    }
    return a;
}

/*beamforming/liner_radon变换 传递的参数：
nz（处理数据体Z方向采样点）、nx（处理数据体X方向采样点）、
nf（数据体频率域变换的频率采样）、np（Radon变换倾角采样点）、
par4（RLS_beamforming的迭代次数）、par6（我也忘了这参数干嘛的）、
dz（数据Z方向采样间隔）、dx（数据X方向采样间隔）、
df（数据变换后频率采样间隔）、dp（数据变换后倾角采样间隔）、
p0（数据变换后的中心倾角值）、data（原始数据）、
realdataTP（beamforming/liner_radon变换得到的数据实部）、
realrebuild（beamforming/liner_radon反变换重建的数据实部）、
realdatafft（beamforming/liner_radon变换得到的中间矩阵实部，用于检查算法）、
dataTP（beamforming/liner_radon变换得到的复数数据体）、
rebuild（beamforming/liner_radon反变换重建的复数数据体）、
dig_n（对角加权值）、
par1,par2,par3,par5（这几个par曾用于调试优化结果）、

allAreal[99],allAimag[99]
（是文件路径，用于存放一个转换矩阵A；在处理多批数据时需要重复计算A，
且A较大不好直接存在内存里，故将其存入文件以备读取使用）
*/
struct linerradon2d
{
    int nz,nx,nf,np,par4,par6,nf1,nf2;
    float dz,dx,df,dp,p0,pcenter;
    fmat data,realdataTP,realrebuildtx,realdatafk;
    cx_fmat datafP,dataTP,rebuildfp,rebuildfx,\
        rebuildtx,datafx,datafk;
    float dig_n,par1,par2,par3,par5;
    char allAreal[99],allAimag[99];

    //cx_fmat Acol1,Arow1;
};

//设置一组常用的初始化参数
void beamforming_parset(int nz,int nx, int nf,\
    struct linerradon2d & par)
{
    par.nx=nx,par.nz=nz,par.nf=nf;
    par.par3=1.0;par.par4=1,par.par6=0;
    par.pcenter=0.0;
    par.dx=(10),par.dz=0.001;
    par.df=1.0/par.dz/par.nf;
    par.np=par.nz,par.dp=2*par.dz/par.nx/par.dx;
    par.p0=-par.dp*int(par.np/2)+par.pcenter;
    par.dig_n=(0.01),par.par1=(0.05),par.par2=(0.03);
    par.par5=(0.01);
    par.nf1=0;par.nf2=par.nf/2;
    par.allAreal[0]='\0';
    strcat(par.allAreal,"./allAreal.bin");
    par.allAimag[0]='\0';
    strcat(par.allAimag,"./allAimag.bin");
}

void beamforming_parupdate(struct linerradon2d & par)
{
    par.p0=-par.dp*int(par.np/2)+par.pcenter;
    par.data.zeros(par.nz,par.nx);
    par.datafk.zeros(par.nf,par.nx);
    par.datafP.zeros(par.nf,par.np);
    par.datafx.zeros(par.nf,par.nx);
    par.dataTP.zeros(par.nz,par.np);
    par.realdataTP.zeros(par.nz,par.np);
    par.realdatafk.zeros(par.nf,par.nx);
    par.rebuildfp.zeros(par.nf,par.np);
    par.rebuildfx.zeros(par.nf,par.nx);
    par.rebuildtx.zeros(par.nz,par.nx);
    par.realrebuildtx.zeros(par.nz,par.nx);
}

void beamforming_cleardata(struct linerradon2d & par)
{
    par.datafk.fill(0.0);
    par.datafP.fill(0.0);
    par.datafx.fill(0.0);
    par.dataTP.fill(0.0);
    par.realdataTP.fill(0.0);
    par.realdatafk.fill(0.0);
    par.rebuildfx.fill(0.0);
    par.rebuildtx.fill(0.0);
    par.realrebuildtx.fill(0.0);
}

//预备计算变换需要的A矩阵，便于后续使用
void linerandontransmat(struct linerradon2d & par)
{
    int i,j,k;
    float w,pi(3.1415926),f,dx(par.dx),\
        dp(par.dp),p1(par.p0);
    int nx(par.nx),np(par.np);
    fmat X(nx,1), P(1,np), forA(nx,np);
    cx_fmat forAw(nx,np), A(nx,np);
    ofstream outreal,outimag;
    outreal.open(par.allAreal);
    outimag.open(par.allAimag);
    for(k=0;k<par.nf;k++)  
    {
        w=2.0*pi*par.df*k;  
        for(i=0;i<nx;i++)
        {
            X(i,0)=(i)*dx-(nx*dx)/2.0;
        }
            
        for(i=0;i<np;i++)
        {
            P(0,i)=i*dp+p1;
        }
        forA=X*P;   

        for(i=0;i<nx;i++)
        {
            for(j=0;j<np;j++)
            {
                (forAw(i,j).real(0.0));
                (forAw(i,j).imag(0.0));
                (forAw(i,j).imag(forA(i,j)*w));
                //forAw(i,j)=forAw(i,j);
            }
        }
        A=exp(forAw);
        datawrite(forA=real(A),nx,np,outreal);
        datawrite(forA=imag(A),nx,np,outimag);
    }
    outreal.close(),outimag.close();
}

//线性Radon变换
void linerradon(struct linerradon2d & par)
{
    int i,j,k,k2;//cout<<"ok"<<endl;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),\
        dp(par.dp),p1(par.p0),dz(par.dz);
    int nx(par.nx),np(par.np),nf(par.nf);

    cx_fmat A(nx,np);

    ofstream outf;
    ifstream infreal,infimag;
    infreal.open(par.allAreal);
    infimag.open(par.allAimag);

    tx_to_fx2d(par);
    
    infimag.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    infreal.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    for(k=par.nf1;k<par.nf2;k++)  
    {
        A.set_real(dataread(nx,np,infreal));
        A.set_imag(dataread(nx,np,infimag));
        par.datafP.row(k)=(A.t()*par.datafx.row(k).st()).st();
    }
    infimag.close(),infreal.close();
    
    for(i=0;i<np;i++)
    {
        par.dataTP.col(i)=ifft(par.datafP.col(i),par.nz);
    }  
    par.realdataTP=real(par.dataTP);
}

//重建数据，重建之前可以对数据按倾角去噪等
void rebuildsignal(struct linerradon2d & par)
{
    int i,j,k,k2;//cout<<"ok"<<endl;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),\
        dp(par.dp),p1(par.p0),dz(par.dz);
    int nx(par.nx),np(par.np),nf(par.nf);

    cx_fmat A(nx,np);

    ofstream outf;
    ifstream infreal,infimag;
    infreal.open(par.allAreal);
    infimag.open(par.allAimag);

    for(i=0;i<np;i++)
    {
        par.rebuildfp.col(i)=fft(par.dataTP.col(i),par.nf);
    } 
    
    infimag.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    infreal.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    for(k=par.nf1;k<par.nf2;k++)  
    {
        A.set_real(dataread(nx,np,infreal));
        A.set_imag(dataread(nx,np,infimag));
        par.rebuildfx.row(k)=(A*par.rebuildfp.row(k).st()).st(); 
    }

    for(i=0;i<nx;i++)
    {
        par.rebuildtx.col(i)=ifft(par.rebuildfx.col(i),par.nz);
    }      
    par.realrebuildtx=real(par.rebuildtx);  
    par.realrebuildtx/=(par.nf/2.0);
    infimag.close(),infreal.close();
}

//RLS_beamforming变换函数
void beamforming_radon(struct linerradon2d & par)
{
    //linerradon(par);
    int i,j,k,k2;
    float w,pi(3.1415926),maxpower;
    float df(par.df),dx(par.dx),dig_n(par.dig_n),\
        dp(par.dp),p1(par.p0),dz(par.dz);
    int nx(par.nx),np(par.np),nf(par.nf);

    fmat radon(np,1),realTP(nf,np);
    cx_fmat S(np,1);
    cx_fmat A(nx,np);
    cx_fmat digA(np,np);
    fmat realdigA(np,np);
    ifstream infreal,infimag;

    linerradon(par);

for(k2=0;k2<par.par4;k2++)
{

    for(j=0;j<np;j++)
    {
        for(i=0;i<nf;i++)
        {
            realTP(i,j)=real(par.dataTP(i,j))\
                *real(par.dataTP(i,j));
        }
    }

    digA.fill(0.0);
    radon.fill(0.0);
    
//根据前一次迭代结果调整获取对角权重
    for(i=0;i<np;i++)
    {
        radon(i,0)=sum(realTP.col(i));
        if(i<np*par.par1)
        radon(i,0)=radon(i,0)*exp(-par.par2*(np*par.par1-i));
        if(i>np*(1-par.par1))
        radon(i,0)=radon(i,0)*exp(-par.par2*(i-np*(1-par.par1)));
    }
    maxpower=max(radon.col(0));
    
    for(i=0;i<np;i++)
    { 
        radon(i,0)=radon(i,0)/maxpower;
    }
    smoothdig(radon,np,par.par6);
    maxpower=max(radon.col(0));
    
    for(i=0;i<np;i++)
    { 
        radon(i,0)=radon(i,0)/maxpower;
        digA(i,i).real(par.par3*dig_n/(radon(i,0)+par.par5));
        digA(i,i).real(dig_n);
    }

///////////////////////////

    infreal.open(par.allAreal);
    infimag.open(par.allAimag);
    infimag.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    infreal.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    for(k=par.nf1;k<par.nf2;k++)
    {
        A.set_real(dataread(nx,np,infreal));
        A.set_imag(dataread(nx,np,infimag));

        //反演求解
        //S=inv_sympd(A.t()*A+digA)*A.t()*par.datafx.row(k).st();
        S=inv(A.t()*A+digA)*A.t()*par.datafx.row(k).st();
        par.datafP.row(k)=S.col(0).st(); 

        for(i=0;i<np;i++)
        {
            radon(i,0)+=abs(par.datafP(k,i));
            if(i<np*par.par1)
            radon(i,0)=radon(i,0)*exp(-par.par2*(np*par.par1-i));
            if(i>np*(1-par.par1))
            radon(i,0)=radon(i,0)*exp(-par.par2*(i-np*(1-par.par1)));
        }
    }
    infimag.close();
    infreal.close();
    //datawrite(radon/=max(radon.col(0)),np,1,"dig.low.bin");

    for(i=0;i<np;i++)
    {
        par.dataTP.col(i)=ifft(par.datafP.col(i),par.nz);
    }            
    par.realdataTP=real(par.dataTP);
}  
}
/////////////////////////
void beamforming_lowanti(struct linerradon2d & par)
{
    //linerradon(par);
    int i,j,k,k2;
    float w,pi(3.1415926),maxpower;
    float df(par.df),dx(par.dx),dig_n(par.dig_n),\
        dp(par.dp),p1(par.p0),dz(par.dz);
    int nx(par.nx),np(par.np),nf(par.nf);

    fmat radon(np,1),realTP(nf,np);
    cx_fmat S(np,1);
    cx_fmat A(nx,np);
    cx_fmat digA(np,np);
    fmat realdigA(np,np);
    ifstream infreal,infimag;
    
{
    digA.fill(0.0);
    radon.fill(0.0);
//由低频到高频反演
    tx_to_fx2d(par);

    infreal.open(par.allAreal);
    infimag.open(par.allAimag);
    infimag.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    infreal.seekg(nx*np*par.nf1*sizeof(float),ios::beg);
    for(k=par.nf1;k<par.nf2;k++)
    {
        maxpower=max(radon.col(0));        
        for(i=0;i<np;i++)
        { 
            digA(i,i).real(par.par3*dig_n/\
                ((radon(i,0)/(maxpower+par.par5))+par.par5));
        }

        A.set_real(dataread(nx,np,infreal));
        A.set_imag(dataread(nx,np,infimag));

        //反演求解
        //S=inv_sympd(A.t()*A+digA)*A.t()*par.datafx.row(k).st();
        S=inv(A.t()*A+digA)*A.t()*par.datafx.row(k).st();
        par.datafP.row(k)=S.col(0).st(); 

        for(i=0;i<np;i++)
        {
            radon(i,0)+=abs(par.datafP(k,i));
            if(i<np*par.par1)
            radon(i,0)=radon(i,0)*exp(-par.par2*(np*par.par1-i));
            if(i>np*(1-par.par1))
            radon(i,0)=radon(i,0)*exp(-par.par2*(i-np*(1-par.par1)));
        }
    }
    infimag.close();
    infreal.close();
    //datawrite(radon/=max(radon.col(0)),np,1,"dig.low.bin");

    for(i=0;i<np;i++)
    {
        par.dataTP.col(i)=ifft(par.datafP.col(i),par.nz);
    }            
    par.realdataTP=real(par.dataTP);
}  
}

//some other function
void smoothdig(fmat & dig,int l,int n)
{
    fmat d1(l,1),d2(l,1);
    d1=dig,d2=dig;
    int j,k;
    for(j=0;j<n;j++)
    {
        for(k=1;k<l-1;k++)
        {
            d2(k,0)=(0.5*d1(k-1,0)+0.5*d1(k+1,0)+d1(k,0))/2;
        }
        d1=d2;
    }
    dig=d2;
}

void tx_to_fx2d(struct linerradon2d & par)
{
    int i;
    for(i=0;i<par.nx;i++)
    {
        par.datafx.col(i)=fft(par.data.col(i),par.nf);
    } 
}

void tx_to_fx2d(cx_fmat datafx, fmat data)
{
    int i;
    for(i=0;i<data.n_cols;i++)
    {
        datafx.col(i)=fft(data.col(i),data.n_rows);
    } 
}

void tx_to_fk2d(struct linerradon2d & par)
{
    int i,j,n1,n2,n11,n22;
    n1=par.datafk.n_rows;
    n11=n1/2;
    n2=par.datafk.n_cols;
    n22=n2/2;
    par.datafk=fft2(par.data,n1,n2);
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i<n11 && j<n22)
            {
                par.realdatafk(n11+i,n22+j)=abs(par.datafk(i,j));
            }
            if(i<n11 && j>n22)
            {
                par.realdatafk(n11+i,j-n22)=abs(par.datafk(i,j));
            }
            if(i>n11 && j<n22)
            {
                par.realdatafk(i-n11,j+n22)=abs(par.datafk(i,j));
            }
            if(i>n11 && j>n22)
            {
                par.realdatafk(i-n11,j-n22)=abs(par.datafk(i,j));
            }
        }
    }    
}

void tx_to_fk2d(fmat realdatafk, cx_fmat datafk, fmat data)
{
    int i,j,n1,n2,n11,n22;
    n1=datafk.n_rows;
    n11=n1/2;
    n2=datafk.n_cols;
    n22=n2/2;
    datafk=fft2(data,n1,n2);
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i<n11 && j<n22)
            {
                realdatafk(n11+i,n22+j)=abs(datafk(i,j));
            }
            if(i<n11 && j>n22)
            {
                realdatafk(n11+i,j-n22)=abs(datafk(i,j));
            }
            if(i>n11 && j<n22)
            {
                realdatafk(i-n11,j+n22)=abs(datafk(i,j));
            }
            if(i>n11 && j>n22)
            {
                realdatafk(i-n11,j-n22)=abs(datafk(i,j));
            }
        }
    }    
}
/*
float Blackman(float n, float N)
{
    float xs, pi(3.1415926);
    xs=0.42-0.5*cos(2*pi*n/N/2)+0.08*cos(4*pi*n/N/2);
    return xs;
}*/
fmat get_blackman_win2d(int n1, int n2, float w1, float w2)
{
    fmat win(n1,n2);
    win.fill(1.0);
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i<=w2)
            {
                n=i;
                n=Blackman(n,w2);
                win(i,j)*=n;
            }
            if(j<=w1)
            {
                n=j;
                n=Blackman(n,w1);
                win(i,j)*=n;
            }
            if(i>=n1-w2-1)
            {
                n=n1-1-i;
                n=Blackman(n,w2);
                win(i,j)*=n;
            }
            if(j>=n2-w1-1)
            {
                n=n2-1-j;
                n=Blackman(n,w1);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//RLS_beamforming变换函数
void beamforming_large(struct linerradon2d & par)
{
    int i,j,k,k2;//cout<<"ok"<<endl;
    float w,pi(3.1415926),maxpower;
    float df(par.df),dx(par.dx),dig_n(par.dig_n),\
        dp(par.dp),p1(par.p0),DT(par.dz);
    int NX(par.nx),Np(par.np),NF(par.nf);
    fmat realTP(NF,Np);
    fmat radon(Np,1);
    fmat radon2(Np,1);
    radon2.fill(0.0);
    fmat data(1,1);
    cx_fmat hessian(Np,Np);
    cx_fmat digA(Np,Np);
    fmat realdigA(Np,Np);
    cx_fmat dataP(NF,Np);
    cx_fmat dataTP(NF,Np);
    par.realdataTP.resize(NF,Np);
    cx_fmat datafft(NF,NX);
    cx_fmat A(NX,Np);
    cx_fmat S(Np,1);
    fmat X(NX,1);
    fmat P(1,Np);
    fmat forA(NX,Np);
    cx_fmat forAw(NX,Np);

    data.copy_size(par.data);
    data=par.data;
    ofstream outf;
    ifstream infreal,infimag;
    infreal.open(par.allAreal);
    infimag.open(par.allAimag);
    float outdata;

//变换是在频率域进行的，需要先获得数据单频片
    for(i=0;i<NX;i++)
    {
        datafft.col(i)=fft(data.col(i),NF);
    } 
    df=1.0/DT/NF;

cout<<"now is running: 1 -- "<<endl;
    
    dataP.fill(0.0);    
    digA.fill(0.0);
    digA.diag()+=dig_n;  
 
    infimag.seekg(0,ios::beg);
    infreal.seekg(0,ios::beg);

//变换只使用正半频率
for(k=0;k<par.nf2;k++)  
    {
     w=2.0*pi*par.df*k;  
        for(i=0;i<NX;i++)
        {
            X(i,0)=(i)*dx-(NX*dx)/2.0;
        }
            
        for(i=0;i<Np;i++)
            {
            P(0,i)=i*dp+p1;
            }
        forA=X*P;   //准备好A矩阵中的p，dx部分
    //  forA.save("A.txt",raw_ascii);
        
        for(i=0;i<NX;i++)
            {
                for(j=0;j<Np;j++)
                {
                    (forAw(i,j).real(0.0));
                    (forAw(i,j).imag(0.0));
                    (forAw(i,j).imag(forA(i,j)*w));
                    //forAw(i,j)=forAw(i,j);
                }
            }
        A=exp(forAw);
        dataP.row(k)=(A.t()*datafft.row(k).st()).st();
    }
    
    for(i=0;i<Np;i++)
    {
        dataTP.col(i)=ifft(dataP.col(i),NF);
    }  
    par.dataTP.copy_size(dataTP);
    par.dataTP=dataTP;
    realTP=real(par.dataTP);
    datawrite(realTP,NF,Np,"radonTP.bin"); //输出线性Radon变换结果

outf.open("dig.bin");
//开始RLS_beamforming迭代
for(k2=0;k2<par.par4;k2++)
{
    for(j=0;j<Np;j++)
    {
        for(i=0;i<NF;i++)
        {
            realTP(i,j)=real(dataTP(i,j))*real(dataTP(i,j));
        }
    }

cout<<"now is running: 2 -- "<<endl;
    digA.fill(0.0);
    
//根据前一次迭代结果调整获取对角权重
    for(i=0;i<Np;i++)
    {
        radon(i,0)=sum(realTP.col(i));
        if(i<Np*par.par1)
        radon(i,0)=radon(i,0)*exp(-par.par2*(Np*par.par1-i));
        if(i>Np*(1-par.par1))
        radon(i,0)=radon(i,0)*exp(-par.par2*(i-Np*(1-par.par1)));
    }
    maxpower=max(radon.col(0));
    
    for(i=0;i<Np;i++)
    { 
        radon(i,0)=radon(i,0)/maxpower;
    }
    smoothdig(radon,Np,par.par6);
    maxpower=max(radon.col(0));
    
    for(i=0;i<Np;i++)
    { 
        radon(i,0)=radon(i,0)/maxpower;
        digA(i,i).real(par.par3*dig_n/(radon(i,0)+par.par5));
        //digA(i,i).imag(par.par3*dig_n/(radon(i,0)+dig_n));
        //digA(i,i).real(par.par3-radon(i,0)+dig_n);
        //digA(i,i).imag(par.par3-radon(i,0)+dig_n);
        //outdata=real(digA(i,i));
        //outf.write((char*)&outdata,sizeof(outdata));
        //outdata=i;
        //outf.write((char*)&outdata,sizeof(outdata));
    }
    datawrite(radon,Np,1,outf);
    maxpower=0;
    for(i=0;i<Np;i++)
    { 
        maxpower+=pow(radon(i,0)-radon2(i,0),2);
    }
    if(sqrt(maxpower)<0.1)
    {
        cout<<"times is: "<<k2<<endl;
        break;
    }
    radon2=radon;

    infimag.seekg(0,ios::beg);
    infreal.seekg(0,ios::beg);
    for(k=0;k<par.nf2;k++)
    {
        w=2.0*pi*par.df*k;  
        for(i=0;i<NX;i++)
        {
            X(i,0)=(i)*dx-(NX*dx)/2.0;
        }
            
        for(i=0;i<Np;i++)
            {
            P(0,i)=i*dp+p1;
            }
        forA=X*P;   //准备好A矩阵中的p，dx部分
    //  forA.save("A.txt",raw_ascii);
        
        for(i=0;i<NX;i++)
            {
                for(j=0;j<Np;j++)
                {
                    (forAw(i,j).real(0.0));
                    (forAw(i,j).imag(0.0));
                    (forAw(i,j).imag(forA(i,j)*w));
                    //forAw(i,j)=forAw(i,j);
                }
            }
        A=exp(forAw);
        hessian=A.t()*A;

        //反演求解
        S=inv(hessian+digA)*A.t()*datafft.row(k).st();

        dataP.row(k)=S.col(0).st();  //???��������ĸ�ֵ�������ģ�����
        //cout<<"now is running: "<<k<<endl;
        /*
        if(k==int(30.0/par.df))
        {
            for(i=0;i<Np;i++)
            {
                for(j=0;j<Np;j++)
                {
                    realdigA(i,j)=pow(hessian(i,j).real(),2)+\
                        pow(hessian(i,j).imag(),2);
                }
            }
            datawrite(realdigA,Np,Np,"digA.bin");
        }*/
    }
    
    for(i=0;i<Np;i++)
    {
        dataTP.col(i)=ifft(dataP.col(i),NF);
    }            
    par.dataTP.copy_size(dataTP);
    par.dataTP=dataTP;
    par.realdataTP=real(dataTP);
    realTP=real(dataTP);
    //datawrite(realTP,NF,Np,"mediumTP.bin");
}
    outf.close();
    infimag.close();
    infreal.close();
}

//_large表示数据较多，会将A存入硬盘，数据少时可以不用该函数
void linerradon_large(struct linerradon2d & par)
{
    int i,j,k,k2;//cout<<"ok"<<endl;
    float w,pi(3.1415926),maxpower;
    float df(par.df),dx(par.dx),dig_n(par.dig_n),\
        dp(par.dp),p1(par.p0),DT(par.dz);
    int NX(par.nx),Np(par.np),NF(par.nf);
    fmat X(NX,1);
    fmat P(1,Np);
    fmat forA(NX,Np);
    cx_fmat forAw(NX,Np);
    fmat realTP(NF,Np);
    fmat radon(Np,1);
    fmat radon2(Np,1);
    radon2.fill(0.0);
    fmat data(1,1);
    cx_fmat digA(Np,Np);
    fmat realdigA(Np,Np);
    cx_fmat dataP(NF,Np);
    cx_fmat dataTP(NF,Np);
    cx_fmat datafft(NF,NX);
    cx_fmat A(NX,Np);
    cx_fmat S(Np,1);

    data.copy_size(par.data);
    data=par.data;
    ofstream outf;
    float outdata;

    for(i=0;i<NX;i++)
    {
        datafft.col(i)=fft(data.col(i),NF);
    } 
    df=1.0/DT/NF;

cout<<"now is running: 1 -- "<<endl;
cout<<dp<<" | "<<p1<<" | "<<Np<<" | "<<endl;
    
    dataP.fill(0.0);    
    digA.fill(0.0);
    digA.diag()+=dig_n;  
 
    for(k=0;k<par.nf2;k++)  
    {
        w=2.0*pi*par.df*k;  
        for(i=0;i<NX;i++)
        {
            X(i,0)=(i)*dx-(NX*dx)/2.0;
        }
            
        for(i=0;i<Np;i++)
            {
            P(0,i)=i*dp+p1;
            }
        forA=X*P;   //准备好A矩阵中的p，dx部分
    //  forA.save("A.txt",raw_ascii);
        
        for(i=0;i<NX;i++)
            {
                for(j=0;j<Np;j++)
                {
                    (forAw(i,j).real(0.0));
                    (forAw(i,j).imag(0.0));
                    (forAw(i,j).imag(forA(i,j)*w));
                    //forAw(i,j)=forAw(i,j);
                }
            }
        A=exp(forAw);
        dataP.row(k)=(A.t()*datafft.row(k).st()).st();
    }
    
    for(i=0;i<Np;i++)
    {
        dataTP.col(i)=ifft(dataP.col(i),NF);
    }  
    par.dataTP.copy_size(dataTP);
    par.dataTP=dataTP;
    par.realdataTP.resize(NF,Np);
    par.realdataTP=real(dataTP);
    datawrite(realTP=real(par.dataTP),NF,Np,"radonTP.bin");
}

//重建较多数据时使用，同上
void rebuildsignal_large(struct linerradon2d & par)
{
    int i,j,k;
    float w,pi(3.1415926);
    float df(par.df),dx(par.dx),dig_n(par.dig_n),\
        dp(par.dp),p1(par.p0),DT(par.dz);
    int NX(par.nx),Np(par.np),NF(par.nf);
    fmat X(NX,1);
    fmat P(1,Np);
    fmat forA(NX,Np);
    cx_fmat datarebuild(NF,NX);
    cx_fmat forAw(NX,Np);
    cx_fmat rebuild_X(NX,1);           
    cx_fmat A(NX,Np);
    cx_fmat dataP(NF,Np);
    cx_fmat rebuild(NF,NX);
    rebuild.fill(0.0);

    df=1.0/DT/NF;
    for(i=0;i<Np;i++)
    {
        dataP.col(i)=fft(par.realdataTP.col(i),NF);
    } 
    /*****beam_forming_rebuild*****/
    cout<<"now is running: 3 -- "<<endl;
    ofstream outf;
    
    for(k=0;k<par.nf2;k++)
    {    
        w=2.0*pi*par.df*k;  
        for(i=0;i<NX;i++)
        {
            X(i,0)=(i)*dx-(NX*dx)/2.0;
        }
            
        for(i=0;i<Np;i++)
            {
            P(0,i)=i*dp+p1;
            }
        forA=X*P;   //准备好A矩阵中的p，dx部分
        for(i=0;i<NX;i++)
            {
                for(j=0;j<Np;j++)
                {
                    (forAw(i,j).real(0.0));
                    (forAw(i,j).imag(0.0));
                    (forAw(i,j).imag(forA(i,j)*w));
                    //forAw(i,j)=forAw(i,j);
                }
            }
        A=exp(forAw);
        rebuild_X=A*dataP.row(k).st();
        rebuild.row(k)=rebuild_X.col(0).st();  
    }

    for(i=0;i<NX;i++)
    {
        datarebuild.col(i)=ifft(rebuild.col(i),NF);
    }      
    par.rebuildtx.copy_size(datarebuild);
    par.rebuildtx=datarebuild; 
    par.realrebuildtx.resize(NF,NX);
    par.realrebuildtx=real(datarebuild);  
}

#endif
