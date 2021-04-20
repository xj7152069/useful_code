/************update in 2020.12.16**************/
/*
    make some improvement, add two par into 
    function getTP; add a new function linerandomtransmat.
    
***********************************************/

#ifndef RADONBEAMFORMING_H_H
#define RADONBEAMFORMING_H_H

#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include<iomanip>
#include<math.h>
using namespace arma;
using namespace std;
#include <xj.c++.include/xjc.h>

fmat CreatePmat(int NX, int Np, float dx, float dp, float p1)
{
    int i;
    fmat X(NX,1);
    fmat P(1,Np);
    fmat forA(NX,Np);

    for(i=0;i<NX;i++)
        {
        X(i,0)=(i)*dx-(NX*dx)/2.0;
        }
        
    for(i=0;i<Np;i++)
        {
        P(0,i)=i*dp+p1;
        }
     forA=X*P;   //׼����A�����е�p��dx����
   //  forA.save("A.txt",raw_ascii);
   return forA;
}

cx_fmat getTP(fmat data, fmat forA, int Np, int NX, int NF, int Nf, \
    float DT=0.0005, float dig_n=0.01, float par1=0.2, float par2=0.03)
{
    int i,j,k;//cout<<"ok"<<endl;
    float w,pi(3.1415926),df,maxpower;
    fmat realTP(NF,Np);
    cx_fmat forAw(NX,Np);
    cx_fmat digA(Np,Np);
    cx_fmat dataP(NF,Np);
    cx_fmat dataTP(NF,Np);
    cx_fmat datafft(NF,NX);
    cx_fmat A(NX,Np);
    cx_fmat S(Np,1);
    fmat radon(Np,1);
    ofstream outf;
    float outdata;

    for(i=0;i<NX;i++)
    {
    datafft.col(i)=fft(data.col(i),NF);
    } 
    df=1.0/DT/NF;

cout<<"now is running: 1 -- "<<endl;
    
    dataP.fill(0.0);    
    digA.fill(0.0);
    digA.diag()+=dig_n;  
 
     for(k=0;k<Nf;k++)  
        {
         w=2.0*pi*(k+1)*df;  
         for(i=0;i<NX;i++)
         {for(j=0;j<Np;j++)
            {
            (forAw(i,j).real(0.0));
            (forAw(i,j).imag(0.0));
            (forAw(i,j).imag(forA(i,j)*w));
            }
         }
         A=exp(forAw);
         dataP.row(k)=(A.t()*datafft.row(k).t()).t();
        }
    
    for(i=0;i<Np;i++)
        {
        dataTP.col(i)=ifft(dataP.col(i),NF);
        }  

    for(j=0;j<Np;j++)
        {for(i=0;i<NF;i++)
                {
                realTP(i,j)=real(dataTP(i,j));
                }
        }
    datawrite(realTP,NF,Np,"media.TP.bin");

    for(j=0;j<Np;j++)
        {for(i=0;i<NF;i++)
                {
                realTP(i,j)=real(dataTP(i,j))*real(dataTP(i,j));
                }
        }
 
cout<<"now is running: 2 -- "<<endl;
    digA.fill(0.0);
    
    for(i=0;i<Np;i++)
        {
        radon(i,0)=sum(realTP.col(i));
        
        if(i<Np*par1)
        radon(i,0)=radon(i,0)*exp(-0.03*(Np*par1-i));
        if(i>Np*(1-par1))
        radon(i,0)=radon(i,0)*exp(-0.03*(i-Np*(1-par1)));
        }
    maxpower=max(radon.col(0));
    outf.open("radon.bin");
    for(i=0;i<Np;i++)
        { 
        radon(i,0)=radon(i,0)/maxpower;
        outdata=radon(i,0);
        outf.write((char*)&outdata,sizeof(outdata));
        outdata=i;
        outf.write((char*)&outdata,sizeof(outdata));
        digA(i,i)=dig_n/(radon(i,0)+dig_n);
        }
    //outf.close();
      
    for(k=0;k<Nf;k++)
        {
         w=2.0*pi*(k+1)*df;  //�����Ƶ��
         for(i=0;i<NX;i++)
         {for(j=0;j<Np;j++)
            {
            (forAw(i,j).real(0.0));
            (forAw(i,j).imag(0.0));
            (forAw(i,j).imag(forA(i,j)*w));
            //forAw(i,j)=forAw(i,j);
            }
         }
         A=exp(forAw);
         S=inv(A.t()*A+digA)*A.t()*datafft.row(k).st();
         dataP.row(k)=S.col(0).st();  //???��������ĸ�ֵ�������ģ�����
         //cout<<"now is running: "<<k<<endl;
        }
    
    for(i=0;i<Np;i++)
        {
        dataTP.col(Np-i-1)=ifft(dataP.col(i),NF);
        }            
    return dataTP;
}

cx_fmat getRebuildData(cx_fmat dataTP, fmat forA, int Np, int NX, int NF, int Nf, float DT=0.0005)
{
    int i,j,k;
    float w,pi(3.1415926),df;
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
    dataP.col(i)=fft(dataTP.col(i),NF);
    } 
    /*****beam_forming_rebuild*****/
    cout<<"now is running: 3 -- "<<endl;
    
     for(k=0;k<Nf;k++)
        {    
         w=2.0*pi*(k+1)*df;  //�����Ƶ��
         for(i=0;i<NX;i++)
         {for(j=0;j<Np;j++)
            {
            (forAw(i,j).real(0.0));
            (forAw(i,j).imag(0.0));
            (forAw(i,j).imag(forA(i,j)*w));
            //forAw(i,j)=forAw(i,j);
            }
         }
         A=exp(forAw);
         rebuild_X=A*dataP.row(k).st();
         rebuild.row(k)=rebuild_X.col(0).st();  //???��������ĸ�ֵ�������ģ�����
        }

    for(i=0;i<NX;i++)
        {
        datarebuild.col(i)=ifft(rebuild.col(i),NF);
        }      
    return datarebuild;      
}

cx_fmat linerandomtransmat(int NX, int Np, float dx, float dp, float p1, float w)
{
    int i,j;

    fmat X(NX,1);
    fmat P(1,Np);
    fmat forA(NX,Np);

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
    cx_fmat forAw(NX,Np);
    cx_fmat A(NX,Np);
    for(i=0;i<NX;i++)
        {for(j=0;j<Np;j++)
            {
                (forAw(i,j).real(0.0));
                (forAw(i,j).imag(0.0));
                (forAw(i,j).imag(forA(i,j)*w));
                //forAw(i,j)=forAw(i,j);
            }
        }
    A=exp(forAw);

   return A;
}


#endif
