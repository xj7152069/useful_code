#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include<iomanip>
#include<math.h>
#include <armadillo>
#include "./include/xjc.h"
using namespace std;
using namespace arma;
float zb(int k, float DT)
{
	float f,pi(3.1415926);
    f=(pi)*(pi)*900*(k*DT-0.04)*\
exp((-pi*pi*900*(k*DT-0.04)*(DT*k-0.04)))\
*(3-2*pi*pi*900*(k*DT-0.04)*(DT*k-0.04));
	return f;
}

int main(int argc, char** argv)
{

    int NZ(1024),NX(200);
    ofstream outf;
    int i,j,k,num(9),f0(100);
    float DT(0.001),outdata;

    fmat A(NZ,NX),B(NZ,NX,fill::zeros);
    A.fill(0.0);
/*
    for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                    if(i>=100+0.5*j)
                        {
                        A(i,j)=A(i,j)+2*zb(int(i-100-0.5*j),DT);
                        }
                    }
            }*/
   
    for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=450)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-450),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        
for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=400-1*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-400+1*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=400+1*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-400-1*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);


for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=300+2*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-300-2*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);

                        for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=700-2*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-700+2*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                       
for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=800-3*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-800+3*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=100+3*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-100-3*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
     /*                   
for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=550+1*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-550-1*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=550-1*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-550+1*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        
for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=300-1.25*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-300+1.25*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=300+1.25*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-300-1.25*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);
                        
for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=600+1.5*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-600-1.5*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);

                        for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        if(i>=600-1.5*j)
                        {
                        A(i,j)=A(i,j)+wavelet01(int(i-600+1.5*j),DT,f0);
                        }}}
                        A=fmatsmooth(A,NZ,NX,num);
                        A=A/(max(max(A.col(NX/2)))-min(min(A.col(NX/2))));
                        B+=A;A.fill(0);

    for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                    if(i>=50)
                        {
                        A(i,j)=A(i,j)+2*zb(int(i-50),DT);
                        }
                    }
            }

     for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                    if(i>=500)
                        {
                        A(i,j)=A(i,j)+2*zb(int(i-500),DT);
                        }
                    }
            }

      for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                    if(i>=sqrt(100*100+10*(j-100)*100*100*(j-100)/3000/3000/DT))
                        {
                        A(i,j)=A(i,j)+zb(int(i-sqrt(100*100+10*(j-100)*100*100*(j-100)/3000/3000/DT)),DT);
                        }
                    }
            }

       for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                    if(i>=sqrt(200*200+8*(j-100)*100*100*(j-100)/3000/3000/DT))
                        {
                        A(i,j)=A(i,j)+zb(int(i-sqrt(200*200+8*(j-100)*100*100*(j-100)/3000/3000/DT)),DT);
                        }
                    }
            }*/ 
A=B;
    for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                        /*
                            if(i>=50)
                            A(i,j)+=B(i-50,j);
                            if(i<(NZ-50))
                            A(i,j)+=B(i+50,j);
                        */
                    }
            }
    outf.open("txz9.bin");
    for(j=0;j<NX;j++)
            {for(i=0;i<NZ;i++)
                    {
                    outdata=A(i,j);
                    outf.write((char*)&outdata,sizeof(outdata));
                    }
            }
    outf.close();

return 0;
}

















































