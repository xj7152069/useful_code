
//输出地下速度模型

#include <iostream>
#include <stdio.h>
#include <fstream>
#include<iomanip>
#include<cmath>
#include "./include/xjc.h"
using namespace std;
float p[1001][400],p2[1001][400];
int main()
{

float phase,pi=3.1415926,dx(5),dz(5);
int i,j,k,n,i1,j1,nz(150),nx(1001);
fmat m1(nz,nx),m2(nz,nx),m3(nz,nx);
fmat base(nx,1);
for(i=0;i<nx;i++)
{
  phase=(float(i)/nx)*pi*2.0;
  base(i,0)=9*sin(phase);
}
/*
m1=dataread(nz,nx,"./model.300_497.bin");
for(i=0;i<nx;i++)
{
  m2.col(i*2)=m1.col(i);
  m2.col(i*2+1)=m1.col(i);
}
for(i=0;i<10;i++)
{
  m2.row(i).fill(0);
}
for(j=0;j<nx;j++)
{
  for(i=1;i<nz;i++)
  {
    if(m2(i-1,j)<=1500 && m2(i,j)>1500)
    cout<<i<<"|";
  }
}
datawrite(m2,nz,nx*2,"./model.300.994.bin");
*/
/*
for(i=0;i<nx;i++)
  {for(j=0;j<nz;j++)
    {  	
      m1(j,i)=0.001;
        if (j>=10)
            m1(j,i)=1.0;
        
        if (j>=(base(i,0)+60))
        {
            m1(j,i)=2.0;
        if ((j-1)<(base(i,0)+60))
        {
            m1(j-1,i)=1.0+(2.000-1.0)*(1-(-(j-1)+(base(i,0)+60)));

        }
        }

        if (j>=(base(i,0)+80))
        {
            m1(j,i)=2.1;
        if ((j-1)<(base(i,0)+80))
        {
            m1(j-1,i)=2.000+(2.100-2.000)*(1-(-(j-1)+(base(i,0)+80)));

        }
        }

        if (j>=(base(i,0)+99))
        {
            m1(j,i)=2.200;
        if ((j-1)<(base(i,0)+99))
        {
            m1(j-1,i)=2.100+(2.200-2.100)*(1-(-(j-1)+(base(i,0)+99)));

        }
        }

        if (j>=140)
            m1(j,i)=2.300;

        if (j>=180)
            m1(j,i)=2.400;

     }
   }
   */
//m2=fmatsmooth(m1,nz,nx,0);
float seabase;
for(i=0;i<nx;i++)
  {for(j=0;j<nz;j++)
    {  	
      m1(j,i)=1; //ro
      m2(j,i)=0; //vp
      m3(j,i)=0; //vs
        if (j>=10+00) //water
            {m1(j,i)=1000;
            m2(j,i)=1500;
            m3(j,i)=0.000;}

        seabase=30+20/(1+exp(0.015*(nx/2-i)));
        if (j>seabase) //seabase
           {m1(j,i)=1200+0*dz*(j-seabase);
            m2(j,i)=1500+0*dz*(j-seabase);
            m3(j,i)=0+1.5*dz*(j-seabase);}
        if (j<=seabase && j>=(seabase-1)) //seabase
            {m1(j,i)=1000+200*(j-seabase+1);}

        seabase=110-0.04*i;
        if (j>seabase){
            m1(j,i)=1800+1.25*dz*(j-80-100)*0;
            m2(j,i)=1500+1.5*dz*(j-80-100)*0;
            m3(j,i)=1000+1.5*dz*(j-50-100)*0;}
        if (j<=seabase && j>=(seabase-1)) //seabase
            {m1(j,i)=m1(j,i)+(1800-m1(j,i))*(j-seabase+1);
            m2(j,i)=m2(j,i)+(1500-m2(j,i))*(j-seabase+1);
            m3(j,i)=m3(j,i)+(1000-m3(j,i))*(j-seabase+1);}
    }}

datawrite(m1,nz,nx,"model.rho.complex2.bin");
datawrite(m2,nz,nx,"model.vp.complex2.bin");
datawrite(m3,nz,nx,"model.vs.complex2.bin");

cout<<"finished"<<endl;

return (0);
}

