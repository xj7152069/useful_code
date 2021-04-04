#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

#include "./include/xjc.h"

int main()
{
    int i,j,k,nz(300),nx(497),nt(5000);
    float dz(6),dx(6),dt(0.0005);
    int sx1(20),sx2(480),sx,sxgep(5);
    int rx1(2),rx2(nx-1),rx,rxgep(1);

    wave2D s;
    caltime_par par;
    par.nvz=nz,par.nvx=nx;
    par.sz=int(s.PML_wide+10);
    par.dvz=dz;par.dvx=dx;
    par.vel=newfmat(par.nvz,par.nvx);
    dataread(par.vel,par.nvz,par.nvx,"../data/model.300_497.bin");
    matsmooth(par.vel,par.vel,nz,nx,25);
    datawrite(par.vel,nz,nx,"../data/model.smooth.bin");

    float ***alltime,***allsuf;
    int numsx=(sx2-sx1)/sxgep;
    alltime=newfmat(nx,nz,nx);
    allsuf=newfmat(numsx,nt,nx);
    char file[99];
    ifstream inf1;
    ofstream out1;

    for(sx=sx1;sx<sx2;sx=sx+sxgep)
    {   
        k=(sx-sx1)/sxgep;
        file[0]='\0';
        //偏移读取的地表数据由波动方程模拟得到
        strcat(file,"../data/suf/suf1.bin");
        strcat(file, numtostr(sx,5));
        dataread(allsuf[k],nt,nx,file);
    }

/////////////////////////////////////////////
/*调用Fortran程序，先计算出所有需要的旅行时场，
存入文件之后读取*/
    out1.open("../data/alltime.bin");
    for(rx=rx1;rx<rx2;rx=rx+1)
    {   
        cout<<"now is caltime: "<<rx<<endl;
        par.sx=rx;

        //封装了Fortran子程序的函数
        //par是caltime_par结构体，用于传参
        //位于./include/xj.ray.f90/caltime.hpp
        caltraveltime2d(par); 

        matcopy(alltime[rx],par.time,nz,nx);
        datawrite(par.time,nz,nx,out1);
    }
    out1.close();
/////////////////////////////////////////////
    inf1.open("../data/alltime.bin");
    for(rx=rx1;rx<rx2;rx=rx+1)
    {   
        //cout<<"now is readtime: "<<rx<<endl;
        dataread(alltime[rx],nz,nx,inf1);
    }
    inf1.close();

    float **image,t,v(3000),pi(3.1415926),v1(1500),v2(4000);
    float xs,r,sz(par.sz*dz),du,u,delay_t(0.05);
    int ut;
    image=newfmat(nz,nx);
    matcopy(image,0.0,nz,nx);
    
//柯西霍夫积分循环
    for(sx=sx1;sx<sx2;sx=sx+sxgep)
    {  
        cout<<"now is run: "<<sx<<endl;
        k=(sx-sx1)/sxgep;
        for(i=par.sz;i<nz;i++)
        {
            for(j=0;j<nx;j++)
            {
                for(rx=rx1;rx<rx2;rx=rx+1)
                {  
                    v=(i-par.sz)*(v2-v1)/(nz-par.sz)+v1;
                    r=sqrt(pow((i-par.sz)*dz,2)+pow((j-sx)*dx,2))\
                        +sqrt(pow((i-par.sz)*dz,2)+pow((j-rx)*dx,2));
                    xs=-(i-par.sz)*dz/r/r/v/2/pi;
                    t=alltime[sx][i][j]+alltime[rx][i][j]+delay_t;
                    ut=round(t/dt);
                    if(ut>1 && ut<nt-2)
                    {
                        du=(allsuf[k][ut+1][rx]-allsuf[k][ut-1][rx])/dt/2;
                        u=allsuf[k][ut][rx]*v/r;
                        //du=0;
                        image[i][j]+=xs*(du+u);
                        //image[i][j]+=allsuf[k][ut][rx];
                    }
                }
            }
        }
    }
    
    datawrite(image,nz,nx,"../data/image.kirchhoffpsdm.bin");
    Laplace(image,image,nz,nx);
    datawrite(image,nz,nx,"../data/image.kirchhoffpsdm.la1.bin");
    Laplace(image,image,nz,nx);
    datawrite(image,nz,nx,"../data/image.kirchhoffpsdm.la2.bin");

    return 0;
}
/*
struct caltime_par
{
    int NS_Z,NS_X,nvz,nvx;
    float DVZ,DVX;
    float **vel,**time;
};*/



