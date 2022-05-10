#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <thread>
#include <future>

using namespace std;
#include "./include/xjc.h"
void obn2d_modeling(int* par1, int* par2, int* par3);
int main()
{
    int i,k,ncpu(25),*par1,*par2;
    par1=new int[ncpu];
    par2=new int[ncpu];
    int sx1(150),sx2(850),dsx(5),s1(sx1/dsx),s2(sx2/dsx);
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=s1;i<=(s2-ncpu);i+=ncpu)
    {
        for(k=0;k<ncpu;k++)
        {
            par1[k]=(i+k)*dsx;
            par2[k]=(i+k+1)*dsx;
            pcal[k]=thread(obn2d_modeling,&par1[k],&par2[k],&dsx);
        }
        for(k=0;k<ncpu;k++)
        {
            pcal[k].join();
        }
    }
    for(i=(s2-ncpu);i<(s2);i++)
    {
        k=i-(s2-ncpu);
        par1[k]=(i)*dsx;
        par2[k]=(i+1)*dsx;
        pcal[k]=thread(obn2d_modeling,&par1[k],&par2[k],&dsx);
    }
    for(i=(s2-ncpu);i<(s2);i++)
    {
        k=i-(s2-ncpu);
        pcal[k].join();
    }
    return 0;
}

void obn2d_modeling(int* par1, int* par2, int* par3)
{
    int nz(200),nx(1001),T(8000),i,j,k,sx(50),\
        sx1(par1[0]),sx2(par2[0]),dsx(par3[0]),sz(10),rz(50);
    float dz(5),dx(5),dt(0.0003),f0(30.0),**movie,**suf,**image,**sr,\
        **smoothmodel,**realmodel;
    char file[99];
    suf=newfmat(T,nx);
    movie=newfmat(nz,nx);

    realmodel=newfmat(nz,nx);
    wave2D s(nz,nx);
    s.cleardata();
    s.dt=dt,s.dx=dx,s.dy=dz;

/////////////////////////////////////////////////////////
    float **basemodel;
    basemodel=newfmat(nz,nx);
    dataread(basemodel,nz,nx,"./model.sigmoid.200.1001.bin");
    matcopy(realmodel,basemodel,nz,nx);
////////////////////////////////////////////////////////
    ofstream out1;
    ifstream inf1;

for(sx=sx1;sx<sx2;sx=sx+dsx)
{
    cout<<"mow is running: "<<sx<<endl;
    file[0]='\0';
    strcat(file,"./data/movie.bin");
    strcat(file,numtostr(sx,5));
    //out1.open(file);

    s.cleardata();
    matcopy(s.model,realmodel,nz,nx);
    s.updatepar();//!!!!
    for(k=0;k<T;k++)
    {
        {
        s.s2[int(sz)][int(sx)]+=100*wavelet01(k,s.dt,f0);
        }
        s.timeslicecal();
        matcopy(suf[k],s.s2[int(rz)],nx);
        if((k%1000)==0)
        {cout<<k<<endl;}
        /*
        if(k%10==0)
        {
            matcopy(movie,s.model,nz,nx);
            matmul(movie,0.0001,nz,nx);
            matadd(movie,s.s2,nz,nx);
            datawrite(movie,s.ny,s.nx,out1);
            //datawrite(s.s2,nz,nx,out1);
        }*/
    }
    //out1.close();
    file[0]='\0';
    strcat(file,"./data/obndata.bin");
    strcat(file,numtostr(sx,5));
    datawrite(suf,T,nx,file);
}

}





