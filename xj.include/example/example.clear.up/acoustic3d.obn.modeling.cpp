/*
    Simulate 3D OBN seismic data by velocity-stress equation.
    Clear up in 2022.11.05, by Xiang Jian, WPI.

*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "./include/xjc.h"

using namespace arma;
using namespace std;
//argv[1]: model-vp path&filename
//argv[2]: model-vs path&filename
//argv[3]: model-rho path&filename
//argv[4]: surface type;free surface,suface=0.0;
//         PML surface, suface=1.0
//run example:
int main(int argc , char *argv[])
{
    int i,j,k;
    //source X-coord in model
    //sx1: begin X-coord;
    //sx2: end X-coord (Take less);
    //dsx: source X-coord gep;
    //ncpu: Number of threads;
    int sx1(50),sx2(251),dsx(1),ncpu(1);
    //modeling par:
    int nz(200),nx(201),ny(201),nt(2000);
    float dx(5.0),dy(5.0),dz(5.0),dt(0.0003),f0(30.0);
    //source and obn Z-coord in model
    int freeSufaceZ(50),sz(50),rz(85),sx(80),sy(80);

    char filemodel[399];

    elastic3D_ARMA obj(nx,ny,nz);
    obj.dz=dz,obj.dx=dx,obj.dy=dy,obj.dt=dt;
    obj.isPMLSurface=0.0;
    obj.nzSampleOfFreeSurface=freeSufaceZ;
    obj.mpar_ro.fill(1000.0);
    obj.mpar_vp.fill(1500.0);
    obj.mpar_vs.fill(0.0);
    
    for(k=rz;k<nz;k++){
        obj.mpar_vp.slice(k).fill(1500.0);
    }
    for(k=rz+50;k<nz;k++){
        obj.mpar_vp.slice(k).fill(2000.0);
    }
    obj.cleardata();
    obj.maxThreadNum=6;

    fcube data3dtzz(nx,ny,nt);
    if(obj.isPMLSurface<0.5){
    for(k=0;k<obj.nzSampleOfFreeSurface;k++){
        obj.mpar_ro.slice(k).fill(1000.0);
        obj.mpar_vp.slice(k).fill(0.0);
        obj.mpar_vs.slice(k).fill(0.0);
    }}
    obj.updatepar();

    ofstream outf1,outf2;
    outf1.open("movie.p.colsy.dat");
    outf2.open("movie.p.rowsx.dat");
    fmat datacol(nx,nz),datacolt(nz,nx);
    fmat datarow(ny,nz),datarowt(nz,ny);
    for(k=0;k<nt;k++){
        obj.acoustic_p(sx,sy,sz)+=1.0*wavelet01(k,obj.dt,f0,1);

        //Using internal parallelism
        TimeSliceCal_acoustic3D_MultiThread(obj);
        //TimeSliceCal_acoustic2D_MultiThread(obj);

        data3dtzz.slice(k)=obj.acoustic_p.slice(rz);

        if(k%50==0){
            datacol=obj.acoustic_p.col(sy);
            datarow=obj.acoustic_p.row(sx);
            datawrite(datacolt=datacol.st(),nz,nx,outf1);
            datawrite(datarowt=datarow.st(),nz,ny,outf2);
        }

        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    outf2.close();
    datawrite3d_bycol_transpose(data3dtzz,nt,nx,"data3d.p.free.dat");

    return 0;
}


