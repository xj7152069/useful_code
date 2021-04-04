#ifndef CALTIME_HPP
#define CALTIME_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>


using namespace std;
using namespace arma;

#include "../xjc.h"
extern"C"
{
    void traveltime_2d_(int* ,int* ,int* ,int* ,\
        float* ,float* ,float* ,float* ,float* ,\
        float* ,float* );
}

struct caltime_par
{
    int sz,sx,nvz,nvx;
    float dvz,dvx;
    float **vel,**time;
};

void caltraveltime2d(struct caltime_par & par)
{
    float depth,dzs,dxs,*vel,*time;
    int i,j,n=par.nvz*par.nvx;
    dzs=par.dvz,dxs=par.dvx;
    depth=dzs*par.nvz;
    vel=new float[n];
    time=new float[n];
    par.time=newfmat(par.nvz,par.nvx);
    matcopy(par.time,0.0,par.nvz,par.nvx);

    for(i=0;i<par.nvx;i++)
    {
        for(j=0;j<par.nvz;j++)
        {
            vel[i*par.nvz+j]=par.vel[j][i];
        }
    }
    
    traveltime_2d_(&par.sz,&par.sx,&par.nvz,&par.nvx,\
        &depth,&par.dvz,&par.dvx,&dzs,&dxs,vel,time);
    
    for(i=0;i<par.nvx;i++)
    {
        for(j=0;j<par.nvz;j++)
        {
            par.time[j][i]=time[i*par.nvz+j];
        }
    }
}

//SUBROUTINE TRAVELTIME_2D(NS_Z,NS_X,nvz,nvx,DEPTH,\
                DVZ,DVX,dzs,dxs,VEL,TIME2)		

#endif