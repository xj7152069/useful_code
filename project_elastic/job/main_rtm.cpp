#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;

#include "../include/xjc.h"
#include "../include/armadillo"
using namespace arma;

int main()
{
    int nz(300),nx(497),T(5000),i,j,k,sx(50),sx1(380),sx2(480);
    float dz(6),dx(6),dt(0.0005),**movie,**suf,**image,**sr,\
        **smoothmodel,**realmodel;
    char file[99];
    //bmodel=newfmat(n1,nx);
    movie=newfmat(nz,nx);
    suf=newfmat(T,nx);
    image=newfmat(nz,nx);
    sr=newfmat(nz,nx);
    smoothmodel=newfmat(nz,nx);
    realmodel=newfmat(nz,nx);
    wave2D s(nz,nx);
    s.cleardata();
    s.dt=dt,s.dx=dx,s.dy=dz;

/////////////////////////////////////////////////////////
    float **basemodel;
    int n1(750);
    basemodel=newfmat(n1,nx);
    dataread(basemodel,n1,nx,"../data/model.marmousi.750_497.bin");
    matcopy(realmodel,basemodel[0][0],nz,nx);
    for(i=50;i<nz;i++)
    {
        matcopy(realmodel[i],basemodel[int((i-50)*3)],nx);
    }
    //dataread(realmodel,nz,nx,"../data/model.marmousi.750_497.bin");
    datawrite(realmodel,nz,nx,"../data/model1.dat");
    matsmooth(smoothmodel,realmodel,s.ny,s.nx,25);
    datawrite(smoothmodel,nz,nx,"../data/model1.smooth.bin");
////////////////////////////////////////////////////////
    ofstream out1;
    ifstream inf1;

for(sx=sx1;sx<sx2;sx=sx+5)
{
    cout<<"mow is running: "<<sx<<endl;
    
    file[0]='\0';
    strcat(file,"../data/movie1.bin");
    //strcat(file,numtostr(sx,5));
    //out1.open(file);
    s.cleardata();
    matcopy(s.model,realmodel,nz,nx);
    matcopy(image,0.0,nz,nx);

    for(k=0;k<T;k++)
    {
        //if(k%10==0)
        {
        s.s2[int(s.PML_wide+10)][int(sx)]+=10*wavelet01(k,s.dt,30.0);
        }
        s.timeslicecal();
        matcopy(suf[k],s.s2[int(s.PML_wide+10)],nx);
        
        //if(k%5==0)
        {
            
            //matcopy(movie,s.model,nz,nx);
            //matmul(movie,0.0001,nz,nx);
            //matadd(movie,s.s2,nz,nx);
            //datawrite(movie,s.ny,s.nx,out1);
            //datawrite(s.s2,nz,nx,out1);
        }
    }
    //out1.close();
    file[0]='\0';
    strcat(file,"../data/suf/suf1.bin");
    strcat(file,numtostr(sx,5));
    datawrite(suf,T,nx,file);
    /*
    s.cleardata();
    matcopy(s.model,smoothmodel,nz,nx);
    inf1.open(file);
    inf1.seekg(-nz*nx*sizeof(float),ios::end);

    for(k=0;k<T-1;k++)
    {
        
        dataread(sr,nz,nx,inf1);
        inf1.seekg(-2*nz*nx*sizeof(float),ios::cur);
        {
        matadd(s.s2[int(s.PML_wide+10)],suf[T-1-k],nx);
        }
        s.timeslicecal();
        matmul(sr,s.s2,nz,nx);
        matadd(image,sr,nz,nx);
    }
    inf1.close();
    file[0]='\0';
    strcat(file,"../data/image1.bin");
    strcat(file,numtostr(sx,5));
    datawrite(image,nz,nx,file);
    */
}

///////////////////////////////////////////////////////////////////
/*
matcopy(image,0.0,nz,nx);
for(sx=sx1;sx<sx2;sx=sx+5)
{
    file[0]='\0';
    strcat(file,"../data/image1.bin");
    strcat(file,numtostr(sx,5));
    dataread(sr,nz,nx,file);
    matadd(image,sr,nz,nx);
    cout<<file<<" || "<<sx<<endl;
}
Laplace(sr,image,nz,nx);
Laplace(image,sr,nz,nx);
datawrite(image,nz,nx,"../data/image1.all.bin");
*/
    return 0;
}





