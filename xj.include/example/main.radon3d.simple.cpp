
#define ARMA_DONT_USE_BLAS //!!!!
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "./include/xjc.h"

using namespace std;
using namespace arma;

int main()
{
    int i,j,k,usecopy(0),i1,j1,k1;
    float fmax(90),frule(50),factor_L2(0.01), factor_L1(10),\
        iterations_num(45), residual_ratio(0.01),ncpu(5);

    //int nz2(3501),nz(3501),nx(200),ny(1),nf(3501);
    int nz(1024),nx(200),ny(200),nf(1024),nz2(nz);
    float npx(41),npy(41),dt(0.001),dx(10),dy(10),dpx(0.00001),dpy(0.000005);

    //int nz(1024),nx(200),ny(1),nf(1024),nz2(3501);
    //int nz(1024),nx(200),ny(1),nf(1024),nz2(nz);
    //float npx(501),npy(1),dt(0.001),dx(10),dy(10),dpx(0.000002),dpy(0.000002);

    //int nz(3000),nx(100),ny(1),nf(nz);
    float x0(-dx*nx/2),y0(-dy*ny/2),px0(-dpx*npx/4),py0(-dpy*npy/3);
    fmat coordx(nx,ny),coordy(nx,ny);
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            coordx(i,j)=dx*i+x0;
            coordy(i,j)=dy*j+y0;
        }}

    ofstream outf1,outf2;
    ifstream inf1,inf2;
    //inf1.open("./data/agcdata.3000.974.bin");
    //inf1.open("./data/txz.three.1024.200.bin");
    //inf1.open("./data/txz9.bin");
    //inf1.open("./data/txz.two.1024.200.bin");
    inf1.open("./data/data3d.two.1024.200.200.bin");

    //inf1.open("./data/data3d.3501.200.200.wiener.bin");
    //outf1.open("./result/data3d.rebuild.bin");
    inf1.seekg(nz*nx*00*sizeof(float),ios::beg);

    fcube data3d(nx,ny,nz),datatp3d(nx,ny,nz),\
        datarecover3d(nx,ny,nz),dataerr3d(nx,ny,nz);
    fmat data,datatp,datarecover,dataerr;
    fmat data_t,datatp_t,datarecover_t,dataerr_t;
    data.zeros(nz,nx);
    datatp.zeros(nz,nx);
    datarecover.zeros(nz,nx);
    dataerr.zeros(nz,nx);
    data_t.zeros(nx,nz);
    datatp_t.zeros(nx,nz);
    datarecover_t.zeros(nx,nz);
    dataerr_t.zeros(nx,nz);

    //data=dataread(nz,nx,inf1);
    dataread3d_bycol_transpose(data3d,nz,nx,inf1);
    data_t=data.t();
//int Beamforming_CG_2D(fmat &tauppanel,fmat &recoverdata,fmat &recovererr,\
 fmat trace, fmat coor, int ns, int ntrace, float dt, int npsr,float psrmin, float dpsr,\
 float fmax,float frule,int ncpu, bool regularization,\
 float factor_L2,float factor_L1, int iterations_num, float residual_ratio);
    //k=Beamforming_CG_2D(datatp_t,datarecover_t,dataerr_t,\
        data_t, coordx, nz, nx,dt, npx, px0, dpx,fmax, frule,\
        ncpu, factor_L2, factor_L1, iterations_num, residual_ratio);
    k=Beamforming_CG_3D(datatp3d,datarecover3d,dataerr3d,\
        data3d,coordx,coordy, nz, nx,ny,dt,\
        npx,px0, dpx,npy,py0, dpy,fmax,frule,\
        ncpu, factor_L2, factor_L1, iterations_num, residual_ratio);
    
    datarecover= datarecover_t.t();
    dataerr= dataerr_t.t();
    datatp= datatp_t.t();
    //datawrite(datarecover,nz,nx,"datarecover.bin");
    //datawrite(dataerr,nz,nx,"dataerr.bin");
    //datawrite(datatp,nz,npx,"datatp.bin");
    datawrite3d_bycol_transpose(dataerr3d,nz,nx,"dataerr.bin");
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,"datarecover.bin");
    datawrite3d_bycol_transpose(datatp3d,nz,npx,"datatp.bin");

    return 0;
}



