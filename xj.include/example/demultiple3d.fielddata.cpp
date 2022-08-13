#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "../include/xjc.h"
#include <omp.h>
#include <armadillo>
using namespace arma;
using namespace std;

int main(int argc, char * argv[])
{
    cout<<"step - <defined variable> is running."<<endl;
    float pi(3.1415926);
    char fileindata[399],filemultiple[399],fileoutdemultiple[399];
    int nx0(atoi(argv[4])),ny0(atoi(argv[5])),nt0(atoi(argv[6])),\
        nx(nx0),ny(ny0),nt(nt0+1000),i,j,k;
    float dx(50),dy(25),dt,sxCoord,syCoord,waterVelocity(1530),\
        cutSlope(1050),frequence,df,fmax(200),fmin,\
        regularL2(0.1),regularL1(0.0),nIteratie(45),residualRatio(0.01);
    int npx(71),npy(71),sxIndex,syIndex,ncpu(atoi(argv[3]));
    float dpx(0.00001),dpy(0.00001),sxCoordmin(1e30),syCoordmin(1e30),\
        sxCoordmax(-sxCoordmin),syCoordmax(-syCoordmin),offset;
    bool sufile(true);
    fcube data3dOrig,data3dCut,data3dMultiple,data3dTaop,\
        data3dRecover,data3dErr,data3dTx1,data3dTx2;
    cx_fcube data3dFx1,data3dFx2;
    fmat data2dtx1,data2dtp1,data2dtx2,data2dtp2;
    fmat data2dtxTrans1,datatpTrans1,data2dtxTrans2,datatpTrans2;
    fmat coordx(nx,ny,fill::zeros),coordy(nx,ny,fill::zeros),\
        coordfold(nx,ny,fill::zeros),seabaseDepth(nx,ny),\
        coordxOrig,coordyOrig,coordfoldOrig,seabaseDepthOrig;
    ifstream inf1,inf2,inf0;
    ofstream outf1,outf2,outf0;
    segyhead2 suHeadSwap;
    segyhead2 **suHeadArray2d;
    segyhead2 **suHeadArray2dSort;
    suHeadArray2d=new segyhead2*[nx0];
    suHeadArray2dSort=new segyhead2*[nx0];
    for(i=0;i<nx0;i++){
        suHeadArray2d[i]=new segyhead2[ny0];
        suHeadArray2dSort[i]=new segyhead2[ny0];
    }
////////////////Read su-data and sort by sx-sy/////////////////
    cout<<"step - <Read su-data> is running."<<endl;
    fileindata[0]='\0';
    strcat(fileindata,argv[1]);
    fileoutdemultiple[0]='\0';
    strcat(fileoutdemultiple,argv[2]);
    filemultiple[0]='\0';
    strcat(filemultiple,"./multiple.code.dat");
    segyhead head;
    head.filename[0]='\0';
    strcat(head.filename,fileindata);
    segyhead_open(head,sufile);
    head.endian='l';
    head.isibm=false;

    data3dOrig.zeros(nx,ny,nt);
    //data3dTx1.zeros(nx,ny,nt);
    for(j=0;j<ny;j++){
        for(i=0;i<nx;i++){
            head.dataraw=segyhead_readonetrace_tofmat(head,head.data);
            suHeadArray2d[i][j]=head.head2;
            suHeadArray2dSort[i][j]=suHeadArray2d[i][j];
            for(k=0;k<min(nt,head.nz);k++){
                //data3dTx1(i,j,k)=head.data(k,0);
                data3dOrig(i,j,k)=head.data(k,0);
            }
            coordx(i,j)=getSuHeadKey(suHeadArray2d[i][j],"gx");
            coordy(i,j)=getSuHeadKey(suHeadArray2d[i][j],"gy");
            seabaseDepth(i,j)=getSuHeadKey(suHeadArray2d[i][j],"swdep");
        }
    }
    sxCoordmin=sum(sum(coordx,1))/nx/ny-dx*(nx-1)/2.0;
    sxCoordmax=sxCoordmin+dx*(nx-1);
    syCoordmin=sum(sum(coordy,1))/nx/ny-dy*(ny-1)/2.0;
    syCoordmax=syCoordmin+dy*(ny-1);
/*
    int nOverfold=0;
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            int sxi,syj;
            float fsxi,fsyj;
            fsxi=getSuHeadKey(suHeadArray2d[i][j],"gx");
            fsyj=getSuHeadKey(suHeadArray2d[i][j],"gy");
            sxi=round((fsxi-sxCoordmin)/dx);
            syj=round((fsyj-syCoordmin)/dy);
            sxi=min(sxi,nx-1);
            syj=min(syj,ny-1);
            sxi=max(sxi,0);
            syj=max(syj,0);
            if(coordfold(sxi,syj)<0.5){
                coordfold(sxi,syj)+=1;
                coordx(sxi,syj)=fsxi;coordy(sxi,syj)=fsyj;
                suHeadArray2dSort[sxi][syj]=suHeadArray2d[i][j];
                data3dOrig(span(sxi,sxi),span(syj,syj),span(0,nt-1))\
                    =data3dTx1(span(i,i),span(j,j),span(0,nt-1));
            }else{
                float fsxi1,fsyj1;
                float fsxi2,fsyj2;
                float ferr1,ferr2;
                fsxi1=((fsxi-sxCoordmin)/dx);
                fsyj1=((fsyj-syCoordmin)/dy);
                fsxi2=((coordx(sxi,syj)-sxCoordmin)/dx);
                fsyj2=((coordy(sxi,syj)-syCoordmin)/dy);
                ferr1=(fsxi1-sxi)*(fsxi1-sxi)+(fsyj1-syj)*(fsyj1-syj);
                ferr2=(fsxi2-sxi)*(fsxi2-sxi)+(fsyj2-syj)*(fsyj2-syj);
                if(ferr1<ferr2){
                    coordx(sxi,syj)=fsxi;coordy(sxi,syj)=fsyj;
                    suHeadArray2dSort[sxi][syj]=suHeadArray2d[i][j];
                    data3dOrig(span(sxi,sxi),span(syj,syj),span(0,nt-1))\
                        =data3dTx1(span(i,i),span(j,j),span(0,nt-1));
                    cout<<"Warning: Overfold: "<<nOverfold<<endl;
                }
                nOverfold++;
            }
        }
    }
    cout<<"Warning: Number of overfold trace is: "<<nOverfold<<endl;
    nOverfold=0;
*/
    dt=getSuHeadKey(suHeadArray2dSort[0][0],"dt")/1000000;
    sxCoord=getSuHeadKey(suHeadArray2dSort[0][0],"sx");
    sxIndex=round((sxCoord-sxCoordmin)/dx);
    syCoord=getSuHeadKey(suHeadArray2dSort[0][0],"sy");
    syIndex=round((syCoord-syCoordmin)/dy);
    coordx=coordx-sxCoord;
    coordy=coordy-syCoord;
    coordfold.fill(10);
    coordxOrig=coordx;
    coordyOrig=coordy;
    coordfoldOrig=coordfold;
    seabaseDepthOrig=seabaseDepth;
    if(ny>3)seabaseDepth=fmatsmooth(seabaseDepth,nx,ny,5);
//////////////////Cut data by offset-slope/////////////////
    cout<<"step - <Cut data by offset-slope> is running."<<endl;
    //data3dCut.zeros(nx,ny,nt);
    for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
        offset=sqrt(coordx(i,j)*coordx(i,j)\
            +coordy(i,j)*coordy(i,j));
        int nCut(floor(offset/dt/cutSlope)+1);
        for(k=0;k<min(nt,nCut);k++){
            //data3dCut(i,j,k)=data3dOrig(i,j,k);
            data3dOrig(i,j,k)=0;
        }
        if(nCut>=nt0)
            coordfold(i,j)=0;
    }}
    datawrite3d_bycol_transpose(data3dOrig,nt0,nx,"data3d.cut.dat");

//////////////////Water coding predicts multiple waves///////////////
if(1){
    cout<<"step - <Predicts multiple waves> is running."<<endl;
    df=1.0/dt/nt;
    data3dFx1.zeros(nx,ny,nt);
    cout<<"  sub-step - <data tx2fx> is running."<<endl;
    tx2fx_3d_thread(data3dFx1,data3dOrig, ncpu);
    cout<<"  sub-step - <Linear-interpolation in fx-dimension> is running."<<endl;
    cxfcubeLinearInterpolation3dByRow(data3dFx1,coordfold);
    cxfcubeLinearInterpolation3dByRow(data3dFx1,coordfold);
    fmatLinearInterpolation3dByRow(coordx);
    fmatLinearInterpolation3dByRow(coordx);
    fmatLinearInterpolation3dByRow(coordy);
    fmatLinearInterpolation3dByRow(coordy);
    fmatLinearInterpolation3dByRow(seabaseDepth);
    fmatLinearInterpolation3dByRow(seabaseDepth);
    cxfcubeLinearInterpolation3dByCol(data3dFx1,coordfold);
    fmatLinearInterpolation3dByCol(coordx);
    fmatLinearInterpolation3dByCol(coordy);
    fmatLinearInterpolation3dByCol(seabaseDepth);
    nx=coordx.n_rows;
    ny=coordx.n_cols;
    data3dFx2.zeros(nx,ny,nt);
    cout<<"  sub-step - <predicts multiple> is running."<<endl;
    getSourceIndes(sxIndex,syIndex,coordx,coordy);
    cout<<"  check for Source-Coord: "<<sxIndex<<"|"<<syIndex<<"|"<<dt<<endl;
// void multiple_code3d(cx_fcube& u2, cx_fcube& u1, fmat& seabase_depth,\
 fmat& coordx_data, fmat& coordy_data, int system_source_ix, int system_source_jy,\
 float water_velocity,float dx, float dy, float df, int fn1, int fn2,\
 int minspacewin,int maxspacewin,float code_pattern, int ncpu, float wavelet_delay=0.0)
    multiple_code3d(data3dFx2, data3dFx1, seabaseDepth,\
        coordx, coordy, coordfold, sxIndex,syIndex,\
        waterVelocity, df, 1*df, fmax/df, 50,200, 0.25, ncpu);
    data3dFx1.zeros(1,1,1);
    //Four times Anti-Linear-interpolation in the spatial dimension
    cout<<"  sub-step - <AntiLinearInterpolation> is running."<<endl;
    fmat coordfoldcopy;
    coordfoldcopy=coordfold;
    cxfcubeAntiLinearInterpolation3dByRow(data3dFx2,coordfoldcopy);
    cxfcubeAntiLinearInterpolation3dByRow(data3dFx2,coordfoldcopy);
    cxfcubeAntiLinearInterpolation3dByCol(data3dFx2,coordfoldcopy);
    nx=coordfoldcopy.n_rows;
    ny=coordfoldcopy.n_cols;
    cout<<"  sub-step - <data fx2tx> is running."<<endl;
    data3dMultiple.zeros(nx,ny,nt);
    fx2tx_3d_thread(data3dMultiple,data3dFx2, ncpu);
    data3dFx2.zeros(1,1,1);
    cout<<"  sub-step - <datawrite3d_bycol_transpose> is running."<<endl;
    datawrite3d_bycol_transpose(data3dMultiple,nt0,nx,filemultiple);
}
////////////////High-slope filter for multiple model//////////////////
if(1){
    cout<<"step - <High-slope filter for multiple model> is running."<<endl;
    dataread3d_bycol_transpose(data3dMultiple,nt0,nx,filemultiple);
    fmat coordfoldcopy=coordfoldOrig;
    //fcubeLinearInterpolation3dByRow(data3dMultiple,coordfoldcopy);
    //fcubeLinearInterpolation3dByRow(data3dMultiple,coordfoldcopy);
    //fcubeLinearInterpolation3dByCol(data3dMultiple,coordfoldcopy);
    //fcubeLinearInterpolation3dByCol(data3dMultiple,coordfoldcopy);
    nx=coordfoldcopy.n_rows;
    ny=coordfoldcopy.n_cols;
    if(ny<50){
        npy=1;npx*=10;dpx/=10;
    }
    int doinv(1);
    data3dTaop.zeros(npx,npy,nt);
    if(doinv>0){
    k=Beamforming_CG_3D(data3dTaop,data3dMultiple,data3dErr,\
        data3dMultiple,coordxOrig,coordyOrig, nt, nx,ny,dt,\
        npx,-npx*dpx/2, dpx,npy,-npy*dpy/2, dpy,fmax,fmax/2,\
        ncpu, regularL2, regularL1,nIteratie,residualRatio);
    }else{
    k=LinerRadon3d(data3dTaop,data3dMultiple,data3dErr,\
        data3dMultiple,coordxOrig,coordyOrig, nt, nx,ny,dt,\
        npx,-npx*dpx/2, dpx,npy,-npy*dpy/2, dpy,fmax,fmax/2,\
        ncpu, regularL2, regularL1,nIteratie,residualRatio);
    }
    data2dtp1.zeros(npx,nt);
    for(k=0;k<npy;k++){
        data2dtp1=data3dTaop.col(k);
        data2dtp1=get_blackman_downwin2d(data2dtp1,100);
        data2dtp1=get_blackman_upwin2d(data2dtp1,100);
        data3dTaop.col(k)=data2dtp1;
    }
if(npy>20){
    data2dtp1.zeros(npy,nt);
    for(k=0;k<npx;k++){
        data2dtp1=data3dTaop.col(k);
        data2dtp1=get_blackman_downwin2d(data2dtp1,10);
        data2dtp1=get_blackman_upwin2d(data2dtp1,10);
        data3dTaop.col(k)=data2dtp1;
    }
}
    datawrite3d_bycol_transpose(data3dTaop,nt0,npx,"./multiple.tp.bin");

    k= Beamforming_recoverdata_3D(data3dMultiple,data3dTaop,\
        coordxOrig,coordyOrig, nt, nx,ny,dt,\
        npx,-npx*dpx/2, dpx,npy,-npy*dpy/2, dpy,fmax,fmax/2,\
        ncpu);
    //fcubeAntiLinearInterpolation3dByRow(data3dMultiple,coordfoldcopy);
    //fcubeAntiLinearInterpolation3dByRow(data3dMultiple,coordfoldcopy);
    //fcubeAntiLinearInterpolation3dByCol(data3dMultiple,coordfoldcopy);
    //fcubeAntiLinearInterpolation3dByCol(data3dMultiple,coordfoldcopy);
    nx=coordfoldcopy.n_rows;
    ny=coordfoldcopy.n_cols;
    datawrite3d_bycol_transpose(data3dMultiple,nt0,nx,filemultiple);
}
/////////////Adaptive subtraction of multiple waves/////////////////
if(1){
    cout<<"step - <Adaptive subtraction of multiple waves> is running."<<endl;
    data3dMultiple.copy_size(data3dOrig);
    data3dRecover.copy_size(data3dOrig);
    data3dErr.copy_size(data3dOrig);
    cout<<"  sub-step - <dataread3d_bycol_transpose> is running."<<endl;
    dataread3d_bycol_transpose(data3dMultiple,nt0,nx,filemultiple);
    datawrite3d_bycol_transpose(data3dMultiple,nt0,nx,"./multiple.check.dat");
    cout<<"  sub-step - <de-multiple in-line> is running.";
for(int iline=0;iline<ny;iline++){
    demultiple2d par;
    par.dt=dt,par.nt=nt,par.dx=dx,par.nx=nx;
    par.data_remove_wave.zeros(par.nt,par.nx);
    par.data_antiremove_wave.zeros(par.nt,par.nx);
    par.data2d.copy_size(par.data_remove_wave);
    par.datadict.copy_size(par.data_remove_wave);
    par.datad.copy_size(par.data_remove_wave);
    par.datah.copy_size(par.data_remove_wave);
    par.datahd.copy_size(par.data_remove_wave);
    par.nwp=(9),par.nwph=(9),par.nwpd=(9),par.nwphd=(9),\
    par.n1w=(300),par.d1w=(30),par.lmd=(0.001);
    par.datawinbeg.zeros(par.nx,1);
    par.datawinend.copy_size(par.datawinbeg);

    fmat datap,datavz;
    fmat datap_t,datavz_t;
    datap_t.zeros(par.nx,par.nt);
    datavz_t.zeros(par.nx,par.nt);

    //dataread(datavz,filein);
    //dataread(datap,fileincode);
    datap_t=data3dMultiple.col(iline);
    datavz_t=data3dOrig.col(iline);
    datap=datap_t.t();
    datavz=datavz_t.t();
    datap=datap/(datap.max()-datap.min());
    datavz=datavz/(datavz.max()-datavz.min());

/**************8888****************/
omp_set_num_threads(ncpu);
#pragma omp parallel for
    for(k=0;k<par.nx;k++){
        fmat datal(nt,1);
        //cout<<k<<"|";
        par.datah.col(k)=hilbert1D(datal=datap.col(k),nt,dt);
        par.datad.col(k)=diff1D(datal=datap.col(k),nt,dt);
        par.datahd.col(k)=diff1D(datal=par.datah.col(k),nt,dt);
    }
    datawrite(par.datah,"./media/obndata.p2.h.bin");
    datawrite(par.datad,"./media/obndata.p2.d.bin");
    datawrite(par.datahd,"./media/obndata.p2.hd.bin");
    dataread(par.datah,"./media/obndata.p2.h.bin");
    dataread(par.datad,"./media/obndata.p2.d.bin");
    dataread(par.datahd,"./media/obndata.p2.hd.bin");
/**************************************************/
    par.data2d=datavz, par.datadict=datap;
    int nw(par.nwp+par.nwpd+par.nwph+par.nwphd);
    fmat filter(nw,par.data2d.n_cols,fill::zeros);
    AdaptiveRemoveMultiple2d(par,1,ncpu);
    cout<<iline<<"-"<<data3dRecover.n_cols<<"|";
    data3dRecover.col(iline)=par.data_remove_wave.t();
    datavz=par.data_remove_wave-par.data_antiremove_wave;
    data3dErr.col(iline)=datavz.t();
}
    datawrite3d_bycol_transpose(data3dRecover,nt0,nx,"./demultiple.dat");
    datawrite3d_bycol_transpose(data3dErr,nt0,nx,"./multiple.remove.dat");
}
///////////////////////////////////////////////////////////////////
    outf0.open(fileoutdemultiple);
    for(j=0;j<ny;j++){
        for(i=0;i<nx;i++){
            outf0.write((char *)(&suHeadArray2dSort[i][j]), sizeof(head.head2));
            for(k=0;k<head.nz;k++){
                head.data(k,0)=data3dRecover(i,j,k);
            }
            datawrite(head.data,head.nz,1,outf0);
        }
    }
    outf0.close();
    head.infile.close();

    return 0;
}


