
#define ARMA_DONT_USE_OPENMP
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
#include <mpi.h>
/* input par:
Tips: input data format is .su! 
    output data format is .su! 
argv[1]: filein data; input Data-upgoing path&filename;
argv[2]: fileout demultiple; output Data-demultiple path&filename;
argv[3]: nx0; space sampling number of input-data-inline;
argv[4]: ny0; space sampling number of input-data-crossline;
argv[5]: nt0; time sampling number of input-data;
argv[6]: space sampling gep of input-data-inline;
argv[7]: space sampling gep of input-data-crossline;
argv[8]: time sampling gep;
argv[9]: wiener filter Length;
argv[10]: wiener filter Slide Length;
argv[11]: wiener Tikhonov;
argv[12]: number of using thread;
*/
int main(int argc, char * argv[])
{
    int mpiMyProcId,mpiNumProc,\
     dataGroupBegIndex(240),dataGroupNum(4),\
     dataGroupEndIndex(dataGroupNum+dataGroupBegIndex);

//Note: A warning may appears when using MPI;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&mpiMyProcId);
MPI_Comm_size(MPI_COMM_WORLD,&mpiNumProc);
    cout<<"step - <defined variable> is running."<<endl;
    float pi(3.1415926);
    char fileindata[399],fileoutdemultiple[399],nameMpiId[99];
    nameMpiId[0]='\0';
    strcat(nameMpiId,numtostr(mpiMyProcId,5));
    fileindata[0]='\0';
    strcat(fileindata,argv[1]);

    int nx0(atoi(argv[3])),ny0(atoi(argv[4])),nt0(atoi(argv[5])),\
        i,j,k,ncpu(atoi(argv[12])),wienerLength(atoi(argv[9])),\
        wienerTimeSlide(atoi(argv[10]));
    float dx(atof(argv[6])),dy(atof(argv[7])),dt(atof(argv[8])),\
        waterVelocity(1530.0),wienerTikhonov(atof(argv[11])),\
        cutSlopeMax(1450.0), cutSlopeMin(1050.0),fmax(150.0),fmin(0.2);
    int sxIndex,syIndex,nx(nx0),ny(ny0),nt(nt0);
    float sxCoordmin(1e30),syCoordmin(1e30),\
        sxCoordmax(-sxCoordmin),syCoordmax(-syCoordmin);
    bool sufile(true);
    int interTimesByRow(2),interTimesByCol(1);

    ifstream inf1,inf2;
    ofstream outf1,outf2;
    segyhead2 **suHeadArray2dOrig;
    newArray(suHeadArray2dOrig,nx0,ny0);

    segyhead head;
    head.filename[0]='\0';
    strcat(head.filename,fileindata);
    segyhead_open(head,sufile);
    head.endian='l';
    head.isibm=false;

    int begIndex=dataGroupBegIndex+mpiMyProcId;
    for(i=0;i<(dataGroupBegIndex+mpiMyProcId)*nx0*ny0;i++){
        head.infile.seekg(240,ios::cur);
        head.infile.seekg(nt0*sizeof(float),ios::cur);
    }

/************** Data Group processing Loop **************/
while(!head.infile.eof() && begIndex<dataGroupEndIndex){
    fileoutdemultiple[0]='\0';
    strcat(fileoutdemultiple,argv[2]);
    strcat(fileoutdemultiple,numtostr(begIndex,5));
    outf1.open(fileoutdemultiple);
    cout<<"step - <MPI> is running: "<<begIndex<<endl;

/////////////Read su-data and sort by sx-sy/////////////
    segyhead2 **suHeadArray2d;
    fcube data3dOrig,data3dCut,data3dMultiple,data3dTaop,\
        data3dRecover,data3dErr,data3dTx1,data3dTx2;
    cx_fcube data3dFx1,data3dFx2;
    fmat coordx(nx0,ny0,fill::zeros),coordy(nx0,ny0,fill::zeros),\
        coordfold(nx0,ny0,fill::zeros),seabaseDepth(nx0,ny0),\
        coordxOrig,coordyOrig,coordfoldOrig,seabaseDepthOrig;

    cout<<"step - <Read su-data> is running."<<endl;
    data3dOrig.zeros(nx0,ny0,nt0);
    coordfold=readSuDataSortByCoord(suHeadArray2dOrig,data3dOrig,\
        nx0,ny0,nt,dx,dy,head);
    for(i=0;i<nx0;i++){
    for(j=0;j<ny0;j++){
        coordx(i,j)=getSuHeadKey(suHeadArray2dOrig[i][j],"gx");
        coordy(i,j)=getSuHeadKey(suHeadArray2dOrig[i][j],"gy");
        seabaseDepth(i,j)=getSuHeadKey(suHeadArray2dOrig[i][j],"swdep");
    }}
    sxCoordmin=sum(sum(coordx,1)/nx)/ny-dx*(nx-1)/2.0;
    sxCoordmax=sxCoordmin+dx*(nx-1);
    syCoordmin=sum(sum(coordy,1)/nx)/ny-dy*(ny-1)/2.0;
    syCoordmax=syCoordmin+dy*(ny-1);

    float dtHead=getSuHeadKey(suHeadArray2dOrig[0][0],"dt")/1000000.0;
    if(abs(dtHead-dt)>0.0001){
        cout<<"Warning: Time-Sample Gep in Head Not Equal the Input Information?"<<endl;
        cout<<"Warning: The Time-Unit is sec. Please Check it!"<<endl;
        cout<<"Warning: Time-Sample Gep Residual is: "<<abs(dtHead-dt)<<endl;
    }
    float sxCoord=getSuHeadKey(suHeadArray2dOrig[0][0],"sx");
    float syCoord=getSuHeadKey(suHeadArray2dOrig[0][0],"sy");
    coordx=coordx-sxCoord;
    coordy=coordy-syCoord;
    coordfold.fill(1.0);
    coordxOrig=coordx;
    coordyOrig=coordy;
    seabaseDepthOrig=seabaseDepth;
    if(ny>3){seabaseDepth=fmatsmooth(seabaseDepth,nx,ny,1);}
    newArray(suHeadArray2d,nx0,ny0);
    for(i=0;i<nx0;i++){
    for(j=0;j<ny0;j++){
        suHeadArray2d[i][j]=suHeadArray2dOrig[i][j];
    }}

//////////////////Cut data by offset-slope/////////////////
    cout<<"step - <Cut data by offset-slope> is running."<<endl;
    data3dCut=data3dOrig;
    coordfold=cutDataBySlope(data3dCut,coordxOrig,coordyOrig,\
        dt, cutSlopeMin, cutSlopeMax);
    data3dOrig=data3dOrig-data3dCut;
    coordfoldOrig=coordfold;

//////////////Water coding predicts multiple waves///////////
if(1){
    cout<<"step - <Predicts multiple waves> is running."<<endl;
    float df=1.0/dt/nt;
    data3dFx1.zeros(nx,ny,nt);
    cout<<"  sub-step - <data tx2fx> is running."<<endl;
    tx2fx_3d_thread(data3dFx1,data3dCut, ncpu);
    cout<<"  sub-step - <Linear-interpolation in fx-dimension> is running."<<endl;
    for(k=0;k<interTimesByRow;k++){
        cxfcubeLinearInterpolation3dByRow(data3dFx1,ncpu);
        fmatLinearInterpolation2dByRow(coordfold);
        fmatLinearInterpolation2dByRow(coordx);
        fmatLinearInterpolation2dByRow(coordy);
        fmatLinearInterpolation2dByRow(seabaseDepth);
    }
    for(k=0;k<interTimesByCol;k++){
        cxfcubeLinearInterpolation3dByCol(data3dFx1,ncpu);
        fmatLinearInterpolation2dByCol(coordfold);
        fmatLinearInterpolation2dByCol(coordx);
        fmatLinearInterpolation2dByCol(coordy);
        fmatLinearInterpolation2dByCol(seabaseDepth);
    }
    int foldGep;
    foldGep=round(pow(2.0,interTimesByRow));
    for(k=0;k<coordfold.n_rows;k++){
        if(k%foldGep!=0)
            coordfold.row(k).fill(0.0);
    }
    foldGep=round(pow(2.0,interTimesByCol));
    for(k=0;k<coordfold.n_cols;k++){
        if(k%foldGep!=0)
            coordfold.col(k).fill(0.0);
    }
    data3dFx2.copy_size(data3dFx1);
    data3dFx2.fill(0.0);

    cout<<"  sub-step - <predicts multiple> is running."<<endl;
    getSourceIndes(sxIndex,syIndex,coordx,coordy);
    cout<<"  check for Source-Coord: "<<sxIndex<<"|"<<syIndex<<"|"<<dt<<endl;
    multiple_code3d(data3dFx2, data3dFx1, seabaseDepth,\
        coordx, coordy, coordfold, sxIndex,syIndex,\
        waterVelocity, df, fmin/df, fmax/df, 50, 200, 0.25, ncpu);
    data3dFx1.zeros(1,1,1);
    //Four times Anti-Linear-interpolation in the spatial dimension
    cout<<"  sub-step - <AntiLinearInterpolation> is running."<<endl;
    for(k=0;k<interTimesByRow;k++){
        cxfcubeAntiLinearInterpolation3dByRow(data3dFx2);
    }
    for(k=0;k<interTimesByCol;k++){
        cxfcubeAntiLinearInterpolation3dByCol(data3dFx2);
    }
    cout<<"  sub-step - <data fx2tx> is running."<<endl;
    data3dMultiple.copy_size(data3dFx2);
    fx2tx_3d_thread(data3dMultiple,data3dFx2, ncpu);
    data3dFx2.zeros(1,1,1);
    fmat fold;
    fold=cutDataBySlope(data3dMultiple,coordxOrig,coordyOrig,\
        dt, cutSlopeMin, cutSlopeMax);
    cout<<"  sub-step - <datawrite3d_bycol_transpose> is running."<<endl;
}

////////////////High-slope filter for multiple model////////////////
// This procedure is not performed when processing 3D data,
// because of the large amount of computation

/////////////Adaptive subtraction of multiple waves/////////////////
if(1){
    cout<<"step - <Adaptive subtraction of multiple waves> is running."<<endl;
    data3dRecover.copy_size(data3dCut);
    data3dErr.copy_size(data3dCut);
    int dataLength(wienerLength*10),slideLength(wienerTimeSlide);
    slideLength=min(dataLength-wienerLength*2,slideLength);
    for(int iline=0;iline<data3dMultiple.n_cols;iline++){
        fmat data2d,multiple2d,demultiple2d;
        data2d=data3dCut.col(iline);
        multiple2d=data3dMultiple.col(iline);
        data2d=data2d.st();
        multiple2d=multiple2d.st();
        demultiple2d=AdaptiveRemoveMultiple2d(data2d, multiple2d, dt,\
            wienerLength, wienerTikhonov, dataLength, slideLength, ncpu);
        data3dRecover.col(iline)=demultiple2d.st();
        //data2d-=demultiple2d;
        //data3dErr.col(iline)=data2d.t();
    }
}

//////////////////////////////Output Data/////////////////////////////
    data3dRecover+=data3dOrig;
    for(j=0;j<data3dOrig.n_cols;j++){
    for(i=0;i<data3dOrig.n_rows;i++){
        outf1.write((char *)(&suHeadArray2d[i][j]), sizeof(head.head2));
        for(k=0;k<head.nz;k++){
            head.data(k,0)=data3dRecover(i,j,k);
        }
        datawrite(head.data,head.nz,1,outf1);
    }}
    outf1.close();

////////Seeks the File Pointer location in next loop////////
    delete [] suHeadArray2d;
    for(i=0;i<(mpiNumProc)*nx0*ny0;i++){
        head.infile.seekg(240,ios::cur);
        head.infile.seekg(nt0*sizeof(float),ios::cur);
    }
    begIndex+=mpiNumProc;
    if(ifstreamFloatEndOfFile(head.infile)){break;}
}
    head.infile.close();

MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
    return 0;
}


