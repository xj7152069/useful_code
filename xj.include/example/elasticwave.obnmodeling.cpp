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
void obn2d_modeling(elastic2D* A, int* par1, int* par2, int* par3,\
 int parNt, int parSourceDeepth, int parOBNDeepth, int freeSufaceZ,\
 float parWaveletFrequence);
int main(int argc , char *argv[])
{
    int i,j,k,*par1,*par2;
    //source X-coord in model
    //sx1: begin X-coord;
    //sx2: end X-coord (Take less);
    //dsx: source X-coord gep;
    //ncpu: Number of threads;
    int sx1(500),sx2(501),dsx(1),ncpu(1);
    //modeling par:
    int Z(200),X(1001),T(8000),nz(Z),nx(X),nt(T);
    float dx(5),dz(5),dt(0.0003),f0(30.0);
    //source and obn Z-coord in model
    int freeSufaceZ(60),sz(60),rz(95);

    fmat modelvp(nz,nx),modelvs(nz,nx),modelrho(nz,nx);
    char filemodel[399];
    filemodel[0]='\0';
    //model path&filename
    dataread(modelvp,nz,nx,argv[1]);
    dataread(modelvs,nz,nx,argv[2]);
    dataread(modelrho,nz,nx,argv[3]);

    //elastic-wave-2d model class
    elastic2D* A;
    A=new elastic2D[ncpu];
    for(k=0;k<ncpu;k++){
        //initialize, do first
        A[k].initialize(nz,nx);
        A[k].dt=dt,A[k].dx=dx,A[k].dy=dz;
        A[k].suface=atof(argv[4]);
        matcopy(A[k].vp,modelvp,nz,nx);
        //density-model
        matcopy(A[k].ro,modelrho,nz,nx);
        //model-vs
        matcopy(A[k].vs,modelvs,nz,nx);
        for(i=0;i<sz;i++){
        for(j=0;j<nx;j++){
            if(A[k].suface<=0.5){
                //free surface, air density, suface=0.0
                A[k].ro[i][j]=1.0;
            }
            else{
                //PML surface, water density, suface=1.0
                A[k].vp[i][j]=A[k].vp[sz][j];
                A[k].ro[i][j]=A[k].ro[sz][j];
            }
        }}
        //update class, do before modeling
        A[k].updatepar();
    }
/////////////////////////////////////
    par1=new int[ncpu];
    par2=new int[ncpu];
    int s1(sx1/dsx),s2(sx2/dsx);
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=s1;i<(s2-ncpu);i+=ncpu){
        for(k=0;k<ncpu;k++){
            par1[k]=(i+k)*dsx;
            par2[k]=(i+k+1)*dsx;
            pcal[k]=thread(obn2d_modeling,&A[k],&par1[k],&par2[k],&dsx,T,sz,rz,freeSufaceZ,f0);
        }
        for(k=0;k<ncpu;k++){
            pcal[k].join();
        }
    }
    for(i=max((s2-ncpu),s1);i<(s2);i++){
        k=i-(s2-ncpu);
        par1[k]=(i)*dsx;
        par2[k]=(i+1)*dsx;
        pcal[k]=thread(obn2d_modeling,&A[k],&par1[k],&par2[k],&dsx,T,sz,rz,freeSufaceZ,f0);
    }
    for(i=max((s2-ncpu),s1);i<(s2);i++){
        k=i-(s2-ncpu);
        pcal[k].join();
    }
    return 0;
}

void obn2d_modeling(elastic2D* A, int* par1, int* par2, int* par3,\
 int parNt, int parSourceDeepth, int parOBNDeepth, int freeSufaceZ,\
 float parWaveletFrequence)
{
    char filep[99];
    char files[99];
    int Z(A[0].ny),X(A[0].nx),T(parNt),nz(Z),nx(X),nt(T),i,j,k;
    float f0(parWaveletFrequence);
    int sx(0),sz(parSourceDeepth),rz(parOBNDeepth);
////////////////////////////////////////////////
    int *obnz;
    obnz=new int[nx];
    for(i=0;i<nx;i++){
        //The sea floor is flat
        obnz[i]=rz;
        //The sea floor is rugged
        //obnz[i]=40+20/(1+exp(0.1*(nx/2-i)))-0.04*i\
        +10/(1+exp(0.1*(3*nx/4-i)))+10/(1+exp(0.1*(nx/4-i)));
    }
////////////////////////////////////////////////
    ofstream outf1,outf2,outf3,outf4;
    fmat suf(T,X);
    fmat sufp(T,X);
    fmat sufs(T,X);
    fmat uu(nz,nx),uuu(nz,nx);

    outf1.open("./data/movie.dat");
    A->Zsiteofseasuface=freeSufaceZ;
    for(sx=par1[0];sx<par2[0];sx=sx+par3[0]){
        A[0].cleardata();
        for(k=0;k<T;k++)
        {
            A->Txx[sz][sx]+=1000000*wavelet01(k,A[0].dt,f0,1);
            A->Tyy[sz][sx]+=1000000*wavelet01(k,A[0].dt,f0,1);

            //Using internal parallelism
            timeslicecal_T_thread(A[0]);
            //don't use internal parallelism
            //A->timeslicecal_T();

            for(i=0;i<nx;i++){
                sufp(k,i)=A->Tyy[obnz[i]][i];
                sufs(k,i)=A->uz[obnz[i]][i];
                suf(k,i)=A->ux[obnz[i]][i];
            }

    
            if(k%20==0)
            {
                uu=matcopy(A->ro,Z,X);
                uu*=0.0000001;
                uuu=matcopy(A->uz,Z,X);
                datawrite(uu=uu+uuu, Z, X, outf1);
            }
    
            if(k%1000==0)
            {cout<<"now is running : "<<k<<endl;}
        }
        outf1.close();

            filep[0]='\0';
            files[0]='\0';
            strcat(filep,"./data/obndata.p.cs.bin");
            strcat(filep,numtostr(sx,5));
            strcat(files,"./data/obndata.vz.cs.bin");
            strcat(files,numtostr(sx,5));
            datawrite(sufp,T,X,filep);
            datawrite(sufs,T,X,files);

            filep[0]='\0';
            strcat(filep,"./data/obndata.vx.cs.bin");
            strcat(filep,numtostr(sx,5));
            datawrite(suf,T,X,filep);

        cout<<"Have output :"<<sx<<endl;
    }

}





