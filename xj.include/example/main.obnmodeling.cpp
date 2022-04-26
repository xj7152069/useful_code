#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
//#include "./armadillo"
using namespace arma;
using namespace std;

#include "./include/xjc.h"

void obn2d_modeling(elastic2D* A, int* par1, int* par2, int* par3);
int main()
{
    int i,k,ncpu(1),*par1,*par2,sx1(500),sx2(501),dsx(1);
    int Z(150),X(1001),T(5000),nz(Z),nx(X),nt(T);
    float dx(5),dz(5),dt(0.0003),f0(30);

    elastic2D* A;
    A=new elastic2D[ncpu];
    for(k=0;k<ncpu;k++){
        A[k].initialize(nz,nx);
        A[k].dt=dt,A[k].dx=dx,A[k].dy=dz;
        A[k].suface=0;
        dataread(A[k].ro,nz,nx,"./model.rho.complex2.bin");
        dataread(A[k].vp,nz,nx,"./model.vp.complex2.bin");
        dataread(A[k].vs,nz,nx,"./model.vs.complex2.bin");
        matcopy(A[k].vs,0.0,nz,nx);
        A[k].updatepar();
    }
    
/////////////////////////////////////
    par1=new int[ncpu];
    par2=new int[ncpu];
    int s1(sx1/dsx),s2(sx2/dsx);
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=s1;i<(s2-ncpu);i+=ncpu)
    {
        for(k=0;k<ncpu;k++)
        {
            par1[k]=(i+k)*dsx;
            par2[k]=(i+k+1)*dsx;
            pcal[k]=thread(obn2d_modeling,&A[k],&par1[k],&par2[k],&dsx);
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
        pcal[k]=thread(obn2d_modeling,&A[k],&par1[k],&par2[k],&dsx);
    }
    for(i=(s2-ncpu);i<(s2);i++)
    {
        k=i-(s2-ncpu);
        pcal[k].join();
    }
    return 0;
}

void obn2d_modeling(elastic2D* A, int* par1, int* par2, int* par3)
{
    char filep[99];
    char files[99];
    int Z(150),X(1001),T(5000),nz(Z),nx(X),nt(T),i,j,k;
    float dx(5),dz(5),dt(0.0003),f0(30);
    int sx(500),sz(10+00),vspx(600);
////////////////////////////////////////////////
int *obnz;
obnz=new int[nx];
    for(i=0;i<nx;i++)
    {
        //obnz[i]=sz;
        obnz[i]=30+20/(1+exp(0.015*(nx/2-i)));
    }
////////////////////////////////////////////////
    ofstream outf1,outf2,outf3,outf4;
    fmat suf(T,X);
    fmat sufp(T,X);
    fmat sufs(T,X);
    fmat vsp(T,nz);
    //fmat uu(nz,nx);
    //fmat uuu(nz,nx);

    //outf1.open("./data/movie2.bin");

    A->Zsiteofseasuface=sz;
    //sz=130+20/(1+exp(0.015*(nx/2-sx)));
for(sx=par1[0];sx<par2[0];sx=sx+par3[0]){
    A[0].cleardata();
    for(k=0;k<T;k++)
    {
        A->Txx[sz][sx]+=1000000*wavelet01(k,A[0].dt,f0,1);
        A->Tyy[sz][sx]+=1000000*wavelet01(k,A[0].dt,f0,1);

        timeslicecal_T_thread(A[0]);
        //A->timeslicecal_T();

        //uu=matcopy(A.uz,Z,X);
        for(i=0;i<nx;i++)
        {
            sufp(k,i)=A->Tyy[obnz[i]][i];
            sufs(k,i)=A->uz[obnz[i]][i];
            suf(k,i)=A->ux[obnz[i]][i];
        }
        for(i=0;i<nz;i++){
            vsp(k,i)=A->Tyy[i][vspx];
        }

/*
        if(k%20==0)
        {
            uu=matcopy(A.vp,Z,X);
            uu*=0.0000001;
            uuu=matcopy(A.uz,Z,X);
            datawrite(uu=uu+uuu, Z, X, outf1);
        }
*/
        //if(k%1000==0)
        //{cout<<"now is running : "<<k<<endl;}
    }
    //outf1.close();

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

        datawrite(vsp,T,nz,"vsp.p.bin");

    cout<<"Have output :"<<sx<<endl;
}
}





