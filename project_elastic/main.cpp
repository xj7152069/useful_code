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

    int Z(500),X(500),T(3000);
    elastic2D A(Z,X);
    ofstream outf1,outf2,outf3,outf4;
    float **sufs,**sufp,**w,**uu;
    sufs=newfmatcs(T,X,0.0);
    sufp=newfmatcs(T,X,0.0);
    uu=newfmatcs(Z,X,0.0);
    outf1.open("u1movie.bin");
    //outf2.open("u2movie.bin");
    outf3.open("u3movie.bin");
    //outf4.open("u4movie.bin");
    w=gauss_souce_windows(Z,X,50,X/2,5);
    datawrite(w,Z,X,"sourcewin.bin");
    int k,i,j;
    for(i=0;i<X;i++)
    {
        for(j=0;j<200;j++)
        {
            A.vp[j][i]=2500;
            A.vs[j][i]=1800;
        }
        for(j=200;j<350;j++)
        {
            A.vp[j][i]=3500;
            A.vs[j][i]=2300;
        }
        for(j=350;j<Z;j++)
        {
            A.vp[j][i]=4500;
            A.vs[j][i]=3000;
        }
    }
    A.updatepar();
    datawrite(A.vp,Z,X,"model.vp.bin");

    for(k=0;k<T;k++)
    {
        //addsouce(A.ux,Z/2,X/2,5,30,A.dt,k,1);
        for(i=0;i<A.ny;i++)
        {
            for(j=0;j<A.nx;j++)
            {
            A.ux[i][j]+=w[i][j]*wavelet02(k,A.dt,30,1);
            }
        }

        A.timeslicecal_u();
        matcopy(sufp[k],A.data.vpx[50],X);

        //matcopy(sufp[k],A.vpy[50],X);

        if(k%1==0)
        {
            //matsmooth(uu,A.ux,Z,X,0);
            matcopy(uu,A.vp,Z,X);
            matmul(uu,0.0001,Z,X);
            matadd(uu,A.ux,Z,X);
            datawrite(uu, Z, X, outf1);
            //matsmooth(uu,A.uz,Z,X,0);
            //datawrite(A.uz, Z, X, outf2);
            datawrite(A.data.vpx, Z, X, outf3);
            //datawrite(A.data.vsx, Z, X, outf4);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    //outf2.close();
    outf3.close();
    //outf4.close();
    datawrite(sufp,T,X,"sufp.bin");
    cout<<"Have output test movie."<<endl;

    //elastic_test();

    return 0;
}





