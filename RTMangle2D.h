/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

//#include "wave2D.h"
extern float** newfmat(int x1, int x2);
extern float*** newfmat(int x1, int x2, int x3);
template<typename T1> extern void matdelete(T1 **mat, int x1);
template<typename T1> extern void matdelete(T1 ***mat, int x1, int x2);
template <typename T1, typename T2> extern void matcopy(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> extern void matcopy(T1 ***mat1, T2 n, int n1, int n2, int n3);
template <typename T1, typename T2> extern void matcopy(T1 **mat1, T2 **mat2, int nz, int nx);

class OF_2D
{
public:
    float **vx1=NULL, **vx2=NULL, **vz1=NULL, **vz2=NULL, ***wave=NULL;
    float **v_x=NULL, **v_z=NULL; 
    float a,dx,dz,dt;
    float xs[5]={0.083333333,-0.666666667,0.0,0.666666667,-0.083333333};
    int x1,x2,x3;

    OF_2D()
    {
        x1=0;x2=0;x3=0;
        a=1.0,dx=5.0,dz=5.0,dt=0.0005;
        cout<<"Warning: Creat an Empty Object-optical-flow-2D"<<endl;
    }
    OF_2D(int nt, int nz, int nx)
    {
        a=1.0,dx=5.0,dz=5.0,dt=0.0005;
        x1=nt,x2=nz,x3=nx;
        vx1=newfmat(x2,x3);
        vx2=newfmat(x2,x3);
        vz1=newfmat(x2,x3);
        vz2=newfmat(x2,x3);
        v_x=newfmat(x2,x3);
        v_z=newfmat(x2,x3);
        wave=newfmat(x1,x2,x3);
        matcopy(vx1,0.0,x2,x3);
        matcopy(vx1,0.0,x2,x3);
        matcopy(vz1,0.0,x2,x3);
        matcopy(vz1,0.0,x2,x3);
        matcopy(wave,0.0,x1,x2,x3);
    }
    ~OF_2D()
    {
        matdelete(vx1,x2);
        matdelete(vx2,x2);
        matdelete(vz1,x2);
        matdelete(vz2,x2);
        matdelete(v_x,x2);
        matdelete(v_z,x2);
        matdelete(wave,x1,x2);
        vx1=NULL;vx2=NULL;vz1=NULL;vz2=NULL;wave=NULL;
    }

    template<typename T1>
    void addtimeslicecal(T1 **p)
    {
        int k;
        for(k=0;k<x1-1;k++)
        {
            matcopy(wave[k],wave[k+1],x2,x3);
        }
        matcopy(wave[x1-1],p,x2,x3);
    }

    void velocityAverage(float **mat1, float **mat2)
    {
        int i,j;
        for(i=1;i<x2-1;i++)
        {
            for(j=1;j<x3-1;j++)
            {
                mat1[i][j]=mat2[i-1][j-1]*1.0/12+mat2[i-1][j]*1.0/6+mat2[i-1][j+1]*1.0/12\
                    +mat2[i][j-1]*1.0/6+mat2[i][j+1]*1.0/6+mat2[i+1][j-1]*1.0/12\
                    +mat2[i+1][j]*1.0/6+mat2[i+1][j+1]*1.0/12-mat2[i][j];
            }
        }
    }

    void velocityIteration()
    {
        this->velocityAverage(v_x,vx1);
        this->velocityAverage(v_z,vz1);
        int i,j;
        float px,pz,pt;
        for(i=2;i<x2-2;i++)
        {
            for(j=2;j<x3-2;j++)
            {
                pt=wave[x1-1][i][j]*0.5/dt-wave[x1-3][i][j]*0.5/dt;
                px=(xs[0]*wave[x1-2][i][j-2]+xs[1]*wave[x1-2][i][j-1]+\
                    xs[3]*wave[x1-2][i][j+1]+xs[4]*wave[x1-2][i][j+2])/dx;
                pz=(xs[0]*wave[x1-2][i-2][j]+xs[1]*wave[x1-2][i-1][j]+\
                    xs[3]*wave[x1-2][i+1][j]+xs[4]*wave[x1-2][i+2][j])/dz;

                vx2[i][j]=v_x[i][j]-px*(px*v_x[i][j]+pz*v_z[i][j]+pt)/\
                    (a+px*px+pz*pz);
                vz2[i][j]=v_z[i][j]-pz*(px*v_x[i][j]+pz*v_z[i][j]+pt)/\
                    (a+px*px+pz*pz);
            }
        }
        matcopy(vx1,vx2,x2,x3);
        matcopy(vz1,vz2,x2,x3);
    }

    void velocityInitialization()
    {
        matcopy(vx1,0.0,x2,x3);
        matcopy(vz1,0.0,x2,x3);
        matcopy(vx2,0.0,x2,x3);
        matcopy(vz2,0.0,x2,x3);
        matcopy(v_x,0.0,x2,x3);
        matcopy(v_z,0.0,x2,x3);
    }

    void dataClear()
    {
        matcopy(vx1,0.0,x2,x3);
        matcopy(vz1,0.0,x2,x3);
        matcopy(vx2,0.0,x2,x3);
        matcopy(vz2,0.0,x2,x3);
        matcopy(v_x,0.0,x2,x3);
        matcopy(v_z,0.0,x2,x3);
        matcopy(wave,0.0,x1,x2,x3);
    }

};



