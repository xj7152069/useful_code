/*********(version 1.0)***********/
/*注释：
    C++程序模板：

*/
/********************************/

#include <iostream>
using namespace std;

#include "xj_c++.h"
#include "xj_c++_armadillo.h"

int main ()
{
    int k,Z(201),X(201);
    float f;
    ofstream outf1;
    outf1.open("movie.bin");

    wave2D A(Z,X);
    float **b;
    b=new float*[Z];      
    for(k=0;k<Z;k++)  
        {  
        b[k]=new float[X]; 
        } 
    matcopy(b,0.0,Z,X);
/*
    A.dx=4.0;
    A.dy=5.0;
    A.dt=0.0006;
    A.PML_wide=20;
    A.R=100000;
*/
    for(k=0;k<800;k++)
        {
        f=-wavelet02(k,A.dt,35,70);
        A.s2[100][100]=A.s2[100][100]+f;
        A.timeslice();

        if(k%100==0)
            cout<<k<<" && "<<A.s3[100][100]<<endl;
        
/*        if(abs(A.s3[0][0])<1000*abs(f))
            {
            A.s3[0][1]=-1000*abs(f);
            A.s3[0][0]=1000*abs(f);
            }
*/
        datawrite(A.s3, Z, X, outf1);

        if (k==400)
            {
            matcopy(b,A.s3, Z, X);
            datawrite(b, Z, X, "test_waveform.bin" );
            }
        }

    return 0;
}


