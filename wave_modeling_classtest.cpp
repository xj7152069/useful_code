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
    int k,Z(301),X(301);
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
    for(k=0;k<1000;k++)
        {
        f=wavelet01(k,A.dt,35,70);
        A.s2[150][150]=A.s2[150][150]+f;
        A.timeslice();

        if(k%100==0)
            cout<<k<<" && "<<A.s3[150][150]<<endl;
        
        A.s3[0][0]=1.0;
        A.s3[0][1]=-1.0;
        datawrite(A.s3, Z, X, outf1);

        if (k==600)
            {
            matcopy(b,A.s3, Z, X);
            datawrite(b, Z, X, "test_waveform.bin" );
            }
        }

    return 0;
}


