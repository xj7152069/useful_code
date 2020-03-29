/*********(version 1.0)***********/
/*注释：
    C++程序模板：

*/
/********************************/

#include <iostream>
using namespace std;

#include "xj_c++.h"

float zb(int k, float DT)
{
    float pi(3.1415926);
	float f;
    f=(pi)*(pi)*900*(k*DT-0.04)*\
exp((-pi*pi*900*(k*DT-0.04)*(DT*k-0.04)))\
*(3-2*pi*pi*900*(k*DT-0.04)*(DT*k-0.04));
	return f;
}


int main ()
{
    int k;
    float f;

    wave_modeling_2D A(200,200);
    A.dx=5.5;
    A.dy=1.0;
    A.dt=0.0005;
    for(k=0;k<800;k++)
    {
    f=zb(k,A.dt);
    A.s2[100][100]=A.s2[100][100]+f;
    A.wave_modeling_2D_time_pice();

    if(k%100==0)
        cout<<k<<" && "<<A.s3[100][100]<<endl;
    
    if (k==450)
        binary_data_write_2D(A.s3, 200, 200, "test_waveform.bin" );
    }
 
    return 0;
}
