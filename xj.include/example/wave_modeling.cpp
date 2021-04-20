#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

#include "/home/xj/github-xj7152069/useful_code/xj.include/xjc.h"

int main()
{
    int nz(300), nx(500), nt(3000), i,j,k;
    wave2D s(nz,nx),s2(nz,nx);
    for(i=0;i<130;i++)
    {for(j=0;j<nx;j++)s.model[i][j]=2000.0;}
    for(i=130;i<220;i++)
    {for(j=0;j<nx;j++)s.model[i][j]=2500.0;}
    for(i=220;i<300;i++)
    {for(j=0;j<nx;j++)s.model[i][j]=5000.0;}
    s2.setvelocity(2000);
    s.updatepar();
    s2.updatepar();
    float **suf1, **suf2;
    suf1=newfmat(nt,nx);
    suf2=newfmat(nt,nx);
    for(k=0;k<nt;k++)
    {
        s.s2[50][50]+=wavelet01(k,s.dt,25);
        s2.s2[50][50]+=wavelet01(k,s.dt,25);
        s.timeslicecal();
        s2.timeslicecal();

        matcopy(suf1[k],s.s2[50],nx);
        matcopy(suf2[k],s2.s2[50],nx);
        
        if(k%100==0)
        cout<<k<<endl;
    }

    ofstream outf;
    matdec(suf1,suf2,nt,nx);
    datawrite(suf1,nt,nx,"data2.bin");

    outf.open("data2.txt",ios::out);
    for(j=0;j<nx;j++)
            {for(i=0;i<nt;i++)
                    {
                    outf<<suf1[i][j]<<' ';
                    }
                    outf<<endl;
            }
    outf.close();


    return 0;
}





