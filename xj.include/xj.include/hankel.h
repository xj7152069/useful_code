/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/
#ifndef HANKEL_H_H
#define HANKEL_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;


void output_qyzp(fvec & s, int n, char const *str1 )
{
    char str[99];
    strcpy(str, str1);
    int i,j,xs(1);
    double max_s,min_s;
    float outdata(0.0);
    min_s=min(s);
    max_s=max(s);
    cout<<"max="<<max_s<<" ; min="<<min_s<<endl;
    max_s=max_s-min_s;
    while(max_s>100)
    {
        max_s=max_s/2;
        xs=xs*2;
    }
    max_s+=10;
    //fmat qyzp(max_s, n);
    //qyzp.fill(0.0);
    ofstream outf;
    outf.open(str);
    for(j=0;j<=n;j++)
    {
        for(i=0;i<int(max_s);i++)
        {
            if(j>0 && i==int((s(j-1)-min_s)/xs))
                {
                    outdata=10.0;
                    outf.write((char*)&outdata,sizeof(outdata));
                }
            else
                {
                    outdata=0.0;
                    outf.write((char*)&outdata,sizeof(outdata));
                }
        }
    }
    outf.close(); 
    cout<<"n1 = "<<int(max_s)<<endl;
}

cx_fmat mathankel(cx_fmat &mat1, int nh, int k)
{
    int N;
    N=nh;
    cx_fmat hankel(N,N);
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            hankel(i,j).real(real(mat1(k,i+j)));
            hankel(i,j).imag(imag(mat1(k,i+j)));
        }
    }
    return hankel;
}

void ihankelaverage(cx_fmat &hankel, int nh, cx_fmat &ihankel, int k2)
{
    int N;
    N=nh;
    float realavg,imagavg;
    int i,j,k,n;
    for(k=0;k<2*N-1;k++)
    {
        realavg=0.0;
        imagavg=0.0;
        n=0;
        for(i=0;i<=k;i++)
        {
            j=k-i;
            if(i<N && j<N)
            {
                realavg=realavg+real(hankel(i,j));
                imagavg=imagavg+imag(hankel(i,j));
                n=n+1;
            }
        }
        realavg=realavg/n;
        imagavg=imagavg/n;
        ihankel(k2,k).imag(imagavg);
        ihankel(k2,k).real(realavg);
    }
}

void ihankel0(cx_fmat &hankel, int nh, cx_fmat &ihankel, int k2)
{
    int N;
    N=nh;
    float realavg,imagavg;
    int i,j,k,n;
    for(k=0;k<2*N-1;k++)
    {
        if(k<nh)
        {
            ihankel(k2,k).imag(imag(hankel(1,k)));
            ihankel(k2,k).real(real(hankel(1,k)));
        }
        else
        {
            ihankel(k2,k).imag(imag(hankel(nh-1,k-nh+1)));
            ihankel(k2,k).real(real(hankel(nh-1,k-nh+1)));
        }
        
    }
}

#endif
