#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

#include "../xjc.h"

int srand_seed(666);

template <typename T1> \
void de_mean_to_zero1d(T1 * s, int n)
{
    int i;
    float meannum(0);
    for(i=0;i<n;i++)
    {
        meannum+=s[i];
    }
    meannum=meannum/n;
    for(i=0;i<n;i++)
    {
        s[i]=s[i]-meannum;
    }
}

template <typename T1> \
float** get_cov_mat1d(T1 * ss1,T1 * ss2, int n)
{
    int i,j,k;
    float** covmat;
    covmat=newfmat(n,n);
    T1 *s1,*s2;
    s1=new T1[2*n];
    s2=new T1[2*n];
    matcopy(covmat,0.0,n,n);
    //de_mean_to_zero1d(s,n);
    for(i=0;i<n;i++)
    {
        s1[i]=ss1[i],s1[i+n]=ss1[i];
        s2[i]=ss2[i],s2[i+n]=ss2[i];
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n;k++)
            {
                covmat[i][j]+=s1[i+k]*s2[j+k];
            }
        }
    }
    delete []s1;
    delete []s2;
    return covmat;
}


template <typename T1> \
float** get_selfcov_mat1d(T1 * s, int n)
{
    int i,j,k;
    float** covmat;
    covmat=newfmat(n,n);
    T1 *s1,*s2;
    s1=new T1[2*n];
    s2=new T1[2*n];
    matcopy(covmat,0.0,n,n);
    de_mean_to_zero1d(s,n);
    for(i=0;i<n;i++)
    {
        s1[i]=s[i],s1[i+n]=s[i];
        s2[i]=s[i],s2[i+n]=s[i];
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n;k++)
            {
                covmat[i][j]+=s1[i+k]*s2[j+k];
            }
        }
    }
    delete []s1;
    delete []s2;
    return covmat;
}

float getrandonnoise(float N=1.0)
{
    float k,n;
	if(srand_seed==666)
    {
        srand((int)time(0));  /*根据当前时间设置“随机数种子”*/
    }
    else
    {
        srand(srand_seed);
    }   
    k=rand();
    srand_seed=int(k)%1000000000;
    n=(float(int(k)%1000000)*0.000001*N-N/2);
    return n;
}

float getgaussonnoise(float N=1.0, float p=1.0)
{
    float n,k;
    if(srand_seed==666)
    {
        srand((int)time(0));  /*根据当前时间设置“随机数种子”*/
    }
    else
    {
        srand(srand_seed);
    }   
    k=rand();
    srand_seed=int(k)%1000000000;
    //cout<<k<<endl;
    n=(float(int(k)%1000000)*0.000001-0.5);
    n=N*exp(-(n*n/2.0/p));
    return n;
}

#endif


