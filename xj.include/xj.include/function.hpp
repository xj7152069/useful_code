#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

#include <time.h>
#include "../xjc.h"
using namespace std;
using namespace arma;
cx_fmat get_blackman_leftwin2d(cx_fmat win, float w);
cx_fmat get_blackman_rightwin2d(cx_fmat win, float w, float wf);
fmat get_blackman_upwin2d(fmat win, float w);
fmat get_blackman_downwin2d(fmat win, float w);
float Blackman(float n, float N);
///////////////////////////////////////////////////////////////////
float Blackman(float n, float N)
{
    float xs, pi(3.1415926);
    xs=0.42-0.5*cos(2*pi*n/N/2)+0.08*cos(4*pi*n/N/2);
    return xs;
}

cx_fmat get_blackman_leftwin2d(cx_fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(j<=w)
            {
                n=j;
                n=Blackman(n,w);
                (win(i,j)).real(real(win(i,j))*n);
                (win(i,j)).imag(imag(win(i,j))*n);
            }
        } 
    }
    return win;
}
cx_fmat get_blackman_rightwin2d(cx_fmat win, float w, float wf)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(j>=(wf-w) && j<wf)
            {
                n=wf-1-j;
                n=Blackman(n,w);
                (win(i,j)).real(real(win(i,j))*n);
                (win(i,j)).imag(imag(win(i,j))*n);
            }
        } 
    }
    return win;
}
fmat get_blackman_leftwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(j<=w)
            {
                n=j;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}
fmat get_blackman_rightwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(j>=n2-w-1)
            {
                n=n2-1-j;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

fmat get_blackman_downwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i>=n1-w-1)
            {
                n=n1-1-i;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

fmat get_blackman_upwin2d(fmat win, float w)
{
    int n1,n2;
    n1=win.n_rows;
    n2=win.n_cols;
    int i,j;
    float n;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(i<=w)
            {
                n=i;
                n=Blackman(n,w);
                win(i,j)*=n;
            }
        } 
    }
    return win;
}

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


