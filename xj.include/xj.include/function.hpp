#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <thread>
#include <future>
#include <time.h>
#include "../xjc.h"
using namespace std;
using namespace arma;

cx_fmat get_blackman_leftwin2d(cx_fmat win, float w);
cx_fmat get_blackman_rightwin2d(cx_fmat win, float w, float wf);
fmat get_blackman_upwin2d(fmat win, float w);
fmat get_blackman_downwin2d(fmat win, float w);
float Blackman(float n, float N);
fmat hilbert1D(fmat s, int n, float dt);
fmat diff1D(fmat s, int n, float dt);
///////////////////////////////////////////////////////////////////

class agc2d{
public:
    fmat data,agcelem,dataagc;
    fcube data3d,agcelem3d;
    int nx,ny,nz,ncpu;
    float stablepar,winnum;
    int agcwx,agcwy,agcwz;
    agc2d(){}
    ~agc2d(){}

agc2d(int n1, int n2, int n3=0){
    nx=n1,ny=n2,nz=n3,ncpu=1;
    agcwx=50,agcwy=2,agcwz=2;
    winnum=(agcwx*2+1)*(agcwy*2+1);
    stablepar=1000;
    if(n3==0){
        this->dataagc.zeros(nx,ny);
        this->data.zeros(nx,ny);
        this->agcelem.zeros(nx,ny);
    }
    else if(n3>0){
        this->data3d.zeros(nx,ny,nz);
        this->agcelem3d.zeros(nx,ny,nz);
    }    
}

void get_agc2d_pthread(agc2d* agc, int n1, int n2)
{
    int i,j,i1,j1;
    float agcpow;
    for(i=n1;i<n2;i++)
    {
    for(j=agc->agcwy;j<agc->ny-agc->agcwy;j++)
    {
        agcpow=0;
        for(i1=-agc->agcwx;i1<=agc->agcwx;i1++)
        {
        for(j1=-agc->agcwy;j1<=agc->agcwy;j1++)
        {
            //agcpow+=abs(data(i+i1,j+j1));
            agcpow+=(data(i+i1,j+j1)*data(i+i1,j+j1));
        }}
        agc->agcelem(i,j)=sqrt(agcpow);
        agcpow=0;
    }}
}

void get_agc2d(fmat & data2, fmat & data)
{
    int i,j,i1,j1;
    float agcpow,maxpow(1.0/this->stablepar);
    for(i=this->agcwx;i<this->nx-this->agcwx;i++)
    {
    for(j=this->agcwy;j<this->ny-this->agcwy;j++)
    {
        agcpow=0;
        for(i1=-this->agcwx;i1<=this->agcwx;i1++)
        {
        for(j1=-this->agcwy;j1<=this->agcwy;j1++)
        {
            agcpow+=(data(i+i1,j+j1)*data(i+i1,j+j1));
        }}
        this->agcelem(i,j)=sqrt(agcpow/this->winnum);
        agcpow=0;
    }}
    for(i=0;i<this->agcwx;i++)
    {
        this->agcelem.row(i)=this->agcelem.row(this->agcwx);
        this->agcelem.row(this->nx-1-i)=this->agcelem.row(this->nx-1-this->agcwx);
    }
    for(i=0;i<this->agcwy;i++)
    {
        this->agcelem.col(i)=this->agcelem.col(this->agcwy);
        this->agcelem.col(this->ny-1-i)=this->agcelem.col(this->ny-1-this->agcwy);
    }

    this->agcelem=this->agcelem/(this->agcelem.max()-this->agcelem.min());
    cout<<this->agcelem.max()<<"|"<<this->agcelem.min()<<endl;
    this->agcelem=this->agcelem+maxpow*(this->agcelem.max()-this->agcelem.min());
    //this->agcelem=this->agcelem+maxpow*(this->agcelem.max());
    this->agcelem=1.0/this->agcelem;
    for(i=0;i<this->nx;i++)
    {
    for(j=0;j<this->ny;j++)
    {
        data2(i,j)=this->agcelem(i,j)*data(i,j);
        this->dataagc(i,j)=this->agcelem(i,j)*data(i,j);
    }}
}

void de_agc2d(fmat & data, fmat & agcelem)
{
    int i,j;
    for(i=0;i<this->nx;i++)
    {
    for(j=0;j<this->ny;j++)
    {
        data(i,j)=data(i,j)/agcelem(i,j);
    }}
}

void add_agc2d(fmat & data, fmat & agcelem)
{
    int i,j;
    for(i=0;i<this->nx;i++)
    {
    for(j=0;j<this->ny;j++)
    {
        data(i,j)=data(i,j)*agcelem(i,j);
    }}
}

void removeagc3d()
{
    int i,j,k;
    for(i=0;i<this->nx;i++)
    {
    for(j=0;j<this->ny;j++)
    {
    for(k=0;k<this->nz;k++)
    {
        this->data3d(i,j,k)=this->data3d(i,j,k)/this->agcelem3d(i,j,k);
    }}}
}

};

void getagc3d_pthread(int coordy, agc2d * par)
{
    int i,j,k,i1,j1,k1;
    float agcpow;
    for(i=par[0].agcwx;i<par[0].nx-par[0].agcwx;i++)
    {
    for(j=coordy;j<(coordy+1);j++)
    {
    for(k=par[0].agcwz;k<par[0].nz-par[0].agcwz;k++)
    {
        agcpow=0;
        for(i1=-par[0].agcwx;i1<=par[0].agcwx;i1++)
        {
        for(j1=-par[0].agcwy;j1<=par[0].agcwy;j1++)
        {
        for(k1=-par[0].agcwz;k1<=par[0].agcwz;k1++)
        {
            agcpow+=abs(par[0].data3d(i+i1,j+j1,k+k1));
        }}}
        par[0].agcelem3d(i,j,k)=agcpow;
        agcpow=0;
    }}}
    for(i=0;i<par[0].agcwx;i++)
    {
        par[0].agcelem3d.row(i)=par[0].agcelem3d.row(par[0].agcwx);
        par[0].agcelem3d.row(par[0].nx-1-i)=par[0].agcelem3d.row(par[0].nx-1-par[0].agcwx);
    }
    for(i=0;i<par[0].agcwz;i++)
    {
        par[0].agcelem3d.slice(i)=par[0].agcelem3d.slice(par[0].agcwz);
        par[0].agcelem3d.slice(par[0].nz-1-i)=par[0].agcelem3d.slice(par[0].nz-1-par[0].agcwz);
    }
}

void getagc3d(agc2d & par)
{
    int coordy,numcpu(par.ncpu),i,j,k;
    agc2d * ppar;
    ppar=&par;
    if(numcpu>(par.ny-par.agcwy-par.agcwy))
    {numcpu=(par.ny-par.agcwy-par.agcwy);}
    thread *pcal;
    pcal=new thread[numcpu];

    int sx1(par.agcwy),sx2(par.ny-par.agcwy),\
        dsx(1),s1(sx1/dsx),s2(sx2/dsx);
    for(i=s1;i<(s2-numcpu);i+=numcpu)
    {
        for(k=0;k<numcpu;k++)
        {
            coordy=(i+k)*dsx;
            pcal[k]=thread(getagc3d_pthread,coordy,ppar);
        }
        for(k=0;k<numcpu;k++)
        {
            pcal[k].join();
        }
    }
    for(i=(s2-numcpu);i<(s2);i++)
    {
        k=i-(s2-numcpu);
        coordy=(i)*dsx;
        pcal[k]=thread(getagc3d_pthread,coordy,ppar);
    }
    for(i=(s2-numcpu);i<(s2);i++)
    {
        k=i-(s2-numcpu);
        pcal[k].join();
    }

    float minpow(1.0/par.stablepar);
    for(i=0;i<par.agcwy;i++)
    {
        par.agcelem3d.col(i)=par.agcelem3d.col(par.agcwy);
        par.agcelem3d.col(par.ny-1-i)=par.agcelem3d.col(par.ny-1-par.agcwy);
    }
    par.agcelem3d=par.agcelem3d/(par.agcelem3d.max()-par.agcelem3d.min());
    cout<<par.agcelem3d.max()<<"|"<<par.agcelem3d.min()<<endl;
    par.agcelem3d=par.agcelem3d+minpow*(par.agcelem3d.max()-par.agcelem3d.min());
    par.agcelem3d=1.0/par.agcelem3d;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
    for(k=0;k<par.nz;k++)
    {
        par.data3d(i,j,k)=par.agcelem3d(i,j,k)*par.data3d(i,j,k);
    }}}
}

/*
void get_agc2d_thread(agc2d & agc)
{
    int i,j,i1,j1;
    float agcpow,maxpow(1.0/agc.stablepar);
    thread *pcal;
    pcal=new thread[agc.ncpu];
    for(i=0;i<agc.ncpu;i++){

        pcal[i]=thread(get_agc2d_pthread,&par,i1,j2);
    }

    for(i=0;i<agc.agcwx;i++)
    {
        agc.agcelem.row(i)=this->agcelem.row(this->agcwx);
        this->agcelem.row(this->nx-1-i)=this->agcelem.row(this->nx-1-this->agcwx);
    }
    for(i=0;i<this->agcwy;i++)
    {
        this->agcelem.col(i)=this->agcelem.col(this->agcwy);
        this->agcelem.col(this->ny-1-i)=this->agcelem.col(this->ny-1-this->agcwy);
    }

    this->agcelem=this->agcelem/(this->agcelem.max()-this->agcelem.min());
    cout<<this->agcelem.max()<<"|"<<this->agcelem.min()<<endl;
    this->agcelem=this->agcelem+maxpow*(this->agcelem.max()-this->agcelem.min());
    this->agcelem=1.0/this->agcelem;
    for(i=0;i<this->nx;i++)
    {
    for(j=0;j<this->ny;j++)
    {
        data2(i,j)=this->agcelem(i,j)*data(i,j);
        this->dataagc(i,j)=this->agcelem(i,j)*data(i,j);
    }}
}
*/
///////////////////////////////////////////////////////////////////

fmat get_agc2d(fmat & data2, fmat & data, int agcwx=50, int agcwy=5)
{
    int i,j,i1,j1,nx(data.n_rows),ny(data.n_cols);
    fmat agcelem;
    agcelem.copy_size(data);
    float agcpow,maxpow(0.001);
    for(i=agcwx;i<nx-agcwx;i++)
    {
    for(j=agcwy;j<ny-agcwy;j++)
    {
        agcpow=0;
        for(i1=-agcwx;i1<=agcwx;i1++)
        {
        for(j1=-agcwy;j1<=agcwy;j1++)
        {
            agcpow+=(data(i+i1,j+j1)*data(i+i1,j+j1));
        }}
        agcelem(i,j)=sqrt(agcpow);
        agcpow=0;
    }}
    for(i=0;i<agcwx;i++)
    {
        agcelem.row(i)=agcelem.row(agcwx);
        agcelem.row(nx-1-i)=agcelem.row(nx-1-agcwx);
    }
    for(i=0;i<agcwy;i++)
    {
        agcelem.col(i)=agcelem.col(agcwy);
        agcelem.col(ny-1-i)=agcelem.col(ny-1-agcwy);
    }

    agcelem=agcelem/(agcelem.max());
    agcelem=agcelem+maxpow;
    //this->agcelem=this->agcelem+maxpow*(this->agcelem.max());
    agcelem=1.0/agcelem;
    for(i=0;i<nx;i++)
    {
    for(j=0;j<ny;j++)
    {
        data2(i,j)=agcelem(i,j)*data(i,j);
    }}
    return agcelem;
}
void de_agc2d(fmat & data, fmat & agcelem)
{
    int i,j;
    int nx(data.n_rows),ny(data.n_cols);
    for(i=0;i<nx;i++)
    {
    for(j=0;j<ny;j++)
    {
        data(i,j)=data(i,j)/agcelem(i,j);
    }}
}
fmat diff1D(fmat s, int n, float dt)
{
    fmat ss(n,1);
    int i,j;
    for(i=1;i<n-1;i++)
    {
        ss(i,0)=(s(i+1,0)-s(i-1,0))/dt/2;
    }

    return ss;
}

fmat hilbert1D(fmat s, int n, float dt)
{
    float *h, pi(3.1415926);
    int i,j,z,hn;
    hn=3/dt;
    z=int(hn/2);
    h=new float[hn];
    for(i=0;i<hn;i++)
    {
        if((i-z)!=0)
        {
            h[i]=1/(pi*dt*(i-z));
        }
        else
        {
            h[i]=0.0;
        }   
    }

    float *h2, *s2, *sh, *sh2;
    s2=new float[hn+n];
    h2=new float[hn+n];
    sh=new float[hn+n];
    sh2=new float[hn+n];
    for(i=0;i<hn+n;i++)
    {
        s2[i]=0.0;
        h2[i]=0.0;
        sh[i]=0.0;
        sh2[i]=0.0;
        if(i<n)
        {
            s2[i]=s(i,0);
        }
        if(i<hn)
        {
            h2[i]=h[i];
        }
    }

    for(i=0;i<hn+n;i++)
    {
        for(j=0;j<=i;j++)
        {
            sh[i]+=s2[i-j]*h2[j];
        }
        if(i>=(z))
        {
            sh2[i-z]=sh[i]/(1.0/dt);
        }
    }
    fmat ss(n,1);
    for(i=0;i<n;i++)
    {
        ss(i,0)=sh2[i];
    }
    delete [] h2;
    delete [] s2;
    delete [] sh;
    delete [] sh2;
    delete [] h;
    return ss;
}

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
cx_fmat get_blackman_downwin2d(cx_fmat win, float w)
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
                (win(i,j)).real(real(win(i,j))*n);
                (win(i,j)).imag(imag(win(i,j))*n);
            }
        } 
    }
    return win;
}

cx_fmat get_blackman_upwin2d(cx_fmat win, float w)
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
////////////function: tx2fx or fx2tx/////////////
void tx2fx_3d_pthread(cx_fcube *data3d_out,fcube *data3d,int i)
{
    int j,n1(data3d[0].n_rows),n2(data3d[0].n_cols),n3(data3d[0].n_slices);
    cx_fmat datafx(n2,n3);
    fmat data(n2,n3);
    data=data3d[0].row(i);
    for(j=0;j<n2;j++)
    {
        datafx.row(j)=fft(data.row(j),n3);
    }
    data3d_out[0].row(i)=datafx;
}
void tx2fx_3d_thread(cx_fcube &data3d_out,fcube &data3d, int ncpu)
{
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices);
    ncpu=min(ncpu,n1);
    cx_fcube *pdata3d_out=&data3d_out;
    fcube *pdata3d=&data3d;
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=0;i<(n1-ncpu);i+=ncpu){
        for(k=0;k<ncpu;k++){
            pcal[k]=thread(tx2fx_3d_pthread,pdata3d_out,pdata3d,i+k);
        }
        for(k=0;k<ncpu;k++){
            if(pcal[k].joinable()){
                pcal[k].join();}
        }
    } 
    for(i=(n1-ncpu);i<(n1);i++){
        k=i-(n1-ncpu);
        pcal[k]=thread(tx2fx_3d_pthread,pdata3d_out,pdata3d,i);
    }
    for(i=(n1-ncpu);i<(n1);i++){
        k=i-(n1-ncpu);
        if(pcal[k].joinable()){
            pcal[k].join();}
    }
    delete [] pcal;
}

void fx2tx_3d_pthread(fcube *data3d_out,cx_fcube *data3d,int i)
{
    int j,n1(data3d[0].n_rows),n2(data3d[0].n_cols),n3(data3d[0].n_slices);
    cx_fmat data(n2,n3),datafx(n2,n3);
    datafx=data3d[0].row(i);
    for(j=0;j<n2;j++)
    {
        data.row(j)=ifft(datafx.row(j),n3);
    }
    data3d_out[0].row(i)=real(data);
}
void fx2tx_3d_thread(fcube &data3d_out,cx_fcube &data3d, int ncpu)
{
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices);
    ncpu=min(ncpu,n1);
    fcube *pdata3d_out=&data3d_out;
    cx_fcube *pdata3d=&data3d;
    thread *pcal;
    pcal=new thread[ncpu];
    
    for(i=0;i<(n1-ncpu);i+=ncpu){
        for(k=0;k<ncpu;k++){
            pcal[k]=thread(fx2tx_3d_pthread,pdata3d_out,pdata3d,i+k);
        }
        for(k=0;k<ncpu;k++){
            if(pcal[k].joinable()){
                pcal[k].join();}
        }
    } 
    for(i=(n1-ncpu);i<(n1);i++){
        k=i-(n1-ncpu);
        pcal[k]=thread(fx2tx_3d_pthread,pdata3d_out,pdata3d,i);
    }
    for(i=(n1-ncpu);i<(n1);i++){
        k=i-(n1-ncpu);
        if(pcal[k].joinable()){
            pcal[k].join();}
    }
    delete [] pcal;
}
void medianfilter(float **mat, int n1, int n2, int w1=25, int w2=1, float lmd=5.0)
{
    fmat data(n1,n2),datapow(n1,n2),datapowfilter(n1,n2);
    int i,j,k1,k2,k3,num((1+2*w1)*(1+2*w2));
    float *win,avera;
    fmat datawin(num,1);

    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            data(i,j)=mat[i][j];
        }
    }
    datapow=abs(data);
    datapowfilter=datapow;
    for(i=w1;i<n1-w1;i++){
    for(j=w2;j<n2-w2;j++){
        k3=0;
    for(k1=-w1;k1<=w1;k1++){
    for(k2=-w2;k2<=w2;k2++){
        datawin(k3,0)=datapow(i+k1,j+k2);
        k3++;
    }}
    avera=sum(sum(datawin))/num;
    if(datapow(i,j)>=lmd*avera)
        datapowfilter(i,j)=avera;
    }}
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            mat[i][j]=mat[i][j]*datapowfilter(i,j)/(datapow(i,j)+0.00001);
        }
    }
}
void medianfilter(fmat & mat, int n1, int n2, int w1=25, int w2=1, float lmd=5.0)
{
    fmat data(n1,n2),datapow(n1,n2),datapowfilter(n1,n2);
    int i,j,k1,k2,k3,num((1+2*w1)*(1+2*w2));
    float *win;
    fmat datawin(num,1),avera(1,1);

    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            data(i,j)=mat(i,j);
        }
    }
    datapow=abs(data);
    datapowfilter=datapow;
    for(i=w1;i<n1-w1;i++){
    for(j=w2;j<n2-w2;j++){
        k3=0;
    for(k1=-w1;k1<=w1;k1++){
    for(k2=-w2;k2<=w2;k2++){
        datawin(k3,0)=datapow(i+k1,j+k2);
        k3++;
    }}
    avera=sum(sum(datawin,1),0)/float(num);
    if(datapow(i,j)>=lmd*avera(0,0))
        datapowfilter(i,j)=avera(0,0);
    }}
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++){
            mat(i,j)=mat(i,j)*datapowfilter(i,j)/(datapow(i,j)+0.00001);
        }
    }
}
#endif


