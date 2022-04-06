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
void multiple_code(fmat& u2, fmat& u1, fmat& green,\
 float df, int fn1=0, int fn2=2048);
///////////////////////////////////////////////////////////////////
fmat single_trace_dewave(fmat & data, fmat & dict,\
    fmat &dataph, fmat &datapd, fmat &dataphd,\
    int nt,int nx,int nwp,int nwph,int nwpd,int nwphd,\
    int nwt,int dnwt,float zzh=0.001)
{
    int i,j,k,wt1,kwt,nw(nwp+nwph+nwpd+nwphd);
    fmat mat1(nwt,nw),mat0(nwt,nw),datap3(nt,nx,fill::zeros),\
        matq(nw,1),mat2(nw,1),matd(nwt,1),matD(nw,nw);
    fmat digmat(nw,nw,fill::zeros),datal(nt,1);
    fmat datap(nt,nx),datavz(nt,nx);

////////////////////////////////////////////////////////////////////
    
    digmat.diag()+=1;
    datap=dict/(dict.max()-dict.min());
    datavz=data/(data.max()-data.min());
    dataph=dataph/(dataph.max()-dataph.min());
    datapd=datapd/(datapd.max()-datapd.min());
    dataphd=dataphd/(dataphd.max()-dataphd.min());

    for(k=0;k<nx;k++){
        cout<<"tracl = "<<k<<endl;
        for(kwt=0;kwt<=(nt-nwt);kwt+=dnwt){
            mat1.fill(0.0);

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwp;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=datap(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwp;j++){
                mat1.col(j-kwt)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwph;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=dataph(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwph;j++){
                mat1.col(j-kwt+nwp)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwpd;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=datapd(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwpd;j++){
                mat1.col(j-kwt+nwp+nwph)=mat0.col(j-kwt);
            }

            mat0.fill(0.0);
            for(i=kwt;i<kwt+nwt;i++){
                for(j=kwt;j<kwt+nwphd;j++){
                    if((i-j)>=0)
                        mat0(i-kwt,j-kwt)=dataphd(i-j+kwt,k);
                }
            }
            for(j=kwt;j<kwt+nwphd;j++){
                mat1.col(j-kwt+nwp+nwph+nwpd)=mat0.col(j-kwt);
            }
            ///////////
            
            for(i=0;i<nwt;i++){
                matd(i,0)=datavz(i+kwt,k);
            }
            matD=mat1.t()*mat1;
            
            digmat.fill(0);
            digmat.diag()+=(matD.max()*zzh+zzh);
            //cout<<digmat(1,1)<<endl;
            matq=inv(matD+digmat)*mat1.t()*matd;
            matd=mat1*matq;
            for(i=0;i<nwt;i++){
                datap3(i+kwt,k)+=datavz(i+kwt,k)-matd(i,0);
            }
        }
    }
    for(i=0;i<nt;i++){
        for(j=0;j<nx;j++){
            if(isnan(abs(datap3(i,j))))
                datap3(i,j)=0;
        }
    }

    return datap3;
}

void multiple_code(fmat& u2, fmat& u1, fmat& green,\
 float df, int fn1, int fn2)
{
    int n1(u1.n_rows),n2(u1.n_cols);
    int i,j,k;
    float w,pi(3.1415926);
    cx_fmat codemat(n2,n2,fill::zeros),onecol(n1,1),\
        u1cx(n2,n1),u2cx(n2,n1,fill::zeros);
    cx_float a;
    a.real(0.0);

    for(k=0;k<n2;k++){
        onecol.col(0)=fft(u1.col(k),n1);
        u1cx.row(k)=onecol.col(0).st();
        //u2cx.row(k)=fft(u2.col(k).t(),n1);
    }
    
    for(k=fn1;k<fn2;k++){
        w=-2.0*pi*df*k;
        cout<<"k="<<k<<" | "<<"w="<<w<<endl;
        for(i=0;i<n2;i++){
            for(j=0;j<n2;j++){
                a.imag(w*green(i,j));
                codemat(i,j)=exp(a);
            }
        }
        u2cx.col(k)=codemat*u1cx.col(k);
    }

    for(k=0;k<n2;k++){
        //u1cx.row(k)=fft(u1.col(k).t(),n1);
        onecol.col(0)=u2cx.row(k).st();
        u2.col(k)=real(ifft(onecol.col(0),n1));
    }
}

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

void removeagc3d(agc2d & par)
{
    int i,j,k;
    for(i=0;i<par.nx;i++)
    {
    for(j=0;j<par.ny;j++)
    {
    for(k=0;k<par.nz;k++)
    {
        par.data3d(i,j,k)=par.data3d(i,j,k)/par.agcelem3d(i,j,k);
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


