/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/
#ifndef RTMANGLE2D_H_H
#define RTMANGLE2D_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

template<typename T1, typename T2> float anglecal_z1_x0(T1 z, T2 x, float err=0.000000000001);
extern float** newfmat(int x1, int x2);
extern float*** newfmat(int x1, int x2, int x3);
template <typename T1> extern void matdelete(T1 **mat, int x1);
template <typename T1> extern void matdelete(T1 ***mat, int x1, int x2);
template <typename T1, typename T2> extern void matcopy(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> extern void matcopy(T1 ***mat1, T2 n, int n1, int n2, int n3);
template <typename T1, typename T2> extern void matcopy(T1 **mat1, T2 **mat2, int nz, int nx);
fmat windowsZ(int Z, int X, float thetaN);
void fft2dwindows(float **ang, float **win, int x1, int x2, float theta, float thetaN);
fmat fft2dwintransform(float **win, int Z, int X);
float angletransform(float a);
void getangle1(float **angle1, int Z, int X);
void getangle2(float **angle1, int Z, int X);
void fft2dwindows_Gabor(float **ang, float **win, int x1, int x2, float theta, float thetaN);
void fft2dwindows_Blackman(float **ang, float **win, int x1, int x2, float theta, float thetaN);
float Gabor(float da,float a);
////////////////////////////////////////////////////////////////////////////
template<typename T1, typename T2>
float anglecal_z1_x0(T1 z, T2 x, float err)
{
    float z0(1.0), x0(0.0), pi(3.1415926), xs, theta, c, l;
    xs=180/pi;
    l=sqrt(z*z+x*x);
    if(l>err)
    {
        c=(z0*z+x0*x)/l;
        theta=acos(c)*xs;
        
        if(x<0.0)
            {theta=-theta;}
            
    }
    else
    {
        theta=0.0;
    }
    return theta;
}

////////////////////////////////////////////////////////
class angle_gather2D
{
private:
    float da, fa, recva;

public:
    float ***DA, ***FA, ***recvA, a_beg, a_gep;
    float **ps, **pr;
    int x1, x2, x3;

    angle_gather2D()
    {
        x1=0;x2=0;x3=0;
        da=0;fa=0;
        cout<<"Warning: Create an Empty object!"<<endl;
    }
    angle_gather2D(int na, int nz, int nx, float A_BEG=-180, float A_GEP=1.0)
    {
        x1=na;x2=nz;x3=nx;
        ps=newfmat(x2,x3);
        pr=newfmat(x2,x3);
        da=0;fa=0;recva=0;
        a_beg=A_BEG;
        a_gep=A_GEP;
    }
    ~angle_gather2D()
    {
        matdelete(ps,x2);
        matdelete(pr,x2);
    }

    template<typename T1, typename T2>
    void theatcal_S(T1 **svz, T2 **svx)
    {
        int i,j;
        for(i=0;i<x2;i++)
        {
            for(j=0;j<x3;j++)
            {
                ps[i][j]=anglecal_z1_x0(svz[i][j],svx[i][j]);
            }
        }
    }

    template<typename T1, typename T2>
    void theatcal_R(T1 **rvz, T2 **rvx)
    {
        int i,j;
        for(i=0;i<x2;i++)
        {
            for(j=0;j<x3;j++)
            {
                pr[i][j]=anglecal_z1_x0(rvz[i][j],rvx[i][j]);
            }
        }
    }

    template<typename T1, typename T2>
    void addFA(T1 **swave, T2 **rwave)
    {
        if(fa==0)
        {
            FA=newfmat(x3,x2,x1);
            fa=1;
            matcopy(FA, 0.0, x3, x2, x1);
        }
        int i,j,k,ang;
        for(i=0;i<x2;i++)
        {
            for(j=0;j<x3;j++)
            {
                ang=int((ps[i][j]-pr[i][j]-a_beg)/a_gep);
                if(ang<x1 && ang>0)
                {
                    FA[j][i][ang]+=swave[i][j]*rwave[i][j];
                }
            }
        }
    }

    template<typename T1, typename T2>
    void addDA(T1 **swave, T2 **rwave)
    {
        if(da==0)
        {
            DA=newfmat(x3,x2,x1);
            da=1;
            matcopy(DA, 0.0, x3, x2, x1);
        }
        int i,j,k,ang;
        for(i=0;i<x2;i++)
        {
            for(j=0;j<x3;j++)
            {
                ang=int(((ps[i][j]+pr[i][j])/2.0-a_beg)/a_gep);
                if(ang<x1 && ang>0)
                {
                    DA[j][i][ang]+=swave[i][j]*rwave[i][j];
                }
            }
        }
    }

    template<typename T1, typename T2>
    void addrecvA(T1 **swave, T2 **rwave)
    {
        if(recva==0)
        {
            recvA=newfmat(x3,x2,x1);
            recva=1;
            matcopy(recvA, 0.0, x3, x2, x1);
        }
        int i,j,k,ang;
        for(i=0;i<x2;i++)
        {
            for(j=0;j<x3;j++)
            {
                ang=int((pr[i][j]-a_beg)/a_gep);
                //cout<<ang<<" | ";
                if(ang<x1 && ang>=0)
                {
                    recvA[j][i][ang]+=swave[i][j]*rwave[i][j];
                }
            }
        }
    }

};

//////////////////////////////////////////////////

class fft2d_anglegather2d
{
public:
    float ***DA, ***FA, da, fa, a_beg, recva_beg, a_gep, recva_gep, a2_err, ***recvA;
    float ***win, *winangle, **ps, **angle1, recvid;
    int x1, x2, x3, x2fft, x3fft, nwin, recva_num;
    fmat *wintrans;

    fft2d_anglegather2d()
    {
        cout<<"Warning: Create an Empty object!"<<endl;
        recvid=0;
    }
    fft2d_anglegather2d(int na, int nz, int nx, int nzfft, int nxfft, float A_BEG=-180, float A_GEP=1.0)
    {
        x1=na;x2=nz;x3=nx;
        x2fft=nzfft;x3fft=nxfft;
        a_beg=A_BEG;
        a_gep=A_GEP;
        da=0;fa=0;
        win=newfmat(x1,x2fft,x3fft);
        //anglewave=newfmat(x1,x2fft,x3fft);
        ps=newfmat(x2,x3);
        FA=newfmat(x3,x2,x1);
        matcopy(FA, 0.0, x3, x2, x1);
        DA=newfmat(x3,x2,x1);
        matcopy(DA, 0.0, x3, x2, x1);
        recvid=0;
    }
    ~fft2d_anglegather2d()
    {
        cout<<"you delete an object-fft2d_anglegather2d"<<endl;
    }

    void createanglewin(int Nwin, float A_err ,float *Winangle)
    {
        nwin=Nwin;
        winangle=Winangle;
        a2_err=A_err;
        win=newfmat(nwin,x2fft,x3fft);
        angle1=newfmat(x2fft,x3fft);
        getangle2(angle1,x2fft,x3fft);
        float a2;
        int k;
        ofstream outf;
        outf.open("win.bin");
        this->wintrans=new fmat[nwin];
        for(k=0;k<nwin;k++)
        {
            matcopy(win[k],0.0,x2fft,x3fft);
            //a2=angletransform(winangle[k]);
            //fft2dwindows_Gabor(angle1,win[k],x2fft,x3fft,winangle[k],a2_err);
            fft2dwindows_Blackman(angle1,win[k],x2fft,x3fft,winangle[k],a2_err);
            this->wintrans[k].resize(x2fft,x3fft);
            this->wintrans[k]=fft2dwintransform(win[k],x2fft,x3fft);
            datawrite(wintrans[k],x2fft,x3fft,outf);
        }
        outf.close();
    }

    void creatRecvAngle(int na, float anglebeg, float anglegep)
    {
        recva_num=na;
        recvA=newfmat(x3,x2,recva_num);
        matcopy(recvA, 0.0, x3, x2, recva_num);
        recva_beg=anglebeg;
        recva_gep=anglegep;
        recvid=1;
        cout<<"You have create a Recv-Angle-Gather"<<endl;
    }

    template<typename T1, typename T2>
    void theatcal_S(T1 **svz, T2 **svx)
    {
        int i,j;
        for(i=0;i<x2;i++)
        {
            for(j=0;j<x3;j++)
            {
                ps[i][j]=anglecal_z1_x0(svz[i][j],svx[i][j]);
            }
        }
    }

    template<typename T1>
    void addFA_DA(T1 **swave, cx_fmat rwave)  //S-wave:pyt ; R-wave:fft2d
    {
        int i,j,k,ang;float imagpower;
        cx_fmat wfft(x2fft,x3fft),wfftwin(x2fft,x3fft),wifft(x2fft,x3fft);
        fmat copyreal(x2fft,x3fft),copyimag(x2fft,x3fft),copy(x2fft,x3fft);
        wfft=fft2(rwave,x2fft,x3fft);
        copyreal=real(wfft);
        copyimag=imag(wfft);
        ofstream outf;
        //outf.open("fft2danglewin.bin");
        for(k=0;k<nwin;k++)
        {   
            copy=matmul(copyreal,this->wintrans[k],x2fft,x3fft);
            wfftwin.set_real(copy);
            copy=matmul(copyimag,this->wintrans[k],x2fft,x3fft);
            wfftwin.set_imag(copy);
            wifft=ifft2(wfftwin,x2fft,x3fft);
            copy=real(wifft);
            //datawrite(copy,x2fft,x3fft,outf);
            for(i=0;i<x2;i++)
            {
                for(j=0;j<x3;j++)
                {
                    imagpower=swave[i][j]*copy(i,j);
                    ang=int((ps[i][j]-winangle[k]-a_beg)/a_gep);
                    if(ang<x1 && ang>0)
                    {
                        FA[j][i][ang]+=imagpower;
                    }
                    ang=int(((ps[i][j]+winangle[k])/2.0-a_beg)/a_gep);
                    if(ang<x1 && ang>0)
                    {
                        DA[j][i][ang]+=imagpower;
                    }
                    if (recvid==1)
                    {
                        recvA[j][i][k]+=imagpower;
                    }
                }
            }
        }
        //outf.close();
    }

    template<typename T1>
    void add_recvAngleGather(T1 **swave, cx_fmat rwave)  //S-wave:pyt ; R-wave:fft2d
    {
        int i,j,k,ang;float imagpower;
        cx_fmat wfft(x2fft,x3fft),wfftwin(x2fft,x3fft),wifft(x2fft,x3fft);
        fmat copyreal(x2fft,x3fft),copyimag(x2fft,x3fft),copy(x2fft,x3fft);
        wfft=fft2(rwave,x2fft,x3fft);
        copyreal=real(wfft);
        copyimag=imag(wfft);
        ofstream outf;
        //outf.open("fft2danglewin.bin");
        for(k=0;k<nwin;k++)
        {   
            copy=matmul(copyreal,this->wintrans[k],x2fft,x3fft);
            wfftwin.set_real(copy);
            copy=matmul(copyimag,this->wintrans[k],x2fft,x3fft);
            wfftwin.set_imag(copy);
            wifft=ifft2(wfftwin,x2fft,x3fft);
            copy=real(wifft);
            //datawrite(copy,x2fft,x3fft,outf);
            for(i=0;i<x2;i++)
            {
                for(j=0;j<x3;j++)
                {
                    imagpower=swave[i][j]*copy(i,j);
                    recvA[j][i][k]+=imagpower;
                }
            }
        }
        //outf.close();
    }

    void cleardata()
    {
        matcopy(FA, 0.0, x3, x2, x1);
        matcopy(DA, 0.0, x3, x2, x1);
        if(recvid!=0)
        {
            matcopy(recvA, 0.0, x3, x2, recva_num);
        }
    }
};

//////////////////////////////////////////////////
class OF_2D
{
public:
    float **vx1=NULL, **vx2=NULL, **vz1=NULL, **vz2=NULL, ***wave=NULL;
    float **v_x=NULL, **v_z=NULL, **zz2d; 
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
        if(nt<5)
        {
            cout<<"Error: The data is not enough!"<<endl;
        }
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
        matcopy(vx2,0.0,x2,x3);
        matcopy(vz1,0.0,x2,x3);
        matcopy(vz2,0.0,x2,x3);
        matcopy(wave,0.0,x1,x2,x3);
        this->velocityInitialization();
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
    void CIG_OF_cal(T1 **wave, int N=5, int n=1)
    {
        int i;
        this->addtimeslicecal(wave);
        //this->velocityInitialization();
        matcopy(vx1,vx2,x2,x3);
        matcopy(vz1,vz2,x2,x3);
        for(i=0;i<N;i++)
        {this->velocityIteration();}
        matsmooth(this->vx2,this->vx1,this->x2,this->x3,n);
        matsmooth(this->vz2,this->vz1,this->x2,this->x3,n);
    }

    void CIG_OF_cal(int N=5, int n=1)
    {
        int i;
        matcopy(vx1,vx2,x2,x3);
        matcopy(vz1,vz2,x2,x3);
        for(i=0;i<N;i++)
        {this->velocityIteration();}
        matsmooth(this->vx2,this->vx1,this->x2,this->x3,n);
        matsmooth(this->vz2,this->vz1,this->x2,this->x3,n);
    }

    template<typename T1>
    void addtimeslicecal(T1 **p)
    {
        int k;
        matcopy(wave[0],p,x2,x3);
        zz2d=wave[0];
        for(k=0;k<x1-1;k++)
        {
            //matcopy(wave[k],wave[k+1],x2,x3);
            wave[k]=wave[k+1];
        }
        wave[x1-1]=zz2d;
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

////////////////////////////////////////////////////////
class Poynting_2D
{
public:
    float **vx=NULL, **vz=NULL, ***wave=NULL, **zz2d;
    float dx,dz,dt;
    float xs[5]={0.083333333,-0.666666667,0.0,0.666666667,-0.083333333};
    int x1,x2,x3;

    Poynting_2D()
    {
        x1=0;x2=0;x3=0;
        dx=5.0,dz=5.0,dt=0.0005;
        cout<<"Warning: Creat an Empty Object-optical-flow-2D"<<endl;
    }
    Poynting_2D(int nt, int nz, int nx)
    {
        dx=5.0, dz=5.0, dt=0.0005;
        x1=nt, x2=nz, x3=nx;
        vx=newfmat(x2,x3);
        vz=newfmat(x2,x3);
        wave=newfmat(x1,x2,x3);
        matcopy(vx,0.0,x2,x3);
        matcopy(vz,0.0,x2,x3);
        matcopy(wave,0.0,x1,x2,x3);
    }
    ~Poynting_2D()
    {
        matdelete(vx,x2);
        matdelete(vz,x2);
        matdelete(wave,x1,x2);
        vx=NULL;vz=NULL;wave=NULL;
    }

    template<typename T1>
    void addtimeslicecal(T1 **p)
    {
        int k;
        matcopy(wave[0],p,x2,x3);
        zz2d=wave[0];
        for(k=0;k<x1-1;k++)
        {
            wave[k]=wave[k+1];
            //matcopy(wave[k],wave[k+1],x2,x3);
        }
        wave[x1-1]=zz2d;
    }

    void velocityCalculate()
    {
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
                vx[i][j]=px*pt;
                vz[i][j]=pz*pt;
            }
        }
    }
    void dataClear()
    {
        matcopy(vx,0.0,x2,x3);
        matcopy(vz,0.0,x2,x3);
        matcopy(wave,0.0,x1,x2,x3);
    }
};

class direction_wavefield_state_transition
{
public:
    struct forA_state_transition
    {
        int nz,nx,na;
        float **velocity,**pz1,**px1,**p_nor,**pz2,**px2;
        float **coord_z1,**coord_z2,**coord_x1,**coord_x2;
        float **err_pz,**err_px,**err_nor,**wave1,**wave2;
        float dz,dx,t_gep,dt,pz_for_or_back_word,px_for_or_back_word;
        float ***angle,a_beg,a_gep;
    }par;

    void forA_state_transition_initialize_fmat(struct forA_state_transition & par);
    void A_state_transition(struct forA_state_transition & par);
    void pfmat_normalize(float **pz, float **px, int nz, int nx, float **p_nor);
    void p_normalize(float & pz, float & px, float & p, float yz=0.000001);
    void cleardata();


    direction_wavefield_state_transition()
    {
        ;
    }
    direction_wavefield_state_transition(int na, int nz, int nx,\
         float A_BEG=-90.0, float A_GEP=1.0)
    {
        par.nz=nz, par.nx=nx, par.na=na;
        par.a_beg=A_BEG,par.a_gep=A_GEP;
        this->forA_state_transition_initialize_fmat(this->par);
    }
    ~direction_wavefield_state_transition()
    {
        ;
    }
};

void direction_wavefield_state_transition::\
cleardata()
{
    matcopy(this->par.coord_x1,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.coord_x2,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.coord_z1,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.coord_z2,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.pz1,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.pz2,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.px1,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.px2,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.p_nor,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.err_pz,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.err_px,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.err_nor,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.wave1,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.wave2,0.0,this->par.nz,this->par.nx);
    matcopy(this->par.angle,0.0,this->par.nx,this->par.nz,this->par.na);
}

void direction_wavefield_state_transition::\
forA_state_transition_initialize_fmat(struct forA_state_transition & par)
{
    par.coord_x1=newfmat(par.nz,par.nx);
    par.coord_x2=newfmat(par.nz,par.nx);
    par.coord_z1=newfmat(par.nz,par.nx);
    par.coord_z2=newfmat(par.nz,par.nx);
    par.velocity=newfmat(par.nz,par.nx);
    par.pz1=newfmat(par.nz,par.nx);
    par.pz2=newfmat(par.nz,par.nx);
    par.px1=newfmat(par.nz,par.nx);
    par.px2=newfmat(par.nz,par.nx);
    par.p_nor=newfmat(par.nz,par.nx);
    par.err_pz=newfmat(par.nz,par.nx);
    par.err_px=newfmat(par.nz,par.nx);
    par.err_nor=newfmat(par.nz,par.nx);
    par.wave1=newfmat(par.nz,par.nx);
    par.wave2=newfmat(par.nz,par.nx);
    par.angle=newfmat(par.nx,par.nz,par.na);
}

void direction_wavefield_state_transition::\
A_state_transition(struct forA_state_transition & par)
{
    int i,j,k,i1,j1,ang;
    float xsz,xsx,theta,sx,sz,d,d1,d2,d3,d4;
    xsz=par.dt*par.t_gep*par.pz_for_or_back_word;
    xsx=par.dt*par.t_gep*par.px_for_or_back_word;
    //cout<<xsz<<"  "<<xsx<<endl;
    pfmat_normalize(par.pz1, par.px1, par.nz, par.nx,par.p_nor);
    for(i=0;i<par.nz;i++)
    {
        for(j=0;j<par.nx;j++)
        {
            par.coord_x1[i][j]=j*par.dx;
            par.coord_z1[i][j]=i*par.dz;
        } 
    }

    for(i=0;i<par.nz;i++)
    {
        for(j=0;j<par.nx;j++)
        {
            par.coord_x2[i][j]=par.coord_x1[i][j]+par.px1[i][j]*xsx*par.velocity[i][j];
            par.coord_z2[i][j]=par.coord_z1[i][j]+par.pz1[i][j]*xsz*par.velocity[i][j];
        }
    }

    for(i=0;i<par.nz;i++)
    {
        for(j=0;j<par.nx;j++)
        {
            j1=int(par.coord_x2[i][j]/par.dx);
            i1=int(par.coord_z2[i][j]/par.dz);

            sx=(par.coord_x2[i][j]/par.dx);
            sz=(par.coord_z2[i][j]/par.dz);
            d1=(sx-j1)*(sx-j1)+(sz-i1)*(sz-i1)+1e-12;
            d2=(sx-j1-1)*(sx-j1-1)+(sz-i1)*(sz-i1)+1e-12;
            d3=(sx-j1)*(sx-j1)+(sz-i1-1)*(sz-i1-1)+1e-12;
            d4=(sx-j1-1)*(sx-j1-1)+(sz-i1-1)*(sz-i1-1)+1e-12;
            d=0,d=1/d1+1/d2+1/d3+1/d4;
            d1=1/d1/d;d2=1/d2/d;d3=1/d3/d;d4=1/d4/d;
            //par.px1[i][j]=par.px1[i][j]*par.p_nor[i][j];
            //par.pz1[i][j]=par.pz1[i][j]*par.p_nor[i][j];
            if(i1>0 && i1<par.nz-1 && j1>0 && j1<par.nx-1)
            {
                par.wave2[i1][j1]+=par.wave1[i][j];
                par.px2[i1][j1]=par.px1[i][j];
                par.pz2[i1][j1]=par.pz1[i][j];
                par.px2[i1][j1+1]=par.px1[i][j];
                par.pz2[i1][j1+1]=par.pz1[i][j];
                par.px2[i1+1][j1]=par.px1[i][j];
                par.pz2[i1+1][j1]=par.pz1[i][j];
                par.px2[i1+1][j1+1]=par.px1[i][j];
                par.pz2[i1+1][j1+1]=par.pz1[i][j];

                theta=anglecal_z1_x0(par.pz1[i][j],par.px1[i][j]);
                ang=round((theta-par.a_beg)/par.a_gep);
                //cout<<ang<<" | "<<endl;
                if(ang>=0 && ang<par.na)
                {
                    par.angle[j1][i1][ang]+=par.wave1[i][j]*d1;
                    par.angle[j1][i1+1][ang]+=par.wave1[i][j]*d3;
                    par.angle[j1+1][i1][ang]+=par.wave1[i][j]*d2;
                    par.angle[j1+1][i1+1][ang]+=par.wave1[i][j]*d4;
                }
            }
        }
    }

}

void direction_wavefield_state_transition::\
pfmat_normalize(float **pz, float **px, int nz, int nx, float **p_nor)
{
    int i,j;
    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
        {
            p_normalize(pz[i][j],px[i][j],p_nor[i][j]);
        } 
    }
}

void direction_wavefield_state_transition::\
p_normalize(float & pz, float & px, float & p, float yz)
{
    p=0;
    p=p+pz*pz;
    p=p+px*px;
    p=sqrt(p);
    if(p>yz)
    {
        pz=pz/p;
        px=px/p;
    }
    else
        p=0;
}

///////////////////////////////////////////////////
fmat windowsZ(int Z, int X, float thetaN)
{
    int i, j;
    float xs, theta;
    fmat w(Z,X);
    w.fill(1.0);
    thetaN=(thetaN/360.0)*2.0*3.1415926;
    for(i=0;i<Z;i++)
    {
        for(j=0;j<X;j++)
        {
            if(i<Z/2 && j<X/2)
            {
                theta=abs(atan(float(j)/i));
                if(theta<=thetaN)
                {
                    xs=Blackman(theta,thetaN);
                    w(i,j)=w(i,j)*xs;
                }
            }
            else if(i<Z/2 && j>=X/2)
            {
                theta=abs(atan(float(X-1-j)/i));
                if(theta<=thetaN)
                {
                    xs=Blackman(theta,thetaN);
                    w(i,j)=w(i,j)*xs;
                }
            }
            else if(i>=Z/2 && j<X/2)
            {
                theta=abs(atan(float(j)/(Z-1-i)));
                if(theta<=thetaN)
                {
                    xs=Blackman(theta,thetaN);
                    w(i,j)=w(i,j)*xs;
                }
            }
            else if(i>=Z/2 && j>=X/2)
            {
                theta=abs(atan(float(X-1-j)/(Z-1-i)));
                if(theta<=thetaN)
                {
                    xs=Blackman(theta,thetaN);
                    w(i,j)=w(i,j)*xs;
                }
            }
        }
    }
    return w;
}

void getangle1(float **angle1, int Z, int X)
{
   int z0,x0,i,j;
   float pz,px;
   z0=Z/2;x0=X/2;
   for(i=0;i<Z;i++)
   {
      for(j=0;j<X;j++)
      {
         pz=i-z0;
         px=j-x0;
         angle1[i][j]=anglecal_z1_x0(px,pz);
      }
   }
   datawrite(angle1,Z,X,"getangle1.bin");
}

void getangle2(float **angle1, int Z, int X)
{
   int z0,x0,i,j;
   float pz,px;
   z0=Z/2;x0=X/2;
   for(i=0;i<Z;i++)
   {
      for(j=0;j<X;j++)
      {
         pz=i-z0;
         px=j-x0;
         angle1[i][j]=anglecal_z1_x0(pz,px);
      }
   }
   datawrite(angle1,Z,X,"getangle2.bin");
}

void fft2dwindows(float **ang, float **win, int x1, int x2, float theta, float thetaN)
{
   int i,j;
   float xs,dtheta;
   for(i=0;i<x1;i++)
   {
      for(j=0;j<x2;j++)
      {
         dtheta=abs(ang[i][j]-theta);
         if(dtheta<=thetaN)
         {
            dtheta=thetaN-dtheta;
            xs=Blackman(dtheta,thetaN);
            win[i][j]=xs;
         }
      }
   }

}

float angletransform(float a)
{
   float a2;
   if(a<=90 && a>=-180)
   {
      a2=-a-90.0;
   }
   else if(a>90 && a<=180)
   {
      a2=-a+270.0;
   }
   return a2;
}

fmat fft2dwintransform(float **win, int Z, int X)
{
   int i,j,x0,z0;
   fmat wintrans(Z,X);
   z0=Z/2;x0=X/2;
   for(i=0;i<Z;i++)
   {
      for(j=0;j<X;j++)
      {
         if(i<z0 && j<x0)
         {
            wintrans(i,j)=win[i+z0][j+x0];
         }
         if(i<z0 && j>=x0)
         {
            wintrans(i,j)=win[i+z0][j-x0];
         }
         if(i>=z0 && j<x0)
         {
            wintrans(i,j)=win[i-z0][j+x0];
         }
         if(i>=z0 && j>=x0)
         {
            wintrans(i,j)=win[i-z0][j-x0];
         }
      }
   }
   return wintrans;
}

void fft2dwindows_Gabor(float **ang, float **win, int x1, int x2, float theta, float thetaN)
{
   int i,j;
   float xs,dtheta;
   for(i=0;i<x1;i++)
   {
      for(j=0;j<x2;j++)
      {
         dtheta=abs(ang[i][j]-theta);
         if(true)
         {
            //dtheta=thetaN-dtheta;
            xs=Gabor(dtheta,thetaN);
            win[i][j]=xs;
         }
      }
   }
}

void fft2dwindows_Blackman(float **ang, float **win, int x1, int x2, float theta, float thetaN)
{
   int i,j;
   float xs,dtheta;
   for(i=0;i<x1;i++)
   {
      for(j=0;j<x2;j++)
      {
         dtheta=abs(ang[i][j]-theta);
         if(dtheta<=thetaN)
         {
            dtheta=thetaN-dtheta;
            xs=Blackman(dtheta,thetaN);
            win[i][j]=xs;
         }
      }
   }
}

float Gabor(float da,float a)
{
   float n,pi(3.1415926);
   n=(1.0/(2*sqrt(pi*a)))*exp(-da*da/4.0/a);
   return n;
}

#endif
