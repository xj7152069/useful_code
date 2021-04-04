/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef WAVE2D_HPP
#define WAVE2D_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "../xjc.h"
/////////////////////////////////////////////////////////////////////////////////////////////
void wave2D_test(int Z, int X);
void wave2D_stabletest(int Z, int X, int T, float dz, float dx, float dt, float **v, int sz, int sx, float sf, float pmlwide, int dmovie);
float wavelet01(int k, float DT, float hz=30.0,float det2=1.0);
template <typename T1> void wavelet01(T1 *w, int N, float DT, float hz=30.0,float det2=1.0);
float wavelet02(int k, float DT, float hz=30.0,float det2=1.0);
template <typename T1> void wavelet02(T1 *w, int N, float DT, float hz=30.0,float det2=1.0);
template <typename T1> float* hilbert1D(T1 *s, int n, float dt);
float Blackman(float n, float N);

class wave2D
{
private:
//����ƫ΢����ɢ���ϵ��
    //float xs2[5]={1.666667/1.463612\
    12,-0.238095/1.463612,0.039683/1.463612,\
    -0.004960/1.463612,0.000317/1.463612};
    float xs2[5]={1.66666639,-0.238095194,3.96825150E-02,-4.96031437E-03,3.17460130E-04};
//һ��ƫ΢����ɢ���ϵ��
    float xs1[10]={-7.93651154E-04,9.92063992E-03,-5.95238246E-02,0.238095269,-0.833333373,\
0.833333373,-0.238095269,5.95238246E-02,-9.92063992E-03,7.93651154E-04};
    float **sx11=NULL, **sx12=NULL, **sx13=NULL; //PML boundary
    float **sxp21=NULL, **sxp22=NULL, **sxp23=NULL; //PML boundary
    float **sx21=NULL, **sx22=NULL; //PML boundary
    float **sx31=NULL,  **sx32=NULL, **sx33=NULL; //PML boundary

public:
    float **p2=NULL; //velocity model
    float ** & model=p2;
    float **s1=NULL, **s2=NULL, **s3=NULL; //time slices, add source to "s2"
    float dx,dy,dt,R;
    int nx,ny,suface,PML_wide;
    
    wave2D();
    wave2D(int x, int y);
    ~wave2D();

    void setvelocity(float v=3000.0);
    void timeslicecal(int restart=0);
    void timeslicecopy();
    void cleardata();
};
///////////////////////////////////////////////////////////////////////////////////////////

wave2D::wave2D()
{
    nx=0;ny=0;
    dx=5.0;dy=5.0;dt=0.0005;
    PML_wide=30;suface=1;R=1000;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

wave2D::wave2D(int z, int x)
{
    nx=x;ny=z;
    dx=5.0;dy=5.0;dt=0.0005;
    PML_wide=30;suface=1;R=1000;
    int i,j;
    s1=new float*[ny];s2=new float*[ny];s3=new float*[ny];  
    sx11=new float*[ny];sx12=new float*[ny];sx13=new float*[ny];    
    sxp21=new float*[ny];sxp22=new float*[ny];sxp23=new float*[ny];    
    sx21=new float*[ny];sx22=new float*[ny];     
    sx31=new float*[ny];sx32=new float*[ny];sx33=new float*[ny];

    for(j=0;j<ny;j++)  
        {  
        s1[j]=new float[nx];s2[j]=new float[nx];s3[j]=new float[nx];   
        sx11[j]=new float[nx];sx12[j]=new float[nx];sx13[j]=new float[nx];    
        sxp21[j]=new float[nx];sxp22[j]=new float[nx];sxp23[j]=new float[nx];  
        sx21[j]=new float[nx];sx22[j]=new float[nx];
        sx31[j]=new float[nx];sx32[j]=new float[nx];sx33[j]=new float[nx];   
        }

    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            s1[j][i]=0.0;s2[j][i]=0.0;s3[j][i]=0.0;
            sx11[j][i]=0.0;sx12[j][i]=0.0;sx13[j][i]=0.0;
            sxp21[j][i]=0.0;sxp22[j][i]=0.0;sxp23[j][i]=0.0;
            sx21[j][i]=0.0;sx22[j][i]=0.0;
            sx31[j][i]=0.0;sx32[j][i]=0.0;sx33[j][i]=0.0;
            }
        }

    p2=new float*[ny];      
    for(j=0;j<ny;j++)  
        {  
        p2[j]=new float[nx];
        } 

    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            p2[j][i]=3000;
            }
        }
}

wave2D::~wave2D()
{
    int i,j;

    for(int i=0;i<ny;i++)  
        delete []p2[i]; 
    delete []p2;

    for(int i=0;i<ny;i++)  
       {
        delete []s1[i];delete []s2[i];delete []s3[i];
        delete []sx11[i];delete []sx12[i];delete []sx13[i];
        delete []sxp21[i];delete []sxp22[i];delete []sxp23[i];
        delete []sx21[i];delete []sx22[i];
        delete []sx31[i];delete []sx32[i];delete []sx33[i];
       }
    delete []s1;delete []s2;delete []s3;
    delete []sx11;delete []sx12;delete []sx13;
    delete []sxp21;delete []sxp22;delete []sxp23;
    delete []sx21;delete []sx22;
    delete []sx31;delete []sx32;delete []sx33;
    p2=NULL;s1=NULL;s2=NULL;s3=NULL;
    sx11=NULL;sx12=NULL;sx13=NULL;
    sxp21=NULL;sxp22=NULL;sxp23=NULL;
    sx21=NULL;sx22=NULL;
    sx31=NULL;sx32=NULL;sx33=NULL;
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void wave2D::setvelocity(float v)
{
    int i,j;
    for(i=0;i<ny;i++)
        {
        for(j=0;j<nx;j++)
            {
            this->p2[i][j]=v;
            }
        }
}

void wave2D::cleardata()
{
    //this->setvelocity(0.0);
    int i,j;
    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            s1[j][i]=0.0;s2[j][i]=0.0;s3[j][i]=0.0;
            sx11[j][i]=0.0;sx12[j][i]=0.0;sx13[j][i]=0.0;
            sxp21[j][i]=0.0;sxp22[j][i]=0.0;sxp23[j][i]=0.0;
            sx21[j][i]=0.0;sx22[j][i]=0.0;
            sx31[j][i]=0.0;sx32[j][i]=0.0;sx33[j][i]=0.0;
            }
        }
    cout<<"All Matrix data has been clear!"<<endl;
}

void wave2D::timeslicecal(int restart)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    float fdx,fdy,fddx,fddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,n,t3,t4;
    float u1(0),u2(0),u(0),ux(0),uy(0);
    float DT2,DT3,DX2,DY2,mo2;
//PML�߽�����պ���d(x),�䳣��ϵ������
    float C_X, C_Y;

    float **sx11_in=this->sx11, **sx12_in=this->sx12, **sx13_in=this->sx13; //PML boundary
    float **sxp21i=this->sxp21, **sxp22i=this->sxp22, **sxp23i=this->sxp23; //PML boundary
    float **sx21_in=this->sx21, **sx22_in=this->sx22; //PML boundary
    float **sx31_in=this->sx31, **sx32_in=this->sx32, **sx33_in=this->sx33; //PML boundary
    float **p2_in=this->p2, **swap,*xs1_in=this->xs1,*xs2_in=this->xs2; //velocity model and swap
    float **s1_in=this->s1, **s2_in=this->s2, **s3_in=this->s3; //time slices, add source to "s2"

    DX=this->dx,DY=this->dy,DT=this->dt,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    DT2=DT*DT,DT3=DT2*DT,DX2=DX*DX,DY2=DY*DY;
    t5=float(Y)/X;
    C_Y=log(R)*3/2/(xshd)/(xshd)/(xshd)/DY2/DY;
    C_X=log(R)*3/2/(xshd)/(xshd)/(xshd)/DX2/DX;


    for(i=5;i<X-5;i++)
        {

        if(i>=0.5*(X))
            {t3=1;}
        else
            {t3=-1;}

        for(j=5;j<Y-5;j++)
            {
        //����ϵ����ö���ƫ΢�ֵ���ɢ����
            for(n=0;n<5;n++)
                {  
                u=u+2*xs2_in[n];
                u1=u1+xs2_in[n]*(s2_in[j-n-1][i]+s2_in[j+n+1][i]);
                u2=u2+xs2_in[n]*(s2_in[j][i-n-1]+s2_in[j][i+n+1]);
                }

            snx1=0.0;snx2=0.0;
            sny1=0.0;sny2=0.0;
            fdx=0;fddx=0;fdy=0;fddy=0;

            if(suface_PML==1)
                {
                if(i>=X-xshd-5 && j<t5*i && j>-t5*i+Y)	
                    snx1=i-(X-xshd-5);
                if(i<=xshd+5 && j>t5*i && j<-t5*i+Y)		
                    snx2=xshd+5-i;
                if(j>=Y-xshd-5 && j>=t5*i && j>=-t5*i+Y)		
                    sny1=j-(Y-xshd-5);
                if(j<=xshd+5 && j<=t5*i && j<=-t5*i+Y)		
                    sny2=xshd+5-j;		
                }  //��������
            else
                {
                if(i>=X-xshd-5 && j<=t5*i)	
                    snx1=i-(X-xshd-5);
                if(i<=xshd+5 && j<=-t5*i+Y)		
                    snx2=xshd+5-i;
                if(j>=Y-xshd-5 && j>t5*i && j>-t5*i+Y)		
                    sny1=j-(Y-xshd-5);
                if(j<=xshd+5)		
                    sny2=0;		
                }

            if(sny1 !=0)
                {
                fdy=p2_in[j][i]*C_Y*sny1*sny1*DY2;
                fddy=p2_in[j][i]*C_Y*2*sny1*DY;	
                }
            if(sny2 !=0)
                {
                fdy=p2_in[j][i]*C_Y*sny2*sny2*DY2;
                fddy=p2_in[j][i]*C_Y*2*sny2*DY;	
                }
            if(snx1 !=0)
                {
                fdx=p2_in[j][i]*C_X*snx1*snx1*DX2;
                fddx=p2_in[j][i]*C_X*2*snx1*DX;	
                }
            if(snx2 !=0)
                {
                fdx=p2_in[j][i]*C_X*snx2*snx2*DX2;
                fddx=p2_in[j][i]*C_X*2*snx2*DX;	
                }

            if(j>=0.5*(Y))
                {t4=1;}
            else
                {t4=-1;}

            mo2=p2_in[j][i]*p2_in[j][i];

        /*****************the code from wcl fortran***************************
			!equation 1
			uz_a3(iiz, iix) = 2.0*uz_a2(iiz, iix) - uz_a1(iiz, iix) + &
				dt2*(v2*z2_deri - 2.0*a*(uz_a2(iiz, iix) - uz_a1(iiz, iix))/dt - a*a*uz_a2(iiz, iix)) 
			!equation 2
			uz_bp3(iiz, iix) = 2.0*uz_bp2(iiz, iix) - uz_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*z1_deri - 2.0*a*(uz_bp2(iiz, iix) - uz_bp1(iiz, iix))/dt - a*a*uz_bp2(iiz, iix)) 
			uz_b2(iiz, iix) = uz_b1(iiz, iix) + dt*(uz_bp2(iiz, iix) - a*uz_b1(iiz, iix))
			!equation 3
			uz_c3(iiz, iix) = 2.0*uz_c2(iiz, iix) - uz_c1(iiz, iix) + dt2*v2*x2_deri
			!update
			u3(inz, inx) = uz_a3(iiz, iix) + uz_b2(iiz, iix) + uz_c3(iiz, iix)
            */

            if(snx1!=0 || snx2!=0)
                {
            //����ϵ�����һ��ƫ΢�ֵ���ɢ����
                for(n=0;n<10;n++)
                    {
                    if(n<5)
                        {
                        ux=ux+s2_in[j][i+n-5]*xs1_in[n];
                        uy=uy+s2_in[j+n-5][i]*xs1_in[n];
                        }
                    else
                        {
                        ux=ux+s2_in[j][i+n-4]*xs1_in[n];
                        uy=uy+s2_in[j+n-4][i]*xs1_in[n];
                        }
                    }

                //equation 1
                sx11_in[j][i]=mo2*DT2*(u2-u*s2_in[j][i])*(1.0/(DX2))\
                -fdx*fdx*DT2*sx12_in[j][i]+(2*sx12_in[j][i]\
                -sx13_in[j][i])+DT*(2*fdx*(sx13_in[j][i]-sx12_in[j][i]));

                //equation 2
			    sxp21i[j][i] = 2.0*sxp22i[j][i] - sxp23i[j][i]\
                +DT2*(-1.0*mo2*fddx*(1.0/(DX))*(ux*t3) - 2.0*fdx*(sxp22i[j][i] \
                - sxp23i[j][i])/DT - fdx*fdx*sxp22i[j][i]);
			    sx21_in[j][i] = sx22_in[j][i] + DT*(sxp22i[j][i] - fdx*sx22_in[j][i]) ;

                //equation 3
                sx31_in[j][i]=DT2*mo2*(1.0/(DY2))*(u1-u*s2_in[j][i])\
                +2*sx32_in[j][i]-sx33_in[j][i];
                }
            else if(sny1!=0 || sny2!=0)
                {
            //����ϵ�����һ��ƫ΢�ֵ���ɢ����
                for(n=0;n<10;n++)
                    {
                    if(n<5)
                        {
                        ux=ux+s2_in[j][i+n-5]*xs1_in[n];
                        uy=uy+s2_in[j+n-5][i]*xs1_in[n];
                        }
                    else
                        {
                        ux=ux+s2_in[j][i+n-4]*xs1_in[n];
                        uy=uy+s2_in[j+n-4][i]*xs1_in[n];
                        }
                    }

                //equation 1
                sx11_in[j][i]=mo2*DT2*(u1-u*s2_in[j][i])*(1.0/(DY2))\
                -fdy*fdy*DT2*sx12_in[j][i]+(2*sx12_in[j][i]\
                -sx13_in[j][i])+DT*(2*fdy*(sx13_in[j][i]-sx12_in[j][i]));

                //equation 2 : ��������ƫ΢��,�轫����Ϊһ��ƫ΢��(p)�Ķ��׵�����ɢ���
                sxp21i[j][i] = 2.0*sxp22i[j][i] - sxp23i[j][i]\
                +DT2*(-1.0*mo2*fddy*(1.0/(DY))*(uy*t4) - 2.0*fdy*(sxp22i[j][i] \
                - sxp23i[j][i])/DT - fdy*fdy*sxp22i[j][i]);
                sx21_in[j][i] = sx22_in[j][i] + DT*(sxp22i[j][i] - fdy*sx22_in[j][i]) ;

                //equation 3
                sx31_in[j][i]=DT2*mo2*(1.0/(DX2))*(u2-u*s2_in[j][i])\
                +2*sx32_in[j][i]-sx33_in[j][i];
                }

            if(snx1==0 && snx2==0 && sny1==0 && sny2==0)
            {
                s3_in[j][i]=(mo2*DT2*(u2-u*s2_in[j][i])*(1.0/(DX2))\
                +DT2*mo2*(1.0/(DY2))*(u1-u*s2_in[j][i]))+2*s2_in[j][i]-s1_in[j][i];
            }
            else
            {
                s3_in[j][i]=sx11_in[j][i]+sx21_in[j][i]+sx31_in[j][i];
            }
            u1=0,u2=0,u=0,ux=0,uy=0;
            }
        }

//ע�ⲻҪֱ�Ӷ��徲̬��ָ��(�����)���Խ�����ַ�ķ�ʽ����ʱ��Ƭ,�ڶ��߳�������ܻ������,
//�Ʋ��ڶ��߳���,���ܻ��ڱ���կ���й��þ�̬�����Ĵ���,�Ӷ������ڴ�����
//���ڶ���̲��л᲻���������������δ������֤  (!!ע��!!)
    for(j=0;j<Y;j++)
    {
    for(i=0;i<X;i++)
        {
        s1_in[j][i]=s2_in[j][i];s2_in[j][i]=s3_in[j][i];
        sx13_in[j][i]=sx12_in[j][i],sx12_in[j][i]=sx11_in[j][i];
        sx33_in[j][i]=sx32_in[j][i],sx32_in[j][i]=sx31_in[j][i];
        sx22[j][i]=sx21[j][i];
        sxp23i[j][i]=sxp22i[j][i],sxp22i[j][i]=sxp21i[j][i];
        }
    }

}

void wave2D::timeslicecopy()
{
    ;
}

void wave2D_test(int Z, int X)
{
    wave2D A(Z,X);
    ofstream outf1;
    outf1.open("testmovie.bin");
    int T(3000);

    int k;
    float *f;
    f=new float[T];
    wavelet01(f,T,A.dt);
    for(k=0;k<T;k++)
    {
        A.s2[Z/2][X/2]=A.s2[Z/2][X/2]+f[k];
        A.timeslicecal();
        {
        datawrite(A.s2, Z, X, outf1);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    cout<<"Have output test movie."<<endl;
}


void wave2D_stabletest(int Z, int X, int T, float dz, float dx, float dt, float **v, int sz, int sx, float sf, float pmlwide, int dmovie)
{
    wave2D A(Z,X);
    ofstream outf1;
    outf1.open("testmovie.bin");
    A.dx=dx,A.dy=dz,A.dt=dt;
    matcopy(A.p2,v,Z,X);
    A.PML_wide=pmlwide;

    int k;
    float *f;
    f=new float[T];
    wavelet01(f,T,dt,sf);
    for(k=0;k<T;k++)
    {
        A.s2[sz][sx]=A.s2[sz][sx]+f[k];
        A.timeslicecal();
        if(k%dmovie==0)
        {
        datawrite(A.s2, Z, X, outf1);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    cout<<"Have output test movie."<<endl;
}
///////////////////////////////////////////////////////////

float wavelet02(int k, float DT, float hz,float det2)
{
    float pi(3.1415926);
	float f,det;
    det=det2*0.05*(30.0/hz);

    f=(pi)*(pi)*hz*hz*(k*DT-det)*\
    exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
    *(3.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));

	return f;
}

template<typename T1>
void wavelet02(T1 *w, int N, float DT, float hz,float det2)
{
    float pi(3.1415926);
	float f,det;
    det=det2*0.05*(30.0/hz);
    int k;

    for(k=0;k<N;k++)
    {
        f=(pi)*(pi)*hz*hz*(k*DT-det)*\
        exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
        *(3.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));
        w[k]=f;
    }
}

float wavelet01(int k, float DT, float hz,float det2)
{
    float pi(3.1415926);
	float f,det;
    det=det2*0.05*(30.0/hz);

    f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
    *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));

	return f;
}

template<typename T1>
void wavelet01(T1 *w, int N, float DT, float hz,float det2)
{
    float pi(3.1415926);
	float f,det;
    det=det2*0.05*(30.0/hz);
    int k;

    for(k=0;k<N;k++)
    {
        f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
        *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));
        w[k]=f;
    }
}

template<typename T1>
float* hilbert1D(T1 *s, int n, float dt)
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
            s2[i]=s[i];
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
    return sh2;
}

float Blackman(float n, float N)
{
    float xs, pi(3.1415926);
    xs=0.42-0.5*cos(2*pi*n/N/2)+0.08*cos(4*pi*n/N/2);
    return xs;
}

#endif

