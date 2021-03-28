/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef ELASTIC2D_HPP
#define ELASTIC_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "../xjc.h"
/////////////////////////////////////////////////////////////////////////////////////////////
float** newfmatcs(int x1, int x2, float f)
{
    float **p,**p2;
    int j;
    p=new float*[x1];      
    for(j=0;j<x1;j++)  
        {  
        p[j]=new float[x2];
        }
    p2=p;
    p=NULL;
    matcopy(p2,f,x1,x2); 
    return p2;
}

struct elastic_wave_data
{
    float **vxx=NULL, **vxy=NULL, **vyx=NULL, **vyy=NULL; //PML boundary
    float **vxx2=NULL, **vxy2=NULL, **vyx2=NULL, **vyy2=NULL; //PML boundary
    float **txx=NULL, **txy=NULL, **tyy=NULL; //PML boundary
    float **txxx=NULL, **txyx=NULL, **tyyx=NULL; //PML boundary
    float **txxy=NULL, **txyy=NULL, **tyyy=NULL; //PML boundary
    float **txxx2=NULL, **txyx2=NULL, **tyyx2=NULL; //PML boundary
    float **txxy2=NULL, **txyy2=NULL, **tyyy2=NULL; //PML boundary
    float **vpx12=NULL, **vpx1=NULL;
    float **vpx22=NULL, **vpx2=NULL;
};

class elastic2D
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
    void cal(float** ut2 ,float **ut1,float **u, float **m,const char & label);

public:
    elastic_wave_data data;
    float **vp=NULL,**vs=NULL,**lmd=NULL,**uz=NULL,**ux=NULL,\
        **ro1=NULL,**miu=NULL,**ro=NULL,**mo=NULL,**mo1=NULL,\
        **Txy=NULL,**Tyy=NULL,**Txx=NULL,vp1,vs1,roo; //velocity model
    float dx,dy,dt,R,dt2,t3;
    int nx,ny,suface,PML_wide;
    
    elastic2D();
    elastic2D(int x, int y);
    ~elastic2D();

    //void timeslicecal(int restart=0);
    void timeslicecal();
    void cleardata();
    void updatepar();
};
///////////////////////////////////////////////////////////////////////////////////////////

elastic2D::elastic2D()
{
    nx=0;ny=0;
    dx=5.0;dy=5.0;dt=0.001;
    dt2=dt;
    PML_wide=30;suface=1;R=1000;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

elastic2D::elastic2D(int z, int x)
{
    nx=x;ny=z;t3=1;
    dx=5.0;dy=5.0;dt=0.001;dt2=dt;
    PML_wide=30;suface=1;R=1000;
    vp1=3000,vs1=2000,roo=2000;
    
    data.txx=newfmatcs(ny,nx,0.0),data.txxx2=newfmatcs(ny,nx,0.0),data.txxx=newfmatcs(ny,nx,0.0),\
    data.txxy2=newfmatcs(ny,nx,0.0),data.txxy=newfmatcs(ny,nx,0.0),data.txy=newfmatcs(ny,nx,0.0),\
    data.txyx2=newfmatcs(ny,nx,0.0),data.txyx=newfmatcs(ny,nx,0.0),data.txyy2=newfmatcs(ny,nx,0.0),\
    data.txyy=newfmatcs(ny,nx,0.0),data.tyy=newfmatcs(ny,nx,0.0),data.tyyx2=newfmatcs(ny,nx,0.0),\
    data.tyyx=newfmatcs(ny,nx,0.0),data.tyyy2=newfmatcs(ny,nx,0.0),data.tyyy=newfmatcs(ny,nx,0.0),\
    data.vxx2=newfmatcs(ny,nx,0.0),data.vxx=newfmatcs(ny,nx,0.0),data.vxy2=newfmatcs(ny,nx,0.0),\
    data.vxy=newfmatcs(ny,nx,0.0),data.vyx2=newfmatcs(ny,nx,0.0),data.vyx=newfmatcs(ny,nx,0.0),\
    data.vyy2=newfmatcs(ny,nx,0.0),data.vyy=newfmatcs(ny,nx,0.0);
    data.vpx12=newfmatcs(ny,nx,0.0),data.vpx1=newfmatcs(ny,nx,0.0);
    data.vpx22=newfmatcs(ny,nx,0.0),data.vpx2=newfmatcs(ny,nx,0.0);

    vs=newfmatcs(ny,nx,vs1),vp=newfmatcs(ny,nx,vp1),ro=newfmatcs(ny,nx,roo);
    miu=newfmatcs(ny,nx,vs1*vs1*roo),lmd=newfmatcs(ny,nx,(vp1*vp1*roo-2*vs1*vs1*roo));
    mo=newfmatcs(ny,nx,vp1*vp1*roo),ro1=newfmatcs(ny,nx,1.0/roo); 
    ux=newfmatcs(ny,nx,0.0),uz=newfmatcs(ny,nx,0.0);
    mo1=newfmatcs(ny,nx,(lmd[0][0]+2*miu[0][0])/(2*lmd[0][0]+2*miu[0][0])/ro[0][0]); 
    Txx=newfmatcs(ny,nx,0.0),Tyy=newfmatcs(ny,nx,0.0),Txy=newfmatcs(ny,nx,0.0);
}

elastic2D::~elastic2D()
{
    matdelete(vs,ny),matdelete(vp,ny),matdelete(lmd,ny),\
    matdelete(miu,ny),matdelete(ro,ny),matdelete(ro1,ny),matdelete(mo,ny);
    vp=NULL,vs=NULL,lmd=NULL,miu=NULL,ro=NULL,ro1=NULL,mo=NULL;

    matdelete(data.txx,ny),matdelete(data.txxx2,ny),matdelete(data.txxx,ny),\
    matdelete(data.txxy2,ny),matdelete(data.txxy,ny),matdelete(data.txy,ny),\
    matdelete(data.txyx2,ny),matdelete(data.txyx,ny),matdelete(data.txyy2,ny),\
    matdelete(data.txyy,ny),matdelete(data.tyy,ny),matdelete(data.tyyx2,ny),\
    matdelete(data.tyyx,ny),matdelete(data.tyyy2,ny),matdelete(data.tyyy,ny),\
    matdelete(data.vxx2,ny),matdelete(data.vxx,ny),matdelete(data.vxy2,ny),\
    matdelete(data.vxy,ny),matdelete(data.vyx2,ny),matdelete(data.vyx,ny),\
    matdelete(data.vyy2,ny),matdelete(data.vyy,ny);
    matdelete(ux,ny),matdelete(uz,ny);
    matdelete(data.vpx22,ny),matdelete(data.vpx1,ny);
    matdelete(data.vpx12,ny),matdelete(data.vpx2,ny);
    matdelete(mo1,ny); 
    matdelete(Txx,ny),matdelete(Tyy,ny),matdelete(Txy,ny);
    uz=NULL,ux=NULL,mo1=NULL;
    Txy=NULL,Txx=NULL,Tyy=NULL;
    data.vxx=NULL, data.vxy=NULL, data.vyx=NULL, data.vyy=NULL; //PML boundary
    data.vxx2=NULL, data.vxy2=NULL, data.vyx2=NULL, data.vyy2=NULL; //PML boundary
    data.txx=NULL, data.txy=NULL, data.tyy=NULL; //PML boundary
    data.txxx=NULL, data.txyx=NULL, data.tyyx=NULL; //PML boundary
    data.txxy=NULL, data.txyy=NULL, data.tyyy=NULL; //PML boundary
    data.txxx2=NULL, data.txyx2=NULL, data.tyyx2=NULL; //PML boundary
    data.txxy2=NULL, data.txyy2=NULL, data.tyyy2=NULL; //PML boundaryt
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

void elastic2D::cleardata()
{
    //this->setvelocity(0.0);
    matdelete(data.txx,ny),matdelete(data.txxx2,ny),matdelete(data.txxx,ny),\
    matdelete(data.txxy2,ny),matdelete(data.txxy,ny),matdelete(data.txy,ny),\
    matdelete(data.txyx2,ny),matdelete(data.txyx,ny),matdelete(data.txyy2,ny),\
    matdelete(data.txyy,ny),matdelete(data.tyy,ny),matdelete(data.tyyx2,ny),\
    matdelete(data.tyyx,ny),matdelete(data.tyyy2,ny),matdelete(data.tyyy,ny),\
    matdelete(data.vxx2,ny),matdelete(data.vxx,ny),matdelete(data.vxy2,ny),\
    matdelete(data.vxy,ny),matdelete(data.vyx2,ny),matdelete(data.vyx,ny),\
    matdelete(data.vyy2,ny),matdelete(data.vyy,ny);
    matdelete(ux,ny),matdelete(uz,ny);
    uz=NULL,ux=NULL;
    data.vxx=NULL, data.vxy=NULL, data.vyx=NULL, data.vyy=NULL; //PML boundary
    data.vxx2=NULL, data.vxy2=NULL, data.vyx2=NULL, data.vyy2=NULL; //PML boundary
    data.txx=NULL, data.txy=NULL, data.tyy=NULL; //PML boundary
    data.txxx=NULL, data.txyx=NULL, data.tyyx=NULL; //PML boundary
    data.txxy=NULL, data.txyy=NULL, data.tyyy=NULL; //PML boundary
    data.txxx2=NULL, data.txyx2=NULL, data.tyyx2=NULL; //PML boundary
    data.txxy2=NULL, data.txyy2=NULL, data.tyyy2=NULL; //PML boundaryt
    cout<<"All Matrix data has been clear!"<<endl;

    data.txx=newfmatcs(ny,nx,0.0),data.txxx2=newfmatcs(ny,nx,0.0),data.txxx=newfmatcs(ny,nx,0.0),\
    data.txxy2=newfmatcs(ny,nx,0.0),data.txxy=newfmatcs(ny,nx,0.0),data.txy=newfmatcs(ny,nx,0.0),\
    data.txyx2=newfmatcs(ny,nx,0.0),data.txyx=newfmatcs(ny,nx,0.0),data.txyy2=newfmatcs(ny,nx,0.0),\
    data.txyy=newfmatcs(ny,nx,0.0),data.tyy=newfmatcs(ny,nx,0.0),data.tyyx2=newfmatcs(ny,nx,0.0),\
    data.tyyx=newfmatcs(ny,nx,0.0),data.tyyy2=newfmatcs(ny,nx,0.0),data.tyyy=newfmatcs(ny,nx,0.0),\
    data.vxx2=newfmatcs(ny,nx,0.0),data.vxx=newfmatcs(ny,nx,0.0),data.vxy2=newfmatcs(ny,nx,0.0),\
    data.vxy=newfmatcs(ny,nx,0.0),data.vyx2=newfmatcs(ny,nx,0.0),data.vyx=newfmatcs(ny,nx,0.0),\
    data.vyy2=newfmatcs(ny,nx,0.0),data.vyy=newfmatcs(ny,nx,0.0);
}

void elastic2D::updatepar()
{
    int i,j;
    dt2=dt;
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            miu[i][j]=vs[i][j]*vs[i][j]*ro[i][j];
            lmd[i][j]=ro[i][j]*(vp[i][j]*vp[i][j]-2.0*vs[i][j]*vs[i][j]);
            mo[i][j]=ro[i][j]*vp[i][j]*vp[i][j];
            ro1[i][j]=1.0/ro[i][j];
        }
    }
}

void elastic2D::cal(float** ut2 ,float **ut1,float **u, float **m,const char & label)
{
    float DX,DY,DT,xshd;
    int X,Y,suface_PML;
    float fdx,fdy,fddx,fddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,i1,j1,n,t4;
    float du1(0),du2(0),du(0),dux(0),duy(0),duxy(0);
    float DT2,DT3,DX2,DY2,mo2;
    float C_X, C_Y;
    float *xs1_in=this->xs1,*xs2_in=this->xs2; 

    DX=this->dx,DY=this->dy,DT=this->dt2,xshd=this->PML_wide;
    X=this->nx,Y=this->ny,suface_PML=this->suface;
    DT2=DT*DT,DT3=DT2*DT,DX2=DX*DX,DY2=DY*DY;
    t5=float(DY)/DX;
    C_Y=log(R)*3/2/(xshd)/(xshd)/(xshd)/DY2/DY;
    C_X=log(R)*3/2/(xshd)/(xshd)/(xshd)/DX2/DX;

    if(label=='x')
    {
    for(j=5;j<Y-5;j++)
    {
        for(i=5;i<xshd;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=C_X*this->vp[j][i]*DX*(xshd-i)*DX*(xshd-i)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=xshd;i<X-xshd;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(i=X-xshd;i<X-5;i++)
        {  
            du=0;
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j][i+j1-5]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j][i+j1-4]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DX;
            du1=C_X*this->vp[j][i]*DX*(i-X+xshd)*DX*(i-X+xshd)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
    }
    }
    else if(label=='y')
    {
    for(i=5;i<X-5;i++)
    {
        for(j=5;j<xshd;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=C_Y*this->vp[j][i]*DY*(xshd-j)*DY*(xshd-j)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=xshd;j<Y-xshd;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=0;
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
        for(j=Y-xshd;j<Y-5;j++)
        {
            du=0;  
            for(j1=0;j1<5;j1++)
            {
                du=du+u[j+j1-5][i]*xs1_in[j1]*this->t3;
            }
            for(j1=5;j1<10;j1++)
            {

                du=du+u[j+j1-4][i]*xs1_in[j1]*this->t3;
            }
            du=m[j][i]*du/DY;
            du1=C_Y*this->vp[j][i]*DY*(j-Y+xshd)*DY*(j-Y+xshd)*ut1[j][i];
            ut2[j][i]=(du-du1)*DT+ut1[j][i];
        }
    }
    }
    else
    {
        cout<<"error!"<<endl;
    }
/*
    for(i=5;i<X-5;i++)
    {
        if(i>=0.5*(X))
        {t3=1;}
        else
        {t3=-1;}

        for(j=5;j<Y-5;j++)
        {
        //����ϵ����ö���ƫ΢�ֵ���ɢ����
            snx1=0.0;snx2=0.0;
            sny1=0.0;sny2=0.0;
            fdx=0;fddx=0;fdy=0;fddy=0;*/
/*
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
*/
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
/*
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
*/
            du1=0,du2=0,du=0,dux=0,duy=0,duxy=0;


//ע�ⲻҪֱ�Ӷ��徲̬��ָ��(�����)���Խ�����ַ�ķ�ʽ����ʱ��Ƭ,�ڶ��߳�������ܻ������,
//�Ʋ��ڶ��߳���,���ܻ��ڱ���կ���й��þ�̬�����Ĵ���,�Ӷ������ڴ�����
//���ڶ���̲��л᲻���������������δ������֤  (!!ע��!!)
    for(j=0;j<Y;j++)
    {
        for(i=0;i<X;i++)
        {
            ut1[j][i]=(ut2[j][i]);
        }
    }
}

void elastic2D::timeslicecal()
{
    int i,j;
    //cal(float** ut2 ,float **ut1,float **u, float **m,const char & label)
    this->cal(this->data.txyx2,this->data.txyx,this->uz,this->miu,'x');
    this->cal(this->data.txyy2,this->data.txyy,this->ux,this->miu,'y');
    this->cal(this->data.tyyx2,this->data.tyyx,this->ux,this->lmd,'x');
    this->cal(this->data.tyyy2,this->data.tyyy,this->uz,this->mo,'y');
    this->cal(this->data.txxx2,this->data.txxx,this->ux,this->mo,'x');
    this->cal(this->data.txxy2,this->data.txxy,this->uz,this->lmd,'y');
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->Txx[i][j]=this->data.txxx[i][j]+this->data.txxy[i][j];
            this->Tyy[i][j]=this->data.tyyx[i][j]+this->data.tyyy[i][j];
            this->Txy[i][j]=this->data.txyx[i][j]+this->data.txyy[i][j];
        }
    }
    this->cal(this->data.vyx2,this->data.vyx,this->Txy,this->ro1,'x');
    this->cal(this->data.vyy2,this->data.vyy,this->Tyy,this->ro1,'y');
    this->cal(this->data.vxx2,this->data.vxx,this->Txx,this->ro1,'x');
    this->cal(this->data.vxy2,this->data.vxy,this->Txy,this->ro1,'y');

    this->cal(this->data.vpx12,this->data.vpx1,this->Txx,this->mo1,'x');
    this->cal(this->data.vpx22,this->data.vpx2,this->Tyy,this->mo1,'x');

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            this->uz[i][j]=this->data.vyx[i][j]+this->data.vyy[i][j];
            this->ux[i][j]=this->data.vxx[i][j]+this->data.vxy[i][j];
        }
    }
    //this->t3=this->t3*(-1);
}

void addsouce(float **u, int n1, int n2, int n, float f,float dt,int nt,float d=1)
{
    int i,j;
    float g;
    for(i=n1-n;i<=n1+n;i++)
    {
        for(j=n2-n;j<=n2+n;j++)
        {
            g=wavelet02(nt,dt,f,1.0);        
            u[i][j]=u[i][j]+g*float(exp(-((i-n1)*(i-n1)+(j-n2)*(j-n2))/2.0/n/d)); 
        } 
    }
}

void elastic_test(int dmovie=1)
{
    int Z(500),X(500),T(3000);
    elastic2D A(Z,X);
    ofstream outf1,outf2,outf3,outf4;
    float **uu;
    uu=newfmatcs(Z,X,0.0);
    outf1.open("u1movie.bin");
    outf2.open("u2movie.bin");
    outf3.open("u3movie.bin");
    outf4.open("u4movie.bin");

    int k,i,j;
    for(k=0;k<T;k++)
    {
        /*
        A.ux[Z/2][X/2]=A.ux[Z/2][X/2]+f[k];
        A.ux[Z/2-1][X/2]=A.ux[Z/2-1][X/2]+f[k];
        A.ux[Z/2][X/2-1]=A.ux[Z/2][X/2-1]+f[k];
        A.ux[Z/2][X/2+1]=A.ux[Z/2][X/2+1]+f[k];
        A.ux[Z/2+1][X/2]=A.ux[Z/2+1][X/2]+f[k];
        */
        addsouce(A.ux,Z/2,X/2,5,30,A.dt,k,1);
        A.timeslicecal();

        if(k%dmovie==0)
        {
            matsmooth(uu,A.ux,Z,X,0);
            datawrite(uu, Z, X, outf1);
            matsmooth(uu,A.uz,Z,X,0);
            datawrite(uu, Z, X, outf2);
            for(i=0;i<Z;i++)
            {
                for(j=0;j<X;j++)
                {
                    uu[i][j]=A.data.vpx1[i][j]+A.data.vpx2[i][j];
                }
            }
            datawrite(uu, Z, X, outf3);
            //datawrite(A.data.vyy, Z, X, outf4);
        }
        if(k%100==0)
        {cout<<"now is running : "<<k<<endl;}
    }
    outf1.close();
    outf2.close();
    outf3.close();
    cout<<"Have output test movie."<<endl;
}

#endif

