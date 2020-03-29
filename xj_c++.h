/*********(version 1.0)***********/
/*
xj_c++.h
    c++ head file: some useful function
*/
/********************************/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
using namespace std;
using namespace arma;

class wave_modeling_2D
{
public:
    float **p2=NULL; //velocity model
    float **s1=NULL, **s2=NULL, **s3=NULL; //time pices
    float **sx11=NULL, **sx12=NULL, **sx13=NULL; //PML boundary
    float **sx21=NULL, **sx22=NULL, **sx23=NULL, **sx24=NULL; //PML boundary
    float **sx31=NULL,  **sx32=NULL, **sx33=NULL; //PML boundary

    float xs2[5]={1.666667,-0.238095,0.039683,-0.004960,0.000317};
    float xs1[10]={-0.0007926,0.00991800,-0.0595200,0.238080,-0.833333,\
0.833333,-0.238080,0.0595200,-0.00991800,0.0007926};

    float dx,dy,dt,PML_wide;
    int nx,ny,suface;

    wave_modeling_2D()
    {
    nx=0;
    ny=0;
    cout<<"Warning: you creat a class-wave_modeling_2D without initialization"<<endl;
    }

    wave_modeling_2D(int x, int y)
    {
        nx=x;
        ny=y;
        dx=10.0;dy=10.0;dt=0.001;PML_wide=25;suface=1;
        int i,j;
        s1=new float*[ny+20];  
        s2=new float*[ny+20];   
        s3=new float*[ny+20];  
        sx11=new float*[ny+20];    
        sx12=new float*[ny+20];     
        sx13=new float*[ny+20];     
        sx21=new float*[ny+20];     
        sx22=new float*[ny+20];     
        sx23=new float*[ny+20];     
        sx24=new float*[ny+20];     
        sx31=new float*[ny+20];     
        sx32=new float*[ny+20];     
        sx33=new float*[ny+20];      
        for(j=0;j<ny+20;j++)  
            {  
            s1[j]=new float[nx+20];  
            s2[j]=new float[nx+20];   
            s3[j]=new float[nx+20];   
            sx11[j]=new float[nx+20];   
            sx12[j]=new float[nx+20];   
            sx13[j]=new float[nx+20];    
            sx21[j]=new float[nx+20];   
            sx22[j]=new float[nx+20];   
            sx23[j]=new float[nx+20];   
            sx24[j]=new float[nx+20];   
            sx31[j]=new float[nx+20];   
            sx32[j]=new float[nx+20];    
            sx33[j]=new float[nx+20];   
            }

        for(j=0;j<ny+20;j++)
            {
            for(i=0;i<nx+20;i++)
                {
                s1[j][i]=0.0;s2[j][i]=0.0;s3[j][i]=0.0;
                sx11[j][i]=0.0;sx12[j][i]=0.0;sx13[j][i]=0.0;
                sx21[j][i]=0.0;sx22[j][i]=0.0;sx23[j][i]=0.0;sx24[j][i]=0.0;
                sx31[j][i]=0.0;sx32[j][i]=0.0;sx33[j][i]=0.0;
                }
            }

        p2=new float*[ny];      
        for(j=0;j<y;j++)  
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

    ~wave_modeling_2D()
    {
        int i,j;

        for(int i=0;i<ny;i++)  
            delete []p2[i]; 
        delete []p2;

        for(int i=0;i<ny+20;i++)  
           {delete []s1[i];
            delete []s2[i];
            delete []s3[i];
            delete []sx11[i];
            delete []sx12[i];
            delete []sx13[i];
            delete []sx21[i];
            delete []sx22[i];
            delete []sx23[i];
            delete []sx24[i];
            delete []sx31[i];
            delete []sx32[i];
            delete []sx33[i];
           }//先单独释放第一维中每个数组的内存  
        delete []s1;
        delete []s2;
        delete []s3;
        delete []sx11;
        delete []sx12;
        delete []sx13;
        delete []sx21;
        delete []sx22;
        delete []sx23;
        delete []sx24;
        delete []sx31;
        delete []sx32;
        delete []sx33;
        s1=NULL;
        s2=NULL;
        s3=NULL;
        sx11=NULL;
        sx12=NULL;
        sx13=NULL;
        sx21=NULL;
        sx22=NULL;
        sx23=NULL;
        sx24=NULL;
        sx31=NULL;
        sx32=NULL;
        sx33=NULL;
        cout<<"You delete a class-wave_modeling_2D"<<endl;
    }
/////////////function wave modeling 2D Start/////////////
    void wave_modeling_2D_time_pice()
    {
        float DX,DY,DT,xshd;
        int X,Y,suface_PML;
        DX=dx;
        DY=dy;
        DT=dt;
        xshd=PML_wide;
        X=nx;
        Y=ny;
        suface_PML=suface;
/*
        float **p2=this->p2; //velocity model
        float **s1=this->s1, **s2=this->s2, **s3=this->s3; //time pices
        float **sx11=this->sx11, **sx12=this->sx12, **sx13=this->sx13; //PML boundary
        float **sx21=this->sx21, **sx22=this->sx22, **sx23=this->sx23, **sx24=this->sx24; //PML boundary
        float **sx31=this->sx31,  **sx32=this->sx32, **sx33=this->sx33; //PML boundary
        float *xs2=this->xs2;
        float *xs1=this->xs1;
*/
        
        float dx,dy,ddx,ddy,snx1,sny1,snx2,sny2,t2,t5;
        int i,j,n,i1,j1,t3,t4;
        float u1(0),u2(0),u(0),ux(0),uy(0);
        float C3=3/2/(xshd)/DX*log(10000000)/(xshd)/DX/(xshd)/DX;
        t5=float(Y)/X;

        for(i=5;i<X+5;i++)
            {
            for(j=5;j<Y+5;j++)
                {

                for(n=0;n<5;n++)
                    {  
                    u=u+2*xs2[n];
                    u1=u1+xs2[n]*(s2[j-n-1][i]+s2[j+n+1][i]);
                    u2=u2+xs2[n]*(s2[j][i-n-1]+s2[j][i+n+1]);
                    }

                for(n=0;n<10;n++)
                    {
                    if(n<5)
                        {
                        ux=ux+s2[j][i+n-5]*xs1[n];
                        uy=uy+s2[j+n-5][i]*xs1[n];
                        }
                    else
                        {
                        ux=ux+s2[j][i+n-4]*xs1[n];
                        uy=uy+s2[j+n-4][i]*xs1[n];
                        }
                    }

                snx1=0.0;snx2=0.0;
                sny1=0.0;sny2=0.0;
                dx=0;ddx=0;dy=0;ddy=0;

                if(suface_PML==1)
                    {
                    if(i>=X-xshd+5 && j<=t5*i+5*(1-t5) && j>=-t5*i+Y+5+5*t5)	
                        snx1=i-(X-xshd+5);
                    if(i<=xshd+5 && j>=t5*i+5*(1-t5) && j<=-t5*i+Y+5+t5*5)		
                        snx2=xshd+5-i;
                    if(j>=Y-xshd+5 && j>=t5*i+5*(1-t5) && j>=-t5*i+Y+5+5*t5)		
                        sny1=j-(Y-xshd+5);
                    if(j<=xshd+5 && j<=t5*i+5*(1-t5) && j<=-t5*i+Y+5+5*t5)		
                        sny2=xshd+5-j;		
                    }  //角落处理问题
                else if(suface_PML==0)
                    {
                    if(i>=X-xshd+5 && j<=t5*i+5*(1-t5))	
                        snx1=i-(X-xshd+5);
                    if(i<=xshd+5 && j<=-t5*i+Y+5+t5*5)		
                        snx2=xshd+5-i;
                    if(j>=Y-xshd+5 && j>=t5*i+5*(1-t5) && j>=-t5*i+Y+5+5*t5)		
                        sny1=j-(Y-xshd+5);
                    if(j<=xshd+5)		
                        sny2=0;		
                    }

                if(sny1 !=0)
                    {
                    dy=p2[j-5][i-5]*C3*sny1*sny1*DY*DY;
                    ddy=p2[j-5][i-5]*C3*2*sny1*DY;	
                    }
                if(sny2 !=0)
                    {
                    dy=p2[j-5][i-5]*C3*sny2*sny2*DY*DY;
                    ddy=p2[j-5][i-5]*C3*2*sny2*DY;	
                    }
                if(snx1 !=0)
                    {
                    dx=p2[j-5][i-5]*C3*snx1*snx1*DX*DX;
                    ddx=p2[j-5][i-5]*C3*2*snx1*DX;	
                    }
                if(snx2 !=0)
                    {
                    dx=p2[j-5][i-5]*C3*snx2*snx2*DX*DX;
                    ddx=p2[j-5][i-5]*C3*2*snx2*DX;	
                    }

                if(i>=0.5*(X+5))
                    t3=1;
                else
                    t3=1;

                if(j>=0.5*(Y+5))
                    t4=1;
                else
                    t4=1;

                if(snx1!=0 || snx2!=0)
                    {
                    sx11[j][i]=p2[j-5][i-5]*p2[j-5][i-5]\
                    *DT*DT*(u2-u*s2[j][i])*(1.0/(DX*DX))\
                    -dx*dx*DT*DT*sx12[j][i]+(2*sx12[j][i]\
                    -sx13[j][i])+DT*(2*dx*(sx13[j][i]-sx12[j][i]));

                    sx21[j][i]=(-p2[j-5][i-5]*p2[j-5][i-5]\
                    *ddx*(1.0/(DX))*(ux*t3)-dx*dx*dx*sx22[j][i]\
                    +3*dx*dx*(sx23[j][i]-sx22[j][i])/DT+3*dx*\
                    (2*sx23[j][i]-sx22[j][i]-sx24[j][i])/(DT*DT)\
                    +(3*sx22[j][i]-3*sx23[j][i]+sx24[j][i])\
                    /(DT*DT*DT))*(DT*DT*DT);

                    sx31[j][i]=DT*DT*p2[j-5][i-5]*p2[j-5][i-5]\
                    *(1.0/(DX*DX))*(u1-u*s2[j][i])+2*sx32[j][i]\
                    -sx33[j][i];
                    }
                else
                    {
                    sx11[j][i]=p2[j-5][i-5]*p2[j-5][i-5]\
                    *DT*DT*(u1-u*s2[j][i])*(1.0/(DX*DX))\
                    -dy*dy*DT*DT*sx12[j][i]+(2*sx12[j][i]\
                    -sx13[j][i])+DT*(2*dy*(sx13[j][i]-sx12[j][i]));

                    sx21[j][i]=(-p2[j-5][i-5]*p2[j-5][i-5]\
                    *ddy*(1.0/(DX))*(uy*t4)-dy*dy*dy*sx22[j][i]\
                    +3*dy*dy*(sx23[j][i]-sx22[j][i])/DT\
                    +3*dy*(2*sx23[j][i]-sx22[j][i]-sx24[j][i])\
                    /(DT*DT)+(3*sx22[j][i]-3*sx23[j][i]+sx24[j][i])\
                    /(DT*DT*DT))*(DT*DT*DT);

                    sx31[j][i]=DT*DT*p2[j-5][i-5]*p2[j-5][i-5]\
                    *(1.0/(DX*DX))*(u2-u*s2[j][i])+2*sx32[j][i]\
                    -sx33[j][i];
                    }   

                // if(snx1!=0 || snx2!=0 || sny1!=0 || sny2!=0)
                s3[j][i]=sx11[j][i]+sx21[j][i]+sx31[j][i];
                u1=0,u2=0,u=0,ux=0,uy=0;

                }
            }

        for(j1=0;j1<Y+20;j1++)
            {
            for(i1=0;i1<X+20;i1++)
                {
                s1[j1][i1]=s2[j1][i1];s2[j1][i1]=s3[j1][i1];
                sx13[j1][i1]=sx12[j1][i1],sx12[j1][i1]=sx11[j1][i1];
                sx33[j1][i1]=sx32[j1][i1],sx32[j1][i1]=sx31[j1][i1];
                sx24[j1][i1]=sx23[j1][i1],sx23[j1][i1]=sx22[j1][i1];
                sx22[j1][i1]=sx21[j1][i1];
                }
            }
    }
/////////////function wave modeling 2D END///////////////

};
//////////////////////Class END///////////////////////////

template <typename TT>
void binary_data_read_2D(TT **data_mat, int nz, int nx, const char * filename)
{
   char str[99];
   strcpy(str, filename);

   int i,j;
   TT read_data;
   ifstream infile; 

   infile.open(str,ios::binary);
   if(!infile) cout<<"file open error: "<<str<<endl;
      for(j=0;j<nx;j++)
         {for(i=0;i<nz;i++)
            {
            infile.read((char *)&read_data, sizeof(read_data));  
            data_mat[i][j]=read_data;
            }
         }
      infile.close();
}

void binary_data_read_2D(fmat & data_mat, int nz, int nx, const char * filename)
{
   char str[99];
   strcpy(str, filename);

   int i,j;
   float read_data;
   ifstream infile; 

   infile.open(str,ios::binary);
   if(!infile) cout<<"file open error: "<<str<<endl;
      for(j=0;j<nx;j++)
         {for(i=0;i<nz;i++)
            {
            infile.read((char *)&read_data, sizeof(read_data));  
            data_mat(i,j)=read_data;
            }
         }
      infile.close();
}

template <typename TT>
void binary_data_write_2D(TT **data_mat, int nz, int nx, char const *filename )
{
    char str[99];
    strcpy(str, filename);

    int i,j;
    TT outdata;

   ofstream outf;
   outf.open(str);
   for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outdata=data_mat[i][j];
         outf.write((char*)&outdata,sizeof(outdata));
         }
      }
   outf.close(); 
}

void binary_data_write_2D(fmat & data_mat, int nz, int nx, char const *filename )
{
    char str[99];
    strcpy(str, filename);

    int i,j;
    float outdata(0.0);

   ofstream outf;
   outf.open(str);
   for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outdata=data_mat(i,j);
         outf.write((char*)&outdata,sizeof(outdata));
         }
      }
   outf.close(); 
}







