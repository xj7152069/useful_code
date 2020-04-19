/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef WAVE2D_H_H
#define WAVE2D_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;

///////////////////////////////////////////////////////
//建议还是在指针数组所在的语句块中释放内存;
//若要调用以下函数释放空间，应在原语句块中再次将指针空置。
template<typename T1>
void matdelete(T1 **mat, int x1);

template<typename T1>
void matdelete(T1 ***mat, int x1, int x2);

//数组定义原则：越有可能作为整体被调用的越往后；
//因为对于如三维指针***p而言，*p和**p都是地址数组。
float** newfmat(int x1, int x2);

float*** newfmat(int x1, int x2, int x3);

void wave2Dtest(int Z, int X);

template<typename T1>
void matprint(T1 *mat1, int n);

template<typename T1>
void matprint(T1 **mat1, int nz, int nx);

template<typename T1>
void matprint(T1 ***mat1, int n1, int n2, int n3);

template<typename T1>
char* numtostr(T1 n, int s);

float wavelet01(int k, float DT, float hz=35.0, int delay=70);

float wavelet02(int k, float DT, float hz=35.0, int delay=70);

template <typename T1, typename T2>
void matdec(T1 **mat1, T2 **mat2, int nz, int nx);

template <typename T1, typename T2>
void matadd(T1 **mat1, T2 n, int nz, int nx);

template <typename T1, typename T2>
void matadd(T1 **mat1, T2 **mat2, int nz, int nx);

template <typename T1, typename T2>
void matcopy(T1 **mat1, T2 n, int nz, int nx);

template <typename T1, typename T2>
void matcopy(T1 ***mat1, T2 n, int n1, int n2, int n3);

template <typename T1, typename T2>
void matcopy(T1 **mat1, T2 **mat2, int nz, int nx);

template <typename T1, typename T2>
void matmul(T1 **mat1, T2 n, int nz, int nx);

template <typename T1, typename T2>
void matmul(T1 **mat1, T2 **mat2, int nz, int nx);

template <typename TT>
void dataread(TT **data_mat, int nz, int nx, const char * filename);

template <typename TT>
void dataread(TT **data_mat, int nz, int nx, ifstream &inf);

template <typename TT>
void datawrite(TT **data_mat, int nz, int nx, const char *filename);

template <typename TT>
void datawrite(TT **data_mat, int nz, int nx, ofstream &outf);

class wave2D
{
public:
    float **p2=NULL; //velocity model
    float **s1=NULL, **s2=NULL, **s3=NULL; //time slices, add source to "s2"
    float **sx11=NULL, **sx12=NULL, **sx13=NULL; //PML boundary
    float **sx21=NULL, **sx22=NULL, **sx23=NULL, **sx24=NULL; //PML boundary
    float **sx31=NULL,  **sx32=NULL, **sx33=NULL; //PML boundary

    float xs2[5]={1.666667,-0.238095,0.039683,-0.004960,0.000317};
    float xs1[10]={-0.0007926,0.00991800,-0.0595200,0.238080,-0.833333,\
0.833333,-0.238080,0.0595200,-0.00991800,0.0007926};

    float dx,dy,dt,PML_wide,R;
    int nx,ny,suface;

    wave2D();
    wave2D(int x, int y);
    ~wave2D();

    template<typename T1>
    void setvelocity(T1 v=3000);
    void timeslicecal();
    void timeslicecopy();
    void cleardata();
};

///////////////////////////////////////////////////////
wave2D::wave2D()
{
    nx=0;
    ny=0;
    dx=5.0;dy=5.0;dt=0.0005;PML_wide=25;suface=1;R=10000000;
    cout<<"Warning: Creat an Empty object-wave_modeling_2D"<<endl;
}

wave2D::wave2D(int z, int x)
{
    nx=x;
    ny=z;
    dx=5.0;dy=5.0;dt=0.0005;PML_wide=20;suface=1;R=100000;
    int i,j;
    s1=new float*[ny];  
    s2=new float*[ny];   
    s3=new float*[ny];  
    sx11=new float*[ny];    
    sx12=new float*[ny];     
    sx13=new float*[ny];     
    sx21=new float*[ny];     
    sx22=new float*[ny];     
    sx23=new float*[ny];     
    sx24=new float*[ny];     
    sx31=new float*[ny];     
    sx32=new float*[ny];     
    sx33=new float*[ny];      
    for(j=0;j<ny;j++)  
        {  
        s1[j]=new float[nx];  
        s2[j]=new float[nx];   
        s3[j]=new float[nx];   
        sx11[j]=new float[nx];   
        sx12[j]=new float[nx];   
        sx13[j]=new float[nx];    
        sx21[j]=new float[nx];   
        sx22[j]=new float[nx];   
        sx23[j]=new float[nx];   
        sx24[j]=new float[nx];   
        sx31[j]=new float[nx];   
        sx32[j]=new float[nx];    
        sx33[j]=new float[nx];   
        }

    for(j=0;j<ny;j++)
        {
        for(i=0;i<nx;i++)
            {
            s1[j][i]=0.0;s2[j][i]=0.0;s3[j][i]=0.0;
            sx11[j][i]=0.0;sx12[j][i]=0.0;sx13[j][i]=0.0;
            sx21[j][i]=0.0;sx22[j][i]=0.0;sx23[j][i]=0.0;sx24[j][i]=0.0;
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
       }
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
    cout<<"Delete an object-wave_modeling_2D"<<endl;
}

template<typename T1>
void wave2D::setvelocity(T1 v)
{
    int i,j;
    for(i=0;i<ny;i++)
        {
        for(j=0;j<nx;j++)
            {
            p2[i][j]=v;
            }
        }
}

void wave2D::cleardata()
{
    this->setvelocity(0.0);
    matcopy(s1,0.0,ny,nx);
    matcopy(s2,0.0,ny,nx);
    matcopy(s3,0.0,ny,nx);
    matcopy(sx11,0.0,ny,nx);
    matcopy(sx12,0.0,ny,nx);
    matcopy(sx13,0.0,ny,nx);
    matcopy(sx21,0.0,ny,nx);
    matcopy(sx22,0.0,ny,nx);
    matcopy(sx23,0.0,ny,nx);
    matcopy(sx24,0.0,ny,nx);
    matcopy(sx31,0.0,ny,nx);
    matcopy(sx32,0.0,ny,nx);
    matcopy(sx33,0.0,ny,nx);
    cout<<"All Matrix data has been clear!"<<endl;
}

void wave2D::timeslicecal()
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
 
    float dx,dy,ddx,ddy,snx1,sny1,snx2,sny2,t2,t5;
    int i,j,n,i1,j1,t3,t4;
    float u1(0),u2(0),u(0),ux(0),uy(0);
    float C3=3/2/(xshd)/DX*log(R)/(xshd)/DX/(xshd)/DX;
    t5=float(Y)/X;

    for(i=5;i<X-5;i++)
        {
        for(j=5;j<Y-5;j++)
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
                if(i>=X-xshd-5 && j<=t5*i && j>=-t5*i+Y)	
                    snx1=i-(X-xshd-5);
                if(i<=xshd+5 && j>=t5*i && j<=-t5*i+Y)		
                    snx2=xshd+5-i;
                if(j>=Y-xshd-5 && j>t5*i && j>-t5*i+Y)		
                    sny1=j-(Y-xshd-5);
                if(j<=xshd+5 && j<t5*i && j<-t5*i+Y)		
                    sny2=xshd+5-j;		
                }  //处理角落
            else if(suface_PML==0)
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
                dy=p2[j][i]*C3*sny1*sny1*DY*DY;
                ddy=p2[j][i]*C3*2*sny1*DY;	
                }
            if(sny2 !=0)
                {
                dy=p2[j][i]*C3*sny2*sny2*DY*DY;
                ddy=p2[j][i]*C3*2*sny2*DY;	
                }
            if(snx1 !=0)
                {
                dx=p2[j][i]*C3*snx1*snx1*DX*DX;
                ddx=p2[j][i]*C3*2*snx1*DX;	
                }
            if(snx2 !=0)
                {
                dx=p2[j][i]*C3*snx2*snx2*DX*DX;
                ddx=p2[j][i]*C3*2*snx2*DX;	
                }

            if(i>=0.5*(X))
                t3=1;
            else
                t3=-1;

            if(j>=0.5*(Y))
                t4=1;
            else
                t4=-1;

            if(snx1!=0 || snx2!=0)
                {
                sx11[j][i]=p2[j][i]*p2[j][i]\
                *DT*DT*(u2-u*s2[j][i])*(1.0/(DX*DX))\
                -dx*dx*DT*DT*sx12[j][i]+(2*sx12[j][i]\
                -sx13[j][i])+DT*(2*dx*(sx13[j][i]-sx12[j][i]));

                sx21[j][i]=(-p2[j][i]*p2[j][i]\
                *ddx*(1.0/(DX))*(ux*t3)-dx*dx*dx*sx22[j][i]\
                +3*dx*dx*(sx23[j][i]-sx22[j][i])/DT+3*dx*\
                (2*sx23[j][i]-sx22[j][i]-sx24[j][i])/(DT*DT)\
                +(3*sx22[j][i]-3*sx23[j][i]+sx24[j][i])\
                /(DT*DT*DT))*(DT*DT*DT);

                sx31[j][i]=DT*DT*p2[j][i]*p2[j][i]\
                *(1.0/(DX*DX))*(u1-u*s2[j][i])+2*sx32[j][i]\
                -sx33[j][i];
                }
            else
                {
                sx11[j][i]=p2[j][i]*p2[j][i]\
                *DT*DT*(u1-u*s2[j][i])*(1.0/(DX*DX))\
                -dy*dy*DT*DT*sx12[j][i]+(2*sx12[j][i]\
                -sx13[j][i])+DT*(2*dy*(sx13[j][i]-sx12[j][i]));

                sx21[j][i]=(-p2[j][i]*p2[j][i]\
                *ddy*(1.0/(DX))*(uy*t4)-dy*dy*dy*sx22[j][i]\
                +3*dy*dy*(sx23[j][i]-sx22[j][i])/DT\
                +3*dy*(2*sx23[j][i]-sx22[j][i]-sx24[j][i])\
                /(DT*DT)+(3*sx22[j][i]-3*sx23[j][i]+sx24[j][i])\
                /(DT*DT*DT))*(DT*DT*DT);

                sx31[j][i]=DT*DT*p2[j][i]*p2[j][i]\
                *(1.0/(DX*DX))*(u2-u*s2[j][i])+2*sx32[j][i]\
                -sx33[j][i];
                }   

            // if(snx1!=0 || snx2!=0 || sny1!=0 || sny2!=0)
            s3[j][i]=sx11[j][i]+sx21[j][i]+sx31[j][i];
            u1=0,u2=0,u=0,ux=0,uy=0;
            }
        }
}

void wave2D::timeslicecopy()
{
    int j1,i1;
    for(j1=0;j1<ny;j1++)
    {
    for(i1=0;i1<nx;i1++)
        {
        s1[j1][i1]=s2[j1][i1];s2[j1][i1]=s3[j1][i1];
        sx13[j1][i1]=sx12[j1][i1],sx12[j1][i1]=sx11[j1][i1];
        sx33[j1][i1]=sx32[j1][i1],sx32[j1][i1]=sx31[j1][i1];
        sx24[j1][i1]=sx23[j1][i1],sx23[j1][i1]=sx22[j1][i1];
        sx22[j1][i1]=sx21[j1][i1];
        }
    }
}

void wave2Dtest(int Z, int X, int T)
{
    wave2D A(Z,X);
    ofstream outf1;
    outf1.open("testmovie.bin");

    int k;
    float f;
    for(k=0;k<T;k++)
    {
        f=wavelet01(k,A.dt);
        A.s2[Z/2][X/2]=A.s2[Z/2][X/2]+f;
        A.timeslicecal();
        A.timeslicecopy();
        datawrite(A.s3, Z, X, outf1);
    }
    outf1.close();
    cout<<"Have output test movie."<<endl;
}

///////////////////////////////////////////////////////////

template<typename T1>
char* numtostr(T1 n, int s)
{
    static char *str=NULL;
    if(s>99)
    {
        cout<<"too long for string!"<<endl;
    }
    else
    {
        double num;
        num=n;
        str=new char[199];
        sprintf(str, "%.99f", num);
        str[s]='\0';
    }
    return str;
}

template <typename TT>
void dataread(TT **data_mat, int nz, int nx, const char * filename)
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

template <typename TT>
void dataread(TT **data_mat, int nz, int nx, ifstream &inf)
{
    int i,j;
    TT read_data;

    for(j=0;j<nx;j++)
        {for(i=0;i<nz;i++)
            {
            inf.read((char *)&read_data, sizeof(read_data));  
            data_mat[i][j]=read_data;
            }
        }
}

template <typename TT>
void datawrite(TT **data_mat, int nz, int nx, char const *filename )
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

template <typename TT>
void datawrite(TT **data_mat, int nz, int nx, ofstream &outf)
{
    int i,j;
    TT outdata;

    for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outdata=data_mat[i][j];
         outf.write((char*)&outdata,sizeof(outdata));
         }
      }
}

template <typename T1, typename T2>
void matcopy(T1 **mat1, T2 **mat2, int nz, int nx)
{
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=mat2[i][j];
            }
        }
}

template <typename T1, typename T2>
void matcopy(T1 **mat1, T2 n, int nz, int nx)
{
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=n;
            }
        }
}

template <typename T1, typename T2>
void matcopy(T1 ***mat1, T2 n, int n1, int n2, int n3)
{
    int i,j,k;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            for(k=0;k<n3;k++)
            {
                mat1[i][j][k]=n;
            }
        }
    }
}

template <typename T1, typename T2>
void matadd(T1 **mat1, T2 **mat2, int nz, int nx)
{
    
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=mat1[i][j]+mat2[i][j];
            }
        }
}

template <typename T1, typename T2>
void matadd(T1 **mat1, T2 n, int nz, int nx)
{
    
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=mat1[i][j]+n;
            }
        }
}

template <typename T1, typename T2>
void matdec(T1 **mat1, T2 **mat2, int nz, int nx)
{ 
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=mat1[i][j]-mat2[i][j];
            }
        }
}

template <typename T1, typename T2>
void matmul(T1 **mat1, T2 n, int nz, int nx)
{
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=mat1[i][j]*n;
            }
        }
}

template <typename T1, typename T2>
void matmul(T1 **mat1, T2 **mat2, int nz, int nx)
{
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat1[i][j]=mat1[i][j]*mat2[i][j];
            }
        }
}

template<typename T1>
void matprint(T1 *mat1, int n)
{
    int i;
    for(i=0;i<n;i++)
    {
        cout<<" ("<<mat1[i]<<") ";
    }
    cout<<endl;
}

template<typename T1>
void matprint(T1 **mat1, int nz, int nx)
{
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            cout<<" ("<<mat1[i][j]<<") ";
            }
        cout<<endl;
        }
}

template<typename T1>
void matprint(T1 ***mat1, int n1, int n2, int n3)
{
    int i,j,k;
    for(k=0;k<n1;k++)
    {
        for(i=0;i<n2;i++)
            {
            for(j=0;j<n3;j++)
                {
                cout<<" ("<<mat1[k][i][j]<<") ";
                }
            cout<<endl;
            }
        cout<<endl;
    }
}

float wavelet02(int k, float DT, float hz, int delay)
{
    float pi(3.1415926);
	float f,det;
    det=delay*DT;

    f=(pi)*(pi)*hz*hz*(k*DT-det)*\
    exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
    *(3.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));

	return f;
}

float wavelet01(int k, float DT, float hz, int delay)
{
    float pi(3.1415926);
	float f,det;
    det=delay*DT;

    f=exp((-pi*pi*hz*hz*(k*DT-det)*(DT*k-det)))\
    *(1.0-2.0*pi*pi*hz*hz*(k*DT-det)*(DT*k-det));

	return f;
}

//数组定义原则：越有可能作为整体被调用的越往后；
//因为对于如三维指针***p而言，*p和**p都是地址数组。
float** newfmat(int x1, int x2)
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
    return p2;
}

float*** newfmat(int x1, int x2, int x3)
{
    float ***p,***p2;
    int i,j;
    p=new float**[x1];      
    for(j=0;j<x1;j++)  
        {  
        p[j]=new float*[x2];
        for(i=0;i<x2;i++)
            {
                p[j][i]=new float[x3];
            }
        }
    p2=p;
    p=NULL;
    return p2;
}

//可以释放内存，但无法将原语句块中的指针空置？
template<typename T1>
void matdelete(T1 **mat, int x1)
{
    int i;
    for(i=0;i<x1;i++)
    {
        delete []mat[i];
        mat[i]=NULL;
    }
    delete []mat;
    mat=NULL;
}

template<typename T1>
void matdelete(T1 ***mat, int x1, int x2)
{
    int i,j;
    for(i=0;i<x1;i++)
    {
        for(j=0;j<x2;j++)
        {
            delete []mat[i][j];
            mat[i][j]=NULL;
        }
        delete []mat[i];
        mat[i]=NULL;
    }    
    delete []mat;
    mat=NULL;
}




#endif

