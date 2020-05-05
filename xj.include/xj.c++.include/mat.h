/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef MAT_H_H
#define MAT_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////
//建议还是在指针数组所在的语句块中释放内存;
//若要调用以下函数释放空间，应在原语句块中再次将指针空置。
template<typename T1> void matdelete(T1 **mat, int x1);
template<typename T1> void matdelete(T1 ***mat, int x1, int x2);
//数组定义原则：越有可能作为整体被调用的越往后；
//因为对于如三维指针***p而言，*p和**p都是地址数组。
float** newfmat(int x1, int x2);
float*** newfmat(int x1, int x2, int x3);
template <typename T1> void matprint(T1 *mat1, int n);
template <typename T1> void matprint(T1 **mat1, int nz, int nx);
template <typename T1> void matprint(T1 ***mat1, int n1, int n2, int n3);
template <typename T1> char* numtostr(T1 n, int s);
template <typename T1, typename T2> void matdec(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename T1, typename T2> void matadd(T1 *mat1, T2 *mat2, int n);
template <typename T1, typename T2> void matadd(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> void matadd(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename T1, typename T2> void matcopy(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> void matcopy(T1 ***mat1, T2 n, int n1, int n2, int n3);
template <typename T1, typename T2> void matcopy(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename T1, typename T2> void matcopy(T1 *mat1, T2 *mat2, int n);
template <typename T1, typename T2> void matmul(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> void matmul(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename TT> void dataread(TT **data_mat, int nz, int nx, const char * filename);
template <typename TT> void dataread(TT **data_mat, int nz, int nx, ifstream &inf);
template <typename TT> void datawrite(TT **data_mat, int nz, int nx, const char *filename);
template <typename TT> void datawrite(TT **data_mat, int nz, int nx, ofstream &outf);
template <typename T1, typename T2> void matsmooth(T1 **mat1, T2 **mat0, int x1, int x2, int k=1);
///////////////////////////////////////////////////////////////////////////////////////////

template<typename T1>
char* numtostr(T1 n, int s)
{
    char *str=NULL;
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
void matcopy(T1 *mat1, T2 *mat2, int n)
{
    int i;
    for(i=0;i<n;i++)
        {
        mat1[i]=mat2[i];
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
void matadd(T1 *mat1, T2 *mat2, int n)
{
    int i;
    for(i=0;i<n;i++)
        {
        mat1[i]+=mat2[i];
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

template <typename T1, typename T2>
void matsmooth(T1 **mat1, T2 **mat0, int x1, int x2, int k)
{
    int i,j,n;
    float **mat2;
    mat2=newfmat(x1, x2);
    matcopy(mat2,mat0,x1,x2);
    
    for(n=0;n<k;n++)
    {
        for(i=1;i<x1-1;i++)
        {
            for(j=1;j<x2-1;j++)
            {
                mat1[i][j]=(mat2[i-1][j-1]*1.0/12+mat2[i-1][j]*1.0/6+mat2[i-1][j+1]*1.0/12\
                        +mat2[i][j-1]*1.0/6+mat2[i][j+1]*1.0/6+mat2[i+1][j-1]*1.0/12\
                        +mat2[i+1][j]*1.0/6+mat2[i+1][j+1]*1.0/12+mat2[i][j]*1.0/3)*(3.0/4.0);
            }
        }
        matcopy(mat2,mat1,x1,x2);
    }
    matdelete(mat2, x1);
}

#endif

