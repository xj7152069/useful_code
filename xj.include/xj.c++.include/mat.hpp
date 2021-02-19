/*********(version 1.0)***********/
/*
wave2D.h
    c++ head file: 
*/
/********************************/
#ifndef MAT_HPP
#define MAT_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "../xjc.h"

/////////////////////////////////////////////////////////////////////////////////////////////
//���黹����ָ���������ڵ��������ͷ��ڴ�;
//��Ҫ�������º����ͷſռ䣬Ӧ��ԭ�������ٴν�ָ����á�
template<typename T1> void matdelete(T1 **mat, int x1);
template<typename T1> void matdelete(T1 ***mat, int x1, int x2);
//���鶨��ԭ��Խ�п�����Ϊ���屻���õ�Խ����
//��Ϊ��������άָ��***p���ԣ�*p��**p���ǵ�ַ���顣
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
template <typename TT> void datawrite(TT *data_mat, int n, const char *filename);
template <typename TT> void datawrite(TT *data_mat, int n, ofstream &outf);
template <typename T1, typename T2> void matsmooth(T1 **mat1, T2 **mat0, int x1, int x2, int k=1);
template <typename T1, typename T2> void Laplace(T1 **mat1, T2 **mat0, int x1, int x2);
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
void datawrite(TT **data_mat, int n, char const *filename )
{
    char str[99];
    strcpy(str, filename);
    int i;
    TT outdata;
    ofstream outf;
    outf.open(str);
    for(i=0;i<n;i++)
    {
        outdata=i;
        outf.write((char*)&outdata,sizeof(outdata));
        outdata=data_mat[i];
        outf.write((char*)&outdata,sizeof(outdata));
    }
    outf.close(); 
}
template <typename TT>
void datawrite(TT **data_mat, int n, ofstream &outf)
{
    int i;
    TT outdata;
    for(i=0;i<n;i++)
    {
        outdata=i;
        outf.write((char*)&outdata,sizeof(outdata));
        outdata=data_mat[i];
        outf.write((char*)&outdata,sizeof(outdata));
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

//���鶨��ԭ��Խ�п�����Ϊ���屻���õ�Խ����
//��Ϊ��������άָ��***p���ԣ�*p��**p���ǵ�ַ���顣
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

//�����ͷ��ڴ棬���޷���ԭ�����е�ָ����ã�
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
    matcopy(mat1,mat0,x1,x2);
    
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

template <typename T1, typename T2>
void Laplace(T1 **mat1, T2 **mat0, int x1, int x2)
    {
        int i,j,n;
        float **mat2;
        mat2=newfmat(x1, x2);
        matcopy(mat2,mat0,x1,x2);
        matcopy(mat1,mat0,x1,x2);
        for(i=1;i<x1-1;i++)
        {
            for(j=1;j<x2-1;j++)
            {
                mat1[i][j]=mat2[i-1][j-1]*1.0/12+mat2[i-1][j]*1.0/6+mat2[i-1][j+1]*1.0/12\
                    +mat2[i][j-1]*1.0/6+mat2[i][j+1]*1.0/6+mat2[i+1][j-1]*1.0/12\
                    +mat2[i+1][j]*1.0/6+mat2[i+1][j+1]*1.0/12-mat2[i][j];
            }
        }
    }

#endif

