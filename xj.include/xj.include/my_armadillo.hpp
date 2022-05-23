/*********(version 1.0)***********/
/*
my_armadillo.h
    c++ head file: some useful function
*/
/********************************/
#ifndef MY_ARMADILLO_HPP
#define MY_ARMADILLO_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "../xjc.h"
using namespace std;
using namespace arma;

///////////////////////////////////////////////////////////////////////////////////
fmat matdiv(fmat mat1, fmat mat2, int nz, int nx, float min=0.000000000001);
template <typename T2> void matcopy(T2 **mat2, fmat &data_mat, int nz, int nx);
template <typename TT> fmat matcopy(TT **mat, int nz, int nx);
fmat matmul(fmat mat1, fmat mat2, int nz, int nx);
template <typename T1> fmat matmul(fmat mat1, T1 n, int nz, int nx);
fmat dataread(int nz, int nx, const char * filename);
fmat dataread(int nz, int nx, ifstream &inf);
void datawrite(fmat & data_mat, int nz, int nx, char const *filename);
void datawrite(fmat & data_mat, int nz, int nx, ofstream &outf);
fmat fmatsmooth(fmat mat2, int x1, int x2, int k);
///////////////////////////////////////////////////////////////////////////////////

template <typename T2>
void matcopy(T2 **mat2, fmat &data_mat, int nz, int nx)
{
    
    int i,j;
    for(i=0;i<nz;i++)
        {
        for(j=0;j<nx;j++)
            {
            mat2[i][j]=data_mat(i,j);
            }
        }
}

template <typename TT>
fmat matcopy(TT **mat, int nz, int nx)
{
    fmat a(nz,nx);
    int i,j;

    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
            {
            a(i,j)=float(mat[i][j]);
            }
    }

    return a;
}

fmat matmul(fmat mat1, fmat mat2, int nz, int nx)
{
    fmat a(nz,nx);
    int i,j;

    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
            {
            a(i,j)=mat1(i,j)*mat2(i,j);
            }
    }

    return a;
}

fmat fmatmul(fmat mat1, fmat mat2)
{
    int nz,nx;
    nz=mat1.n_rows;
    nx=mat1.n_cols;
    fmat a(nz,nx);
    int i,j;

    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
            {
            a(i,j)=mat1(i,j)*mat2(i,j);
            }
    }
    return a;
}

fmat matdiv(fmat mat1, fmat mat2, int nz, int nx, float min)
{
    fmat a(nz,nx);
    int i,j;

    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
            {
            a(i,j)=mat1(i,j);
            if(abs(mat2(i,j)>=abs(min)))
                {a(i,j)=mat1(i,j)/mat2(i,j);}
            
            }
    }

    return a;
}

template <typename T1>
fmat matmul(fmat mat1, T1 n, int nz, int nx)
{
    fmat a(nz,nx);
    int i,j;

    for(i=0;i<nz;i++)
    {
        for(j=0;j<nx;j++)
            {
            a(i,j)=mat1(i,j)*n;
            }
    }

    return a;
}

fmat dataread(int nz, int nx, const char * filename)
{
   char str[99];
   strcpy(str, filename);
   fmat data_mat(nz,nx);

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
    return data_mat;
}

void dataread(fmat & data_mat, int nz, int nx, ifstream &inf)
{
    int i,j;
    float read_data;

    for(j=0;j<nx;j++)
    {
        for(i=0;i<nz;i++)
        {
        inf.read((char *)&read_data, sizeof(read_data));  
        data_mat(i,j)=read_data;
        }
    }
}

void dataread3d_bycol_transpose(fcube & data3d, const char * filename)
{
    char str[99];
    strcpy(str, filename);
    ifstream infile;
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices); 
    fmat datacol_transpose(n3,n1);
    fmat datacol(n1,n3);
    float read_data;

    infile.open(str,ios::binary);
    if(!infile) cout<<"file open error: "<<str<<endl;
    for(i=0;i<n2;i++){
    for(j=0;j<n1;j++){
        for(i=0;i<n3;i++){
        inf.read((char *)&read_data, sizeof(read_data));  
        datacol_transpose(i,j)=read_data;
        }
    }
    datacol=datacol_transpose.t();
    data3d.col(k)=datacol;
    }
    infile.close();
}
void dataread3d_bycol_transpose(fcube & data3d, ifstream &inf)
{
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices); 
    fmat datacol_transpose(n3,n1);
    fmat datacol(n1,n3);
    float read_data;

    for(i=0;i<n2;i++){
    for(j=0;j<n1;j++){
        for(i=0;i<n3;i++){
        inf.read((char *)&read_data, sizeof(read_data));  
        datacol_transpose(j,i)=read_data;
        }
    }
    datacol=datacol_transpose.t();
    data3d.col(k)=datacol;
    }
}

fmat dataread(int nz, int nx, ifstream &inf)
{
    int i,j;
    float read_data;
    fmat data_mat(nz,nx);

    for(j=0;j<nx;j++)
        {for(i=0;i<nz;i++)
            {
            inf.read((char *)&read_data, sizeof(read_data));  
            data_mat(i,j)=read_data;
            }
        }
    return data_mat;
}

void datawrite(fmat & data_mat, int nz, int nx, char const *filename)
{
    char str[99];
    strcpy(str, filename);
    int i,j;
    ofstream outf;
    outf.open(str);
    for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         //outdata=float(data_mat(i,j));
         outf.write((char*)&data_mat(i,j),sizeof(float));
         }
      }
    outf.close(); 
}

void datawrite(fmat & data_mat, int nz, int nx, ofstream &outf)
{
    int i,j;
    for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outf.write((char*)&data_mat(i,j),sizeof(float));
         }
      }
}

fmat fmatsmooth(fmat mat2, int x1, int x2, int k)
{
    int i,j,n;
    fmat mat1(x1,x2);
    mat1=mat2;

    for(n=0;n<k;n++)
    {
        for(i=0;i<x1;i++)
        {
            for(j=0;j<x2;j++)
            {
                if(i>=1 && i<x1-1 && j>=1 && j<x2-1)
                {
                    mat1(i,j)=((mat2(i-1,j-1)*1.0/12+mat2(i-1,j)*1.0/6+mat2(i-1,j+1)*1.0/12\
                        +mat2(i,j-1)*1.0/6+mat2(i,j+1)*1.0/6+mat2(i+1,j-1)*1.0/12\
                        +mat2(i+1,j)*1.0/6+mat2(i+1,j+1)*1.0/12+mat2(i,j)*1.0/3))*(3.0/4.0);
                }
                if(i==0)
                {mat1(i,j)=mat1(i+1,j);}
                if(j==0)
                {mat1(i,j)=mat1(i,j+1);}
                if(i==x1-1)
                {mat1(i,j)=mat1(i-1,j);}
                if(j==x2-1)
                {mat1(i,j)=mat1(i,j-1);}
            }
        }
        mat2=mat1;
    }
    return mat1;
}




#endif

