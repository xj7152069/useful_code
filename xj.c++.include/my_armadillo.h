/*********(version 1.0)***********/
/*
my_armadillo.h
    c++ head file: some useful function
*/
/********************************/
#ifndef MY_ARMADILLO_H_H
#define MY_ARMADILLO_H_H

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
using namespace std;
using namespace arma;

///////////////////////////////////////////////////////////////////////////////////
template <typename T2> void matcopy(T2 **mat2, fmat &data_mat, int nz, int nx);
template <typename TT> fmat matcopy(TT **mat, int nz, int nx);
fmat matmul(fmat mat1, fmat mat2, int nz, int nx);
template <typename T1> fmat matmul(fmat mat1, T1 n, int nz, int nx);
void dataread(fmat & data_mat, int nz, int nx, const char * filename);
void dataread(fmat & data_mat, int nz, int nx, ifstream &inf);
void datawrite(fmat & data_mat, int nz, int nx, char const *filename);
void datawrite(fmat & data_mat, int nz, int nx, ofstream &outf);
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

void dataread(fmat & data_mat, int nz, int nx, const char * filename)
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

void dataread(fmat & data_mat, int nz, int nx, ifstream &inf)
{
    int i,j;
    float read_data;

    for(j=0;j<nx;j++)
        {for(i=0;i<nz;i++)
            {
            inf.read((char *)&read_data, sizeof(read_data));  
            data_mat(i,j)=read_data;
            }
        }
}

void datawrite(fmat & data_mat, int nz, int nx, char const *filename)
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

void datawrite(fmat & data_mat, int nz, int nx, ofstream &outf)
{
    int i,j;
    float outdata(0.0);

    for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outdata=data_mat(i,j);
         outf.write((char*)&outdata,sizeof(outdata));
         }
      }
}






#endif

