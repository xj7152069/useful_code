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
fmat fmatdiv(fmat mat1, fmat mat2)
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
            a(i,j)=mat1(i,j);
            if(abs(mat2(i,j)>=0.000000001))
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

inline fmat getsubfmat(fmat & data, int n1x1,int n1x2,int n2x1,int n2x2 ){
    int submatn1(n1x2-n1x1+1),submatn2(n2x2-n2x1+1);
    fmat submat(submatn1,submatn2);
    int i,j;
    for(i=n1x1;i<=n1x2;i++){
    for(j=n2x1;j<=n2x2;j++){
        submat(i-n1x1,j-n2x1)=data(i,j);
    }}
    return submat;
}

inline void copydata2subfmat(fmat & data, fmat & subdata,\
 int n1x1,int n1x2,int n2x1,int n2x2 ){
    int submatn1(n1x2-n1x1+1),submatn2(n2x2-n2x1+1);
    int i,j;
    for(i=n1x1;i<=n1x2;i++){
    for(j=n2x1;j<=n2x2;j++){
        data(i,j)=subdata(i-n1x1,j-n2x1);
    }}
}

inline void adddata2subfmat(fmat & data, fmat & subdata,\
 int n1x1,int n1x2,int n2x1,int n2x2 ){
    int submatn1(n1x2-n1x1+1),submatn2(n2x2-n2x1+1);
    int i,j;
    for(i=n1x1;i<=n1x2;i++){
    for(j=n2x1;j<=n2x2;j++){
        data(i,j)+=subdata(i-n1x1,j-n2x1);
    }}
}

inline void muldata2subfmat(fmat & data, fmat & subdata,\
 int n1x1,int n1x2,int n2x1,int n2x2 ){
    int submatn1(n1x2-n1x1+1),submatn2(n2x2-n2x1+1);
    int i,j;
    for(i=n1x1;i<=n1x2;i++){
    for(j=n2x1;j<=n2x2;j++){
        data(i,j)*=subdata(i-n1x1,j-n2x1);
    }}
}
inline void mul2subfmat(fmat & data, float subdata,\
 int n1x1,int n1x2,int n2x1,int n2x2 ){
    int submatn1(n1x2-n1x1+1),submatn2(n2x2-n2x1+1);
    int i,j;
    for(i=n1x1;i<=n1x2;i++){
    for(j=n2x1;j<=n2x2;j++){
        data(i,j)*=subdata;
    }}
}
/////////////////////////////////////////////////////////
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

void dataread(fmat & data_mat,int nz, int nx, const char * filename)
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
    {
        for(i=0;i<nz;i++)
        {
        inf.read((char *)&read_data, sizeof(read_data));  
        data_mat(i,j)=read_data;
        }
    }
}
void dataread(fmat & data_mat,const char * filename)
{
   char str[99];
   strcpy(str, filename);

   int i,j;
   int nz(data_mat.n_rows);
   int nx(data_mat.n_cols); 
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

void dataread(fmat & data_mat, ifstream &inf)
{
    int nz(data_mat.n_rows);
    int nx(data_mat.n_cols); 
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
    fmat datacol_transpose(n3,n1,fill::zeros);
    fmat datacol(n1,n3);
    float read_data;

    infile.open(str,ios::binary);
    if(!infile) cout<<"file open error: "<<str<<endl;
    for(i=0;i<n2;i++){
    for(j=0;j<n1;j++){
        for(i=0;i<n3;i++){
        infile.read((char *)&read_data, sizeof(read_data));  
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
    fmat datacol_transpose(n3,n1,fill::zeros);
    fmat datacol(n1,n3);
    float read_data;

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
}
void dataread3d_bycol_transpose(fcube & data3d, int nz,int nx,const char * filename)
{
    char str[99];
    strcpy(str, filename);
    ifstream infile;
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices); 
    fmat datacol_transpose(n3,n1,fill::zeros);
    fmat datacol(n1,n3);
    float read_data;

    infile.open(str,ios::binary);
    if(!infile) cout<<"file open error: "<<str<<endl;
    for(k=0;k<n2;k++){
    for(j=0;j<nx;j++){
        for(i=0;i<nz;i++){
        infile.read((char *)&read_data, sizeof(read_data));  
        datacol_transpose(i,j)=read_data;
        }
    }
    datacol=datacol_transpose.t();
    data3d.col(k)=datacol;
    }
    infile.close();
}
void dataread3d_bycol_transpose(fcube & data3d, int nz,int nx,ifstream &inf)
{
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices); 
    fmat datacol_transpose(n3,n1,fill::zeros);
    fmat datacol(n1,n3);
    float read_data;

    for(k=0;k<n2;k++){
    for(j=0;j<nx;j++){
        for(i=0;i<nz;i++){
        inf.read((char *)&read_data, sizeof(read_data));  
        datacol_transpose(i,j)=read_data;
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

void datawrite3d_bycol_transpose(fcube & data3d, int nz, int nx, char const *filename)
{
    char str[99];
    strcpy(str, filename);
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices); 
    fmat datacol_transpose(n3,n1);
    fmat datacol(n1,n3);
    ofstream outf;
    outf.open(str);
    for(k=0;k<n2;k++){
        datacol=data3d.col(k);
        datacol_transpose= datacol.t();
    for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outf.write((char*)&datacol_transpose(i,j),sizeof(float));
         }
      }
    }
    outf.close(); 
}
void datawrite3d_byrow_transpose(fcube & data3d, int nz, int ny, char const *filename)
{
    char str[99];
    strcpy(str, filename);
    int i,j,k;
    int n1(data3d.n_rows),n2(data3d.n_cols),n3(data3d.n_slices); 
    fmat datarow_transpose(n3,n2);
    fmat datarow(n2,n3);
    ofstream outf;
    outf.open(str);
    for(k=0;k<n1;k++){
        datarow=data3d.row(k);
        datarow_transpose=datarow.t();
    for(j=0;j<ny;j++)
      {for(i=0;i<nz;i++)
         {
         outf.write((char*)&datarow_transpose(i,j),sizeof(float));
         }
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
void datawrite(fmat & data_mat,ofstream &outf)
{
    int nz(data_mat.n_rows);
    int nx(data_mat.n_cols); 
    int i,j;
    for(j=0;j<nx;j++)
      {for(i=0;i<nz;i++)
         {
         outf.write((char*)&data_mat(i,j),sizeof(float));
         }
      }
}
void datawrite(fmat & data_mat, char const *filename)
{
    int nz(data_mat.n_rows);
    int nx(data_mat.n_cols); 
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

template <typename T>
inline T matMulCG(T & mat1, T & mat2)
{
//Subfunctions use for solveCG()
    int nz,nx;
    nz=mat1.n_rows;
    nx=mat1.n_cols;

    T a(1,1,fill::zeros);
    int i,j;
    for(i=0;i<nz;i++){
        a+=mat1.row(i)*mat2.row(i).t();
    }
    return a;
}

template <typename T>
bool solveCG(T a, T& x, T b, T W,\
    double residual_ratio, int num)
{
/*
The conjugate gradient method solves the normal equation:
input: 
cx_fmat: a*x=b; The 'x' must be initialized and cannot be zero;
    The 'x' will be updated in subfunctions;
cx_fmat: W; Diagonal matrix W, used for regularization;
float residual_ratio: Iterative residual of the CG-method;
int num: The number of iterations of the CG-method;
    Iterations times is related to the matrix dimension,
    Generally not more than the Dimension of x;
*/
    int ip,jp,in,jn,k,iter(0);
    int np1(x.n_rows),np2(x.n_cols);
    int iterations_num=num;

//This method solve matrix equation: (aH*a+W)*x=aH*b
    b=(a.t()*b);
    a=(a.t()*a+W);

    T gradient_rk,gradient_rk_1,\
        gradient_cg_pk,gradient_cg_pk_1,\
        datatp_k,datatp_k_1,\
        recoverdatatx_uk,A_gradient_cg_pk;
    T sum_num(1,1,fill::zeros),\
        residual_pow(1,1,fill::zeros),\
        residual_k(1,1,fill::zeros);
    T beta_k(1,1),alpha_k(1,1);
    datatp_k.zeros(np1,np2);
    datatp_k_1.copy_size(datatp_k);
    gradient_rk.copy_size(datatp_k);
    gradient_rk_1.copy_size(datatp_k);
    gradient_cg_pk.copy_size(datatp_k);
    gradient_cg_pk_1.copy_size(datatp_k); 
    A_gradient_cg_pk.copy_size(datatp_k);

    iter=0;
    datatp_k=x;
    //Initial gradient
    gradient_rk=b-a*datatp_k;
    gradient_cg_pk=gradient_rk;
    A_gradient_cg_pk=a*gradient_cg_pk;
    residual_pow(0,0)=(residual_ratio*gradient_rk.n_elem);
        
    alpha_k=matMulCG<T>(gradient_rk,gradient_rk);
    alpha_k=alpha_k/(matMulCG<T>(gradient_cg_pk,A_gradient_cg_pk));
    //get new solution
    datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
    //gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;
    gradient_rk_1=b-a*datatp_k_1;

    beta_k=matMulCG<T>(gradient_rk_1,gradient_rk_1);
    beta_k=beta_k/(matMulCG<T>(gradient_rk,gradient_rk));
    gradient_cg_pk_1=gradient_rk_1+beta_k(0,0)*gradient_cg_pk;
    //updata
    datatp_k=datatp_k_1;
    gradient_rk=gradient_rk_1;
    gradient_cg_pk=gradient_cg_pk_1;
    //cal residual_pow
    sum_num.set_real(sum(sum(abs(gradient_rk)),1));
    residual_k=sum_num;

    while(iter<iterations_num && abs(residual_k(0,0))\
        >abs(residual_pow(0,0))){
        iter++;
        //std::cout<<iter<<"|"<<residual_k(0,0)<<std::endl;
        A_gradient_cg_pk=a*gradient_cg_pk;
        
        alpha_k=matMulCG<T>(gradient_cg_pk,A_gradient_cg_pk);
        if(abs(alpha_k(0,0))<=residual_ratio)break;
        alpha_k=matMulCG<T>(gradient_rk,gradient_rk)/alpha_k;
        //get new solution
        datatp_k_1=datatp_k+alpha_k(0,0)*gradient_cg_pk;
        //gradient_rk_1=gradient_rk-alpha_k(0,0)*A_gradient_cg_pk;
        gradient_rk_1=b-a*datatp_k_1;

        beta_k=matMulCG<T>(gradient_rk,gradient_rk);
        if(abs(beta_k(0,0))<=residual_ratio)break;
        beta_k=matMulCG<T>(gradient_rk_1,gradient_rk_1)/beta_k;
        gradient_cg_pk_1=gradient_rk_1+(beta_k(0,0))*gradient_cg_pk;
        //updata
        datatp_k=datatp_k_1;
        gradient_rk=gradient_rk_1;
        gradient_cg_pk=gradient_cg_pk_1;
        //cal residual_pow
        sum_num.set_real(sum(sum(abs(gradient_rk)),1));
        residual_k=sum_num;
    }
    //output result: x
    x=datatp_k;
    std::cout<<iter<<"|"<<abs(residual_k(0,0))<<std::endl;

//Determine if the CG-method converges;
//If it does not converge, it needs to 
// change the initial value to solve again;
    bool convergence;
    if(isnan(abs(residual_k(0,0)))){
            convergence=false;
        }else{
            convergence=true;
        }
    return convergence;
}

template <typename T, typename T1, typename T2,\
     typename T3, typename T4>
T solveCG_real(T1 &a, T2 &x, T3 &b, T4 &W,\
    double residual_ratio, int num)
{
    T xResult,a2,b2,w2;
    xResult.copy_size(x);
    for(int i=0;i<x.n_rows;i++){
    for(int j=0;j<x.n_cols;j++){
        xResult(i,j)=x(i,j);
    }}

    a2.copy_size(a);
    for(int i=0;i<a.n_rows;i++){
    for(int j=0;j<a.n_cols;j++){
        a2(i,j)=a(i,j);
    }}

    b2.copy_size(b);
    for(int i=0;i<b.n_rows;i++){
    for(int j=0;j<b.n_cols;j++){
        b2(i,j)=b(i,j);
    }}

    w2.copy_size(W);
    for(int i=0;i<W.n_rows;i++){
    for(int j=0;j<W.n_cols;j++){
        w2(i,j)=W(i,j);
    }}
    solveCG(a2, xResult, b2, w2, residual_ratio, num);
    return xResult;
}

template <typename T, typename T1, typename T2,\
     typename T3, typename T4>
T solveCG_complex(T1 &a, T2 &x, T3 &b, T4 &W,\
    double residual_ratio, int num)
{
    T xResult,a2,b2,w2;
    xResult.copy_size(x);
    for(int i=0;i<x.n_rows;i++){
    for(int j=0;j<x.n_cols;j++){
        xResult(i,j).real(real(x(i,j)));
        xResult(i,j).imag(imag(x(i,j)));
    }}

    a2.copy_size(a);
    for(int i=0;i<a.n_rows;i++){
    for(int j=0;j<a.n_cols;j++){
        a2(i,j).real(real(a(i,j)));
        a2(i,j).imag(imag(a(i,j)));
    }}

    b2.copy_size(b);
    for(int i=0;i<b.n_rows;i++){
    for(int j=0;j<b.n_cols;j++){
        b2(i,j).real(real(b(i,j)));
        b2(i,j).imag(imag(b(i,j)));
    }}

    w2.copy_size(W);
    for(int i=0;i<W.n_rows;i++){
    for(int j=0;j<W.n_cols;j++){
        w2(i,j).real(real(W(i,j)));
        w2(i,j).imag(imag(W(i,j)));
    }}
    solveCG(a2, xResult, b2, w2, residual_ratio, num);
    return xResult;
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

fmat Laplace(fmat mat0)
    {
        int i,j,n;
        fmat mat2,mat1;
        int x1(mat0.n_rows), x2(mat0.n_cols);
        mat2=mat0;
        mat1=mat0;
        for(i=1;i<x1-1;i++)
        {
            for(j=1;j<x2-1;j++)
            {
                mat1(i,j)=mat2(i-1,j-1)*1.0/12+mat2(i-1,j)*1.0/6+mat2(i-1,j+1)*1.0/12\
                    +mat2(i,j-1)*1.0/6+mat2(i,j+1)*1.0/6+mat2(i+1,j-1)*1.0/12\
                    +mat2(i+1,j)*1.0/6+mat2(i+1,j+1)*1.0/12-mat2(i,j);
            }
        }
        return mat1;
    }


#endif

