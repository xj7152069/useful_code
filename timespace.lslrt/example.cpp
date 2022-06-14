#include <fstream>
#include <string.h>
#include <iostream>
#include "./include/slantstack3d.hpp"
using namespace std;

float** newfmat(int x1, int x2);
template<typename T1>
void matdelete(T1 **mat, int x1);
template <typename TT> void dataread(TT **data_mat, int nz, int nx, const char * filename);
template <typename TT> void dataread(TT **data_mat, int nz, int nx, ifstream &inf);
template <typename TT> void datawrite(TT **data_mat, int nz, int nx, const char *filename);
template <typename TT> void datawrite(TT **data_mat, int nz, int nx, ofstream &outf);

// 2-D Slant-stack in the time - space domain.
void Slantstack_CG_2D(float **recoverdata, float **tauppanel, float **trace,\
 float *coor, int ntrace, float x0, int ns, float dt,\
 int npsr,float psrmin, float dpsr,\
 float factor_L2,float factor_L1, \
 int iterations_num, float residual_ratio);
int main()
{
    int i,j,k,usecopy(0),i1,j1,k1;
    //input 3D data
    string infile;
    int nz2(1024),nz(1024),nx(200),ny(1),nf(0);
    float dt(0.001);
    float **recover,**trace,**tauppanel,*coor;
    ifstream inf1; 
    ofstream outf1; 
    float dx(5),x0(dx*nx/2);
    int npsr(501);
    float dpsr(0.000005);
    float psrmin(-dpsr*npsr/2);

    float **basedata(newfmat(nz2,nx)),\
        **recoverdata(newfmat(nz2,nx)),\
        **datatp(newfmat(nz2,npsr));
    trace=newfmat(nx, nz2);
    recover=newfmat(nx, nz2);
    tauppanel=newfmat(npsr, nz2);
    coor=new float[nx];

    inf1.open("txz.two.1024.200.bin");  //local data
    dataread(basedata,nz2,nx,inf1);
    inf1.close();


    for(i=0;i<nx;i++){
        coor[i]=dx*i-x0;
        for(j=0;j<nz2;j++){
            trace[i][j]=basedata[j][i];
        }
    }
    //medianfilter(trace, nx, nz2);

//void Slantstack_CG_2D(float **recoverdata, float **tauppanel, float **trace,\
float *coor, int ntrace, float x0, int ns, float dt,\
int npsr,float psrmin, float dpsr,\
float factor_L2,float factor_L1, \
int iterations_num, float residual_ratio)
    //Slantstack_CG_2D(recover,tauppanel,trace,\
            coor,nx, x0,nz2,nz2, dt,\
            npsr,psrmin, dpsr,\
            1.0,1000.0, 40,0.1,true);
//void apply_LocalLSLRT_to_CMPgather(int ntrCMP, float *offset_site, int ns, float dt_s,\
 float **gather, float offsetWidth, int nphr, float phrmin, float dphr, int ntau,\
 float **taup2d_lslrt, int mintrace, float factor_L2,float factor_L1, \
 int iterations_num, float residual_ratio);
    //apply_LocalLSLRT_to_CMPgather(nx, coor, nz2, dt,\
        trace, 400, npsr, psrmin, dpsr, nz2,\
        tauppanel, 20, 1.0,100.0, 45,0.1);
//void apply_GlobalLSLRT_to_CMPgather(int ntrCMP, float *offset, int ns, float dt_s,\
 float **gather, int nphr, float phrmin, float dphr, int ntau,float **taup2d_lslrt,\
 float factor_L2,float factor_L1, \
 int iterations_num, float residual_ratio);
 apply_GlobalLSLRT_to_CMPgather(nx, coor, nz2, dt,\
        trace, npsr, psrmin, dpsr, nz2,\
        tauppanel, 1000.0,0.0, 45,0.1);
            

    for(i=0;i<nx;i++){
        for(j=0;j<nz2;j++){
            recoverdata[j][i]=recover[i][j];
    }}
    for(i=0;i<npsr;i++){
        for(j=0;j<nz2;j++){
            datatp[j][i]=tauppanel[i][j];
    }}
    datawrite(recoverdata,nz2,nx,"datarecover.bin");
    //datawrite(recoverdata-=basedata,nz2,nx,"dataerr.bin");
    datawrite(datatp,nz2,npsr,"datatp.bin");

//float **basedata(newfmat(nz2,nx)),\
        **recoverdata(newfmat(nz2,nx)),\
        **datatp(newfmat(nz2,npsr));
    matdelete(trace,nx);
    matdelete(recover,nx);
    matdelete(tauppanel,npsr);
    delete [] coor;

    return 0;
}


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

