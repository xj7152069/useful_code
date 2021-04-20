//预编译
#ifndef XJ.C++_H_H
#define XJ.C++_H_H

/*@Xiang Jian, Used for Geophycise Research
    #include "./xj.c++.include/wave2D.h"
    #include "./xj.c++.include/RTMangle2D.h"
    #include "./xj.c++.include/my_armadillo.h"
*/

#include "./xj.c++.include/mat.hpp"
#include "./xj.c++.include/my_armadillo.hpp"
#include "./xj.c++.include/wave2D.hpp"
//#include "./xj.ray.f90/caltime.hpp"
#include "./xj.c++.include/radon_beamforming.hpp"

/**************************Extern**************************************
//建议还是在指针数组所在的语句块中释放内存;
//若要调用以下函数释放空间，应在原语句块中再次将指针空置。
template<typename T1> void matdelete(T1 **mat, int x1);
template<typename T1> void matdelete(T1 ***mat, int x1, int x2);
//数组定义原则：越有可能作为整体被调用的越往后；
//因为对于如三维指针***p而言，*p和**p都是地址数组。
float** newfmat(int x1, int x2);
float*** newfmat(int x1, int x2, int x3);
void wave2Dtest(int Z, int X);
template <typename T1> void matprint(T1 *mat1, int n);
template <typename T1> void matprint(T1 **mat1, int nz, int nx);
template <typename T1> void matprint(T1 ***mat1, int n1, int n2, int n3);
template <typename T1> char* numtostr(T1 n, int s);
float wavelet01(int k, float DT, float hz=35.0, int delay=70);
template <typename T1> void wavelet01(T1 *w, int N, float DT, float hz=35.0, int delay=70);
float wavelet02(int k, float DT, float hz=35.0, int delay=70);
template <typename T1> void wavelet02(T1 *w, int N, float DT, float hz=35.0, int delay=70);
template <typename T1, typename T2> void matdec(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename T1, typename T2> void matadd(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> void matadd(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename T1, typename T2> void matcopy(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> void matcopy(T1 ***mat1, T2 n, int n1, int n2, int n3);
template <typename T1, typename T2> void matcopy(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename T1, typename T2> void matmul(T1 **mat1, T2 n, int nz, int nx);
template <typename T1, typename T2> void matmul(T1 **mat1, T2 **mat2, int nz, int nx);
template <typename TT> void dataread(TT **data_mat, int nz, int nx, const char * filename);
template <typename TT> void dataread(TT **data_mat, int nz, int nx, ifstream &inf);
template <typename TT> void datawrite(TT **data_mat, int nz, int nx, const char *filename);
template <typename TT> void datawrite(TT **data_mat, int nz, int nx, ofstream &outf);
void matsmooth(float **mat1, float **mat2, int x1, int x2);
template <typename T1> float* hilbert1D(T1 *s, int n, float dt);

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
***************************************************************************************/



/**************************Extern**************************************
template <typename T2> void matcopy(T2 **mat2, fmat &data_mat, int nz, int nx);
template <typename TT> fmat matcopy(TT **mat, int nz, int nx);
fmat matmul(fmat mat1, fmat mat2, int nz, int nx);
template <typename T1> fmat matmul(fmat mat1, T1 n, int nz, int nx);
void dataread(fmat & data_mat, int nz, int nx, const char * filename);
void dataread(fmat & data_mat, int nz, int nx, ifstream &inf);
void datawrite(fmat & data_mat, int nz, int nx, char const *filename);
void datawrite(fmat & data_mat, int nz, int nx, ofstream &outf);
***************************************************************************************/




#endif






