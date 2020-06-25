/**************用�?�异值分解将线性信号分离v2.0*****************/
/*注释�?
    主�?�功能：
        1.按编号�?�取Dip-angle二进制文件；
        2.将信号转换到频率-空间域，以每�?单�?�片构建Hankel矩阵�?
        3.用SVD变换进�?��?�异值分解，选则需要的奇异值重构“Hankel矩阵”；
        4.将重构的矩阵反�?��?�平均，使其成为标准的反对�?�矩阵；
        5.反变换回时间-空间域信号，输出结果�?
*/
/************************************************************/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
using namespace std;
using namespace arma;

#include <xj.c++.h>
extern cx_fmat mathankel(cx_fmat &mat1, int N, int k);
extern void ihankelaverage(cx_fmat &hankel, int N, cx_fmat &ihankel, int k2);

int main(int argc, char** argv)
{
    int x1(1024),x2(200),n2,i,j,k;
    int fx1(1024),fx2,Rank1(0),Rank2(1);
    char file[99];

    cout<<"Rank1=";cin>>Rank1;
    cout<<"Rank2=";cin>>Rank2;

    file[0]='\0';
    strcat(file,"txz1.bin");
    n2=x2;  //Ҫ��������Ϊ��������Ϊ˫�������
    if(x2%2==0)
    {
        x2=x2-1;
    }
    fx2=x2;

    fmat dataraw(x1,n2);
    dataraw.fill(0.0);
    dataraw=dataread(x1,n2,file);   

    fmat data(x1,x2);
    cx_fmat fftdata(fx1,fx2), ihankel(fx1,fx2);
    for(i=0;i<x2;i++)
    {
        data.col(i)=dataraw.col(i);
        fftdata.col(i)=fft(data.col(i),fx1);
    }

    int nh;
    nh=(x2+1)/2;
    fmat zeromat(nh,nh);
    zeromat.fill(0.0);
    cx_fmat hankel(nh,nh), hankel_svd(nh,nh);
    cx_fmat U;
    fvec s;
    cx_fmat V;
    for(k=0;k<fx1;k++)
    {
        hankel=mathankel(fftdata, nh, k);
        svd(U,s,V,hankel);  //SVD
        hankel_svd.set_imag(zeromat);
        hankel_svd.set_real(zeromat);
        for(i=Rank1;i<Rank2;i++)  //rebuild Hankel
            {
            hankel_svd=hankel_svd+U.col(i)*s(i)*(V.col(i).t());
            }
        if(k==40)
            output_qyzp(s,nh,"qyzp.bin");
        if(k%100==0)
        {
            cout<<k<<endl;
        }
        ihankelaverage(hankel_svd, nh, ihankel, k);
    }

    cx_fmat ifftdata(fx1,fx2);
    for(i=0;i<fx2;i++)
    {
        ifftdata.col(i)=ifft(ihankel.col(i),fx1);
    }
    for(i=0;i<x1;i++)
    {
        for(j=0;j<x2;j++)
        {
            dataraw(i,j)=ifftdata(i,j).real();
        }
    }
    datawrite(dataraw,x1,n2,"out.bin");

    return 0;
}
















