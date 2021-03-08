/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/
#include <armadillo>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace arma;
using namespace std;

#include <xj.c++.h>
 
int main ()
{
   int NZ(200),NX(100),NF(1024),Nf(500),Np(2001),ncig(501);
   float dx(10.0),dp(0.000001),p1(-0.001);
   float DT(0.001);
   fmat forA(NX,Np);
   fmat imge(NF,NX);
   fmat data(NZ,NX);
   cx_fmat dataTP(NF,Np);
   cx_fmat datarebuild(NF,NX);
   ifstream inf;
   fmat realtp(NF,Np);
   fmat realrebuild(NF,NX);
   inf.open("./data701/fa.angle.cig.gather250.0");

   int i,j,k;
   k=210;
   inf.seekg(k*NX*NZ*sizeof(float),ios::beg);
   data=dataread(NZ,NX,inf);
   forA=CreatePmat(NX, Np, dx, dp, p1);
   dataTP=getTP(data,forA,Np,NX,NF,Nf,DT);
   realtp=real(dataTP);
   datawrite(realtp,NF,Np,"tp.bin");
   for(i=0;i<Np;i++)
   {
      if(abs(i-500)>=20)
      {dataTP.col(i).fill(0.0);}
   }
   datarebuild=getRebuildData(dataTP,forA,Np,NX,NF,Nf,DT);
   realrebuild=real(datarebuild);
   datawrite(realrebuild,NF,NX,"rebuild.bin");
 
   return 0;
}
