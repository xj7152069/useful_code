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

#include "radonbeamforming.h"
#include <xj.c++.h>
 
int main ()
{
   int NZ(200),NX(181),NF(1024),Nf(500),Np(1001);
   float dx(10.0),dp(0.000001),p1(-0.0005);
   float DT(0.001);
   fmat forA(NX,Np);
   fmat data(NZ,NX);
   cx_fmat dataTP(NF,Np);
   cx_fmat datarebuild(NF,NX);
   data=dataread(NZ,NX,"cdp260.000000");
   int i;
   for(i=0;i<30;i++)
   {
      data.row(i).fill(0.0);
   }
   forA=CreatePmat(NX, Np, dx, dp, p1);
   dataTP=getTP(data,forA,Np,NX,NF,Nf,DT);
   for(i=0;i<Np;i++)
   {
      if(abs(i-500)>=20)
      {dataTP.col(i).fill(0.0);}
   }
   datarebuild=getRebuildData(dataTP,forA,Np,NX,NF,Nf,DT);

   fmat realtp(NF,Np);
   fmat realrebuild(NF,NX);
   realtp=real(dataTP);
   realrebuild=real(datarebuild);
   datawrite(realtp,NF,Np,"tp.bin");
   datawrite(realrebuild,NF,NX,"rebuild.bin");
 
   return 0;
}
