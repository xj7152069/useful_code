/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
#include <xj.c++.h>
using namespace std;
 
int main ()
{
   int i,j,k,nz(250),nx(501),nray(nz*2+nx*2);
   float dx,dz,**d;
   d=newfmat(1,3);
   fmat invm(nz,nx);
   fmat ray(nz,nx);
   fmat A(nray,nz*nx);
   fmat m(nz*nx,1);
   fmat m2(nz,nx);
   fmat t(nray,1);
   //fmat I(nz*nx,nz*nx);
   //I.diag(0.001);
   ifstream inf1,inf2;
   inf1.open("ray.bin");
   inf2.open("pathmat.bin");
   
   for(k=0;k<nray;k++)
   {
      ray=dataread(nz,nx,inf2);
      for(i=0;i<nz;i++)
      {
         for(j=0;j<nx;j++)
         {
            A(k,i*nz+j)=ray(i,j);
         }
      }
      while(1)
      {
         dataread(d,1,3,inf1);
         if(d[0][0]==3.0)
         {
            t(k,0)=d[0][1];
            break;
         }
      }
   }
   cout<<"ok"<<endl;
   m=inv(A.st()*A)*A.st()*t;
   for(i=0;i<nz;i++)
      {
         for(j=0;j<nx;j++)
         {
            m2(i,j)=m(i*nz+j,0);
         }
      }
   datawrite(m,nz,nx,"invmodel.bin");

   return 0;
}
