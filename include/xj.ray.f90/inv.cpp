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
   int i,j,k,nz(60),nx(80),ns(2*(nz-1+nx-1)),nray((nz*2+nx*2-4)*ns);
   float dx,dz,**d;
   d=newfmat(1,3);
   fmat invm(nz,nx);
   fmat ray(nz,nx);
   fmat A(nray,nz*nx);
   fmat m(nz*nx,1);
   fmat m2(nz,nx);
   fmat t(nray,1);
   fmat I(nz*nx,nz*nx);
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
            A(k,i*nx+j)=ray(i,j);
         }
      }
      while(1)
      {
         dataread(d,1,3,inf1);
         if(d[0][0]==3.0)
         {
            t(k,0)=d[0][1];
            if(d[0][1]<=0)
            {
               t(k,0)=0.0;
               A.row(k).fill(0.0);
               cout<<"error? "<<k<<endl;
            }
            //cout<<"ok"<<k<<endl;
            break;
         }
      }
   }
   cout<<"ok"<<endl;
   I=A.st()*A;
   for(k=0;k<nx*nz;k++)
   {
      I(k,k)+=0.001;
   }
   m=inv(I)*A.st()*t;
   cout<<"ok"<<endl;
   for(i=0;i<nz;i++)
      {
         for(j=0;j<nx;j++)
         {
            m2(i,j)=1.0/m(i*nx+j,0);
         }
      }
   datawrite(m2,nz,nx,"invmodel.bin");

   return 0;
}
