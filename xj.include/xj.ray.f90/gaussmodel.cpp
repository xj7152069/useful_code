/*********(version 1.0)***********/
/*注释�??
    C++程序模板�??

*/
/********************************/

#include <iostream>
#include <xj.c++.include/mat.h>
using namespace std;
void addrandonnoise(float** data,int n1, int n2, float N);
void addgaussonnoise(float** data,int n1, int n2, float N, float p);
void addgaussobject(float** data,int nz,int nx,int wz,int wx,float N,float p);

float** denoise_nearwindows(float** data, int nz, int nx,\
   int nwz, int nwx, int nwz2, int nwx2, float pdp);

float** denoise_bilateral(float** data, int nz, int nx,\
   int nwz, int nwx, float pdp, float dz, float dx, float dp);

float** denoise_near(float** data, int nz, int nx,\
   int nwz, int nwx, float pdp);

float** denoise_gauss_normal(float** data, int nz, int nx,\
   int nwz, int nwx, float p);
 
int main ()
{
   int nz(60),nx(80);
   float **data;
   data=newfmat(nz,nx);
   matcopy(data,3000.0,nz,nx);
   addgaussobject(data,nz,nx,30,40,800,300);
   //addgaussobject(data,nz,nx,25,60,1000,80);
   //addgaussobject(data,nz,nx,35,45,-600,100);
   datawrite(data,nz,nx,"model.dat");
   return 0;
}
void addgaussobject(float** data,int nz,int nx,int wz,int wx,float N,float p)
{
   int i,j,k;
   float w;

   for(i=0;i<nz;i++)
   {
      for(j=0;j<nx;j++)
      {
         w=N*exp(-((i-wz)*(i-wz)+(j-wx)*(j-wx))/2/p);
         data[i][j]+=w;
      }
   }

}

float** denoise_nearwindows(float** data, int nz, int nx,\
   int nwz, int nwx, int nwz2, int nwx2, float pdp)
{
   int i,j,k,i1,j1,i2,j2,wz,wx,wz2,wx2,js;
   float W,**w,**data2,p(0),sumpower(0);

   if(nwz%2==0)
   {nwz=nwz+1;}
   if(nwx%2==0)
   {nwx=nwx+1;}
   if(nwz2%2==0)
   {nwz2=nwz2+1;}
   if(nwx2%2==0)
   {nwx2=nwx2+1;}
   w=newfmat(nwz,nwx);
   data2=newfmat(nz,nx);
   matcopy(data2,0.0,nz,nx);
   wz=(nwz-1)/2;
   wx=(nwx-1)/2;
   wz2=(nwz2-1)/2;
   wx2=(nwx2-1)/2;
   W=0.0;

   //datawrite(w,nwz,nwx,"gauss.normal.windows.bin");
   js=0;
   for(i=wz+wz2;i<nz-wz-wz2;i++)
   {
      for(j=wx+wx2;j<nx-wx-wx2;j++)
      {
         p=0;
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               sumpower=0;
               for(i2=0;i2<nwz2;i2++)
               {
                  for(j2=0;j2<nwx2;j2++)
                  {
                     sumpower+=((data[i-wz+i1-wz2+i2][j-wx+j1-wx2+j2]-data[i-wz2+i2][j-wx2+j2])*\
                           (data[i-wz+i1-wz2+i2][j-wx+j1-wx2+j2]-data[i-wz2+i2][j-wx2+j2]));
                  }
               }
               if(p<(sumpower/pdp))
               {p=(sumpower/pdp);}
            }
         }
         
         W=0.0;
         if(p<1)
         {p=1;}
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               sumpower=0;
               for(i2=0;i2<nwz2;i2++)
               {
                  for(j2=0;j2<nwx2;j2++)
                  {
                     sumpower+=((data[i-wz+i1-wz2+i2][j-wx+j1-wx2+j2]-data[i-wz2+i2][j-wx2+j2])*\
                           (data[i-wz+i1-wz2+i2][j-wx+j1-wx2+j2]-data[i-wz2+i2][j-wx2+j2]));
                  }
               }
               w[i1][j1]=exp(-sumpower/p);
               W=W+w[i1][j1];
            }
         }
         js++;
         if(js%100000==0)
         {cout<<"denoise_nearwindows is running:"<<js<<"  p= "<<p<<endl;}
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               data2[i][j]+=w[i1][j1]*data[i-wz+i1][j-wx+j1];
            }
         }
         data2[i][j]=data2[i][j]/W;
      }
   }
   return data2;
}


float** denoise_bilateral(float** data, int nz, int nx,\
   int nwz, int nwx, float pdp, float dz, float dx, float dp)
{
   int i,j,k,i1,j1,wz,wx,js;
   float W,**wd,**w,**data2,p(0);

   if(nwz%2==0)
   {nwz=nwz+1;}
   if(nwx%2==0)
   {nwx=nwx+1;}
   w=newfmat(nwz,nwx);
   wd=newfmat(nwz,nwx);
   data2=newfmat(nz,nx);
   matcopy(data2,0.0,nz,nx);
   wz=(nwz-1)/2;
   wx=(nwx-1)/2;
   W=0.0;

   for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               wd[i1][j1]=exp(-((i1-wz)*(i1-wz)*dz*dz+(j1-wx)*(j1-wx)*dx*dx)/dp);
            }
         }
   datawrite(wd,nwz,nwx,"bilateral.distance.windows.bin");
   js=0;
   for(i=wz;i<nz-wz;i++)
   {
      for(j=wx;j<nx-wx;j++)
      {
         p=0;
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               if(p<((data[i-wz+i1][j-wx+j1]-data[i][j])*\
                     (data[i-wz+i1][j-wx+j1]-data[i][j]))/pdp)
                  {p=((data[i-wz+i1][j-wx+j1]-data[i][j])*\
                     (data[i-wz+i1][j-wx+j1]-data[i][j]))/pdp;}
            }
         }
         
         W=0.0;
         if(p<1)
         {p=1;}
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               w[i1][j1]=exp(-(data[i-wz+i1][j-wx+j1]-data[i][j])*\
                     (data[i-wz+i1][j-wx+j1]-data[i][j])/p);
               W=W+w[i1][j1]*wd[i1][j1];
            }
         }

         js++;
         if(js%100000==0)
         {cout<<"denoise_bilateral is running:"<<js<<"  p= "<<p<<endl;}

         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               data2[i][j]+=w[i1][j1]*wd[i1][j1]*data[i-wz+i1][j-wx+j1];
            }
         }
         data2[i][j]=data2[i][j]/W;
      }
   }
   return data2;
}

float** denoise_near(float** data, int nz, int nx, int nwz, int nwx, float pdp)
{
   int i,j,k,i1,j1,wz,wx,js;
   float W,**w,**data2,p(0);

   if(nwz%2==0)
   {nwz=nwz+1;}
   if(nwx%2==0)
   {nwx=nwx+1;}
   w=newfmat(nwz,nwx);
   data2=newfmat(nz,nx);
   matcopy(data2,0.0,nz,nx);
   wz=(nwz-1)/2;
   wx=(nwx-1)/2;
   W=0.0;

   //datawrite(w,nwz,nwx,"gauss.normal.windows.bin");
   js=0;
   for(i=wz;i<nz-wz;i++)
   {
      for(j=wx;j<nx-wx;j++)
      {
         p=0;
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               if(p<((data[i-wz+i1][j-wx+j1]-data[i][j])*\
                     (data[i-wz+i1][j-wx+j1]-data[i][j]))/pdp)
                  {p=((data[i-wz+i1][j-wx+j1]-data[i][j])*\
                     (data[i-wz+i1][j-wx+j1]-data[i][j]))/pdp;}
            }
         }
         
         W=0.0;
         if(p<1)
         {p=1;}
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               w[i1][j1]=exp(-(data[i-wz+i1][j-wx+j1]-data[i][j])*\
                     (data[i-wz+i1][j-wx+j1]-data[i][j])/p);
               W=W+w[i1][j1];
            }
         }
         js++;
         if(js%100000==0)
         {cout<<"denoise_near is running:"<<js<<"  p= "<<p<<endl;}
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               data2[i][j]+=w[i1][j1]*data[i-wz+i1][j-wx+j1];
            }
         }
         data2[i][j]=data2[i][j]/W;
      }
   }
   return data2;
}


float** denoise_gauss_normal(float** data, int nz, int nx, int nwz, int nwx, float p)
{
   int i,j,k,i1,j1,wz,wx;
   float W,**w,**data2;

   if(nwz%2==0)
   {nwz=nwz+1;}
   if(nwx%2==0)
   {nwx=nwx+1;}
   w=newfmat(nwz,nwx);
   data2=newfmat(nz,nx);
   matcopy(data2,0.0,nz,nx);
   wz=(nwz-1)/2;
   wx=(nwx-1)/2;
   W=0.0;

   for(i=0;i<nwz;i++)
   {
      for(j=0;j<nwx;j++)
      {
         w[i][j]=exp(-((i-wz)*(i-wz)+(j-wx)*(j-wx))/2/p);
         W=W+w[i][j];
      }
   }
   datawrite(w,nwz,nwx,"gauss.normal.windows.bin");
   for(i=wz;i<nz-wz;i++)
   {
      for(j=wx;j<nx-wx;j++)
      {
         for(i1=0;i1<nwz;i1++)
         {
            for(j1=0;j1<nwx;j1++)
            {
               data2[i][j]+=w[i1][j1]*data[i-wz+i1][j-wx+j1];
            }
         }
         data2[i][j]=data2[i][j]/W;
      }
   }
   return data2;
}

void addrandonnoise(float** data,int n1, int n2, float N)
{
   int i,j,k;
   float n;
   for(i=0;i<n1;i++)
   {
      for(j=0;j<n2;j++)
      {
         k=rand();
         n=(float(k%10)*0.1*N-N/2);
         data[i][j]+=n;
      }
   }
}

void addgaussonnoise(float** data,int n1, int n2, float N, float p)
{
   int i,j,k;
   float n;
   for(i=0;i<n1;i++)
   {
      for(j=0;j<n2;j++)
      {
         k=rand();
         n=(float(k%10)*0.1-0.5);
         n=N*exp(-(n*n/2.0/p));
         n=n-N/2.0;
         data[i][j]+=n;
      }
   }
}
