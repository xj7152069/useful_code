/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
using namespace std;
#include <xj.c++.h>
float fmatpxxg(fmat & a1, int n, int minn);
float fmatpxxg2(fmat & a1, int n, int minn);
int main ()
{
   int i,j,i1,j1,sx1,sx2,NZ(200),NX(401),NA(37),minn(10),dsx(40);
   float ***corangle;
   fmat imge(NZ,NX),a1(NZ,NA),a2(NZ,NA),la1(NZ,1),la2(NZ,1),la3(1,NA);
   ifstream infa1,infa2;
   char file[99],filea1[99],filea2[99];
   corangle=newfmat(NX,NZ,NA);
   matcopy(corangle,0.0,NX,NZ,NA);
   file[0]='\0';
   strcat(file,"./data/recv.angle.cig.gather");
   imge.fill(0.0);

for(j=50;j<350;j+=10)
{
    cout<<"now is running: "<<j<<endl;
    filea1[0]='\0';
    strcat(filea1,file);
    strcat(filea1,numtostr(j,5));
    infa1.open(filea1);
    for(i=0;i<NX;i++)
    {
        infa1.seekg(NZ*NA*i*sizeof(float),ios::beg);
        a1=dataread(NZ,NA,infa1);
        a2=dataread(NZ,NA,infa1);
        if(abs(j-i)>dsx)
        {
            for(i1=0;i1<NZ;i1++)
            {
                la3.row(0)=a1.row(i1);
                imge(i1,i)+=fmatpxxg2(la3,NA,minn);
            }
        }
        /*
        for(i1=0;i1<NA;i1++)
        {
            la1.col(0)=a1.col(i1);
            for(j1=0;j1<NA;j1++)
            {
                if(abs(j1-i1)>=10)
                {
                    la2.col(0)=a2.col(j1);
                    la2=matmul(la2,la1,NZ,1);
                    imge.col(i)=imge.col(i)+la2.col(0);
                }
            }
        }
        */
    }
    infa1.close();
}
    datawrite(imge,NZ,NX,"diffs.imge.bin");
    float **mat1,**mat0;
    mat0=newfmat(NZ,NX);
    mat1=newfmat(NZ,NX);
    matcopy(mat0,imge,NZ,NX);
    Laplace(mat1,mat0,NZ,NX);
    datawrite(mat1,NZ,NX,"diffs.imge.la.bin");
    Laplace(mat1,mat1,NZ,NX);
    datawrite(mat1,NZ,NX,"diffs.imge.la2.bin");
    cout<<"ok"<<endl;
    return 0;
}

float fmatpxxg2(fmat & a1, int n, int minn)
{
    fmat a2(1,n);
    int j,i,k;
    float mina(99999),a,b,b2;
    a2=matmul(a1,a1,1,n);
    mina=99999;
    for(i=0;i<n-minn;i++)
    {
        a=0;b=0;
        for(j=i;j<i+minn;j++)
        {
            a+=a2(0,j);
            b+=a1(0,j);
        }
        if(a<mina)
        {
            mina=a;
            b2=b;
        }
    }

    a=0;b=b2;
    for(i=0;i<n;i++)
    {
        a+=a1(0,i);
    }
    a=a*b;
    return a;
}

float fmatpxxg(fmat & a1, int n, int minn)
{
    fmat a2(1,n);
    int j,i,k;
    float mina(99999),a,b;
    a2=matmul(a1,a1,1,n);
    for(i=0;i<n;i++)
    {
        mina=99999;
        for(j=i;j<n;j++)
        {
            if(a2(0,j)<mina)
            {
                mina=a2(0,j);
                k=j;
            }
        }
        a=a2(0,i);
        a2(0,i)=mina;
        a2(0,k)=a;
        a=a1(0,i);
        a1(0,i)=a1(0,k);
        a1(0,k)=a;
    }

    a=0;b=0;
    for(i=0;i<minn;i++)
    {
        a+=a1(0,i);
    }
    for(i=minn;i<n;i++)
    {
        b+=a1(0,i);
    }
    a=a*b;
    return a;
}

