/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
using namespace std;
#include <xj.c++.h>
 
int main ()
{
   int i,j,i1,j1,sx1,sx2,NZ(200),NX(501),NA(37);
   float ***corangle;
   fmat imge(NZ,NX),a1(NZ,NA),a2(NZ,NA),la1(NZ,1),la2(NZ,1);
   ifstream infa1,infa2;
   char file[99],filea1[99],filea2[99];
   corangle=newfmat(NX,NZ,NA);
   matcopy(corangle,0.0,NX,NZ,NA);
   file[0]='\0';
   strcat(file,"./data/recv.angle.cig.gather");
   imge.fill(0.0);

for(j=250;j<260;j+=10)
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
    cout<<"ok"<<endl;
    return 0;
}
