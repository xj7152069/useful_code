/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
using namespace std;

#include "xjc.h"

float Blackman(int n, int N);
fmat windowsZ(int Z, int X, float k);
fmat fmatsmooth(fmat mat2, int x1, int x2, int k);

int main ()
{
    int i,j,k,sx,source_beg(50),source_gep(10),source_end(350);
    int x1(200),x2(401),x3(401*100);
    char file1[99];
    float **fa,**da,**imge;
    float **far,**dar,**imger;
    fa=newfmat(x1,x3);
    da=newfmat(x1,x3);
    imge=newfmat(x1,x2);
    far=newfmat(x1,x3);
    dar=newfmat(x1,x3);
    imger=newfmat(x1,x2);
    matcopy(fa,0.0,x1,x3);
    matcopy(da,0.0,x1,x3);
    matcopy(imge,0.0,x1,x2);

    sx=source_beg-source_gep;
while(sx<source_end)
{
    sx+=source_gep;
    cout<<sx<<endl;
    file1[0]='\0';
    strcat(file1,"./data/fa.angle.cig.gather");
    strcat(file1,numtostr(sx,5));
    dataread(far,x1,x3,file1);
    matadd(fa,far,x1,x3);

    file1[0]='\0';
    strcat(file1,"./data/da.angle.cig.gather");
    strcat(file1,numtostr(sx,5));
    dataread(dar,x1,x3,file1);
    matadd(da,dar,x1,x3);

    file1[0]='\0';
    strcat(file1,"./data/fa.angle.imag");
    strcat(file1,numtostr(sx,5));
    dataread(imger,x1,x2,file1);
    matadd(imge,imger,x1,x2);
}

    datawrite(fa,x1,x3,"./data/fa.angle.cig.gather");
    datawrite(da,x1,x3,"./data/da.angle.cig.gather");
    datawrite(imge,x1,x2,"./data/fa.angle.imge");
    cout<<"finshed"<<endl;
    return 0;
}


