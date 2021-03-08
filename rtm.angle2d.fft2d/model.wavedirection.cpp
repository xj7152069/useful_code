
//输出地下速度模型

#include <iostream>
#include <stdio.h>
#include <fstream>
#include<iomanip>
#include<cmath>
using namespace std;
#include<xj.c++.h>
int main()
{
float **p, **p2;
p=newfmat(1001,200);
p2=newfmat(1001,200);
double pi=3.1415926;
int i,j,k,n,i1,j1;
   
for(i=0;i<401;i++)
  {for(j=0;j<200;j++)
    {  	
       p[i][j]=3000;
        if (j>80)
            p[i][j]=3500;
        
        if ((j+60)*(j+60)+(i-200)*(i-200)>40000 && j>110)
            p[i][j]=4000;
        
        
        if((j-110)*(j-110)+(i-160)*(i-160)<=4)
            p[i][j]=4000;
     }
   }
for(i=0;i<1001;i++)
  {for(j=0;j<200;j++)
    {  	
        if( i>=300 && i<=700)
        {p2[i][j]=p[i-300][j];}
        if( i<300 || i>700)
        {p2[i][j]=p[0][j];}
    }
  }
     
ofstream outf;
outf.open("model.dat");
 for(j=0;j<1001;j++)
   {for(i=0;i<200;i++)
      {
       outf.write((char*)&p2[j][i],sizeof(p2[j][i]));
       }
    }
outf.close();
/*
cout<<"input the smooth times: ";
cin>>n;
for(k=0;k<n;k++)
{
for(i=1;i<400;i++)
  {for(j=1;j<199;j++)
    {  	
    p2[i][j]=0.2*(p[i][j]+p[i+1][j+1]+p[i-1][j-1]+p[i+1][j-1]+p[i-1][j+1]);
    }
  }
for(i=0;i<401;i++)
  {for(j=0;j<200;j++)
    {  	
    p[i][j]=p2[i][j];
    }
  }
}

outf.open("smooth_model.dat");
 for(j=0;j<401;j++)
   {for(i=0;i<200;i++)
      {
       outf.write((char*)&p[j][i],sizeof(p[j][i]));
       }
    }
outf.close();   
*/
cout<<"finished"<<endl;

return (0);
}

