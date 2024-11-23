
#include<xjc.h>

int main(int argc , char *argv[])
{
    string filein;
    filein=argv[1];
    std::streampos filesize=getFileSize(filein.c_str())/sizeof(float);
    cout<<"File size is: "<<filesize<<endl;
    fmat data1=dataread(filesize,1,filein.c_str());
    float s=atof(argv[2]);
    data1*=s;    
    datawrite(data1,argv[3]);

    return 0;
}
