/*
    Read header keyword of seismic data in su-format.
    Clear up in 2022.11.05, by Xiang Jian, WPI.

*/
#define ARMA_DONT_USE_BLAS //!!!!
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "../include/xjc.h"

using namespace std;
using namespace arma;
void radoninv(struct linerradon3d & par, bool recoverdata);
/*
argv[1]: input file name
argv[2]: output file name
argv[3]: name of head-key:
argv[4]: begin of read segy-trace
argv[5]: end of read segy-trace
argv[6]: nz simple num
argv[7]: trace readgep
data-format: 1-.su
*/
int main(int argc , char *argv[])
{
    int i,j,k,usecopy(0),i1,j1,k1,suformat(0);
    string infile;
    int nz(atoi(argv[6])),nx(1),ny(1),nf(0),obn1,obn2,obn;
    float dt;
    float **trace,**tauppanel,*coor;
    char file[399],file1[199],filekey[99],file3[99];
    char fileout[399],fileout1[99];
    ofstream outfhead;
    ifstream inf1,inf2,inf3;
    int num(0),begtrace(atoi(argv[4])),\
        readgep(atoi(argv[7]));
    float endtrace(atof(argv[5]));
    //suformat=atoi(argv[8]);
    suformat=1;

    file1[0]='\0';
    strcat(file1,argv[1]);
    fileout[0]='\0';
    strcat(fileout,argv[2]);

    //strcat(fileout,fileout1);
    strcat(fileout,".head");
    strcat(fileout,argv[3]);
    
    inf1.open(fileout);
    if(inf1){
        inf1.close();
        outfhead.open(fileout,ios::app);
    }else{
        outfhead.open(fileout);
    }
    
    file[0]='\0';
    filekey[0]='\0';
    strcat(filekey,argv[3]);
    strcat(file,file1);

    segyhead head;
    head.filename[0]='\0';
    strcat(head.filename,file);
    if(suformat==1)
        segyhead_open(head,true);
    else
        segyhead_open(head,false);

if(head.infile){
    if(suformat==1){
        head.endian='l';
        head.isibm=false;
    }else{
        segyhead_endianget(head);
    }
    for(i=0;i<begtrace;i++){
        head.infile.seekg(240,ios::cur);
        head.infile.seekg(nz*nx*sizeof(float),ios::cur);
    }

    fmat a(1,1);
    for(i=0;i<endtrace;i++)
    {
        if(head.infile.eof()){
            cout<<i<<endl;
            break;
        }
        int test;
        {
            //head.dataraw=segyhead_readonetrace_tofmat(head,head.data);
            //nz=head.data.n_rows;
            head.infile.read((char *)(&head.head2), sizeof(head.head2));
            if(head.endian=='b'){
                nz=int(getendianchange(head.head2.ns));}
            else{
                nz=int(head.head2.ns);}
            head.infile.seekg(nz*sizeof(float),ios::cur);
            if(head.endian=='b' && i%readgep==0){
                //float gx=getendianchange(head.head2.gx);
                //float gy=getendianchange(head.head2.gy);
                float keynum=getendianchange(getSuHeadKey(head.head2,filekey));
                outfhead.write((char *)(&keynum), sizeof(keynum));
            }
            else if(i%readgep==0){
                //float gx=(head.head2.gx);
                //float gy=(head.head2.gy);
                float keynum=(getSuHeadKey(head.head2,filekey));
                outfhead.write((char *)(&keynum), sizeof(keynum));
            }
        }
    head.infile.read((char *)(&test), sizeof(test));
    if(head.infile.eof())
    {
        cout<<i<<endl;
        break;}
    else{
        head.infile.seekg(-1*sizeof(test),ios::cur);
    }
    }
    head.infile.close();
}
    outfhead.close();
    return 0;
}
