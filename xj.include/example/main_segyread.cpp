/*
    wiener filter: suit for 2d or 3d data
*/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "./include/xjc.h"
using namespace std;
using namespace arma;


int main()
{
    ifstream inf;
    inf.open("./result/data3.wiener.bin");
    char file[99];
    file[0]='\0';
    strcat(file,"/root/pai-test/project3.sgy");
    segyhead head;
    head.filename[0]='\0';
    strcat(head.filename,file);
    segyhead_open(head);
    segyhead_endianget(head);
    int nline(51),keyw(1),num(0);
/*
    for(int k=0;k<1;k++)
    {
        segyhead_readonetrace_tofmat(head,data);
    }*/
    
    ofstream outf;
    outf.open("./result/project3.wiener.segy");

//outf.write((char*)&data_mat(i,j),sizeof(float));
outf.write((char*)(&head.head0), sizeof(head.head0));
outf.write((char*)(&head.head1), sizeof(head.head1));
while(!head.infile.eof())
{
    segyhead_readonetrace_tofmat(head,head.data);
    dataread(head.data,head.data.n_rows,1,inf);
    fmatendianchange(head.data,head.data.n_rows,1);
    if(!head.infile.eof())
    {
    outf.write((char*)(&head.head2), sizeof(head.head2));
    datawrite(head.data,head.data.n_rows,1,outf);
    }
}
cout<<head.infile.eof()<<"|"<<inf.eof()<<endl;
    inf.close();
    head.infile.close();
    outf.close();

    return 0;
}












