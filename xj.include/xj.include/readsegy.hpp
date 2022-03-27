#ifndef READSEGY_HPP
#define READSEGY_HPP

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "../xjc.h"

using namespace std;

unsigned int ibm2ieee (const unsigned int num);
void segyhead_endianget(segyhead & head);
void segyhead_initialize(segyhead & head);
template <typename T1> inline void endianchange(T1 & a);
template <typename T1> inline float getendianchange(T1 a);
inline float floatendianchange(float a);
void segyhead_open(segyhead & head);
void fmatendianchange(fmat & data,int n1,int n2);
void segyhead_readonetrace_tofmat(segyhead & head, fmat & trace);
//void segyhead_readoneline_tofmat(segyhead & head, fmat & data);

///////////////////////////////////////////////////////////////////
void dataread_fromsegy_to3dblockone(segyhead & head, int nline, int keyw)
{
    //keyw==1, use cdp; keyw==2, use trace;
    int i,j,k,nz(3501),ny(15),cdp0,cdp1,trac0,trac1,sx0,sx1;
    ny=nline;
    fmat numnx(ny,1);
    ofstream outf,outf2;
    outf.open("./data1.bin",ios::trunc);
    int nx(0),maxnx(-1),nynum(0),begnx(0);
    ifstream infhead;
    while(1){
        if(head.begnx==0){
            segyhead_readonetrace_tofmat(head,head.data);
            nz=head.data.n_rows;
	        nx++;head.begnx=1;    
            datawrite(head.data,nz,1,outf);
            if(head.endian=='l'){
                trac0=getendianchange(head.head2.cdp);
                cdp0=getendianchange(head.head2.f1);
		        sx0=getendianchange(head.head2.d2);
                }
            else{
                trac0=head.head2.cdp;cdp0=head.head2.f1;
		        sx0=head.head2.d2;
                }
            }
        else if(head.begnx==1){
            segyhead_readonetrace_tofmat(head,head.data);
            nz=head.data.n_rows;
            if(head.endian=='l'){
                trac1=getendianchange(head.head2.cdp);
                cdp1=getendianchange(head.head2.f1);
		        sx1=getendianchange(head.head2.d2);
                }
            else{
                trac1=head.head2.cdp;cdp1=head.head2.f1;
		        sx1=head.head2.d2;
                }
            //cout<<"cdp="<<cdp1<<"|"<<"trac="<<trac1<<"|"<<sx1<<"|"<<nynum<<endl;
            if(keyw==2)
            {cdp1=trac1;cdp0=trac0;}
            else if(keyw==3)
            {cdp1=sx1;cdp0=sx0;}
            if(cdp1==cdp0){
                nx++;
                datawrite(head.data,nz,1,outf);
            }
            else{
                if(maxnx<nx)
                {maxnx=nx;}
                //cout<<nx<<",";
                numnx(nynum,0)=nx;
                nx=1;
                cdp0=cdp1;
		        sx0=sx1;
		        trac0=trac1;
                nynum++;
                if(nynum>=(ny)){
                    head.begnx=2;
                    break;
                }
                else
                {datawrite(head.data,nz,1,outf);}
            } 
        }
        else if(head.begnx==2){
            nz=head.data.n_rows;
            nx++;head.begnx=1;    
            datawrite(head.data,nz,1,outf);
            if(head.endian=='l'){
                trac0=getendianchange(head.head2.cdp);
                cdp0=getendianchange(head.head2.f1);
	            sx0=getendianchange(head.head2.d2);
            }
            else{
                trac0=head.head2.cdp;cdp0=head.head2.f1;
	            sx0=head.head2.d2;
            }
        }
        if(head.infile.eof())
        {break;}
    }
    //head.infile.close();
    outf.close();
    //datawrite(numnx,ny,1,"numofnx.bin");
/////////////////////////////
    ifstream inf;
    fmat data;
    data.resize(nz,maxnx);
    inf.open("./data1.bin");
    outf.open("./data1.oneblock.bin",ios::trunc);
    outf2.open("./data1.nx.bin",ios::app);
    for(i=0;i<ny;i++)
    {
        data.fill(0.0);
        dataread(data,nz,numnx(i,0),inf);
        datawrite(data,nz,maxnx,outf);
    }
    datawrite(numnx,ny,1,outf2);
    inf.close();
    outf.close();
    outf2.close();
    head.nz=nz,head.nx=maxnx;
    head.datanx.zeros(ny,1);
    head.datanx=numnx;
    //cout<<"nz="<<head.nz<<"|"<<"nx="<<head.nx<<endl;
}

/*
void segyhead_readoneline_tofmat(segyhead & head, fmat & data)
{
    int nz,nx,k;
    nz=getendianchange(head.head1.hns);
    nx=getendianchange(head.head1.nart)+\
        getendianchange(head.head1.ntrpr);
    cout<<nz<<","<<nx<<endl;
    fmat trace;
    data.zeros(nz,nx);
    for(k=0;k<nx;k++)
    {
    if(!head.infile.eof())
    {
        segyhead_readonetrace_tofmat(head,trace);
        data.col(k)=trace.col(0);
    }
    }
    //head.infile.read((char *)(&head.head0), sizeof(head.head0));
    //head.infile.read((char *)(&head.head1), sizeof(head.head1));
}
*/
void segyhead_readonetrace_tofmat(segyhead & head, fmat & trace)
{
    int nz;
    head.infile.read((char *)(&head.head2), sizeof(head.head2));
    if(head.endian=='l')
        nz=int(getendianchange(head.head1.hns));
    else
        nz=int(head.head1.hns);
    trace.zeros(nz,1);

    if(head.isibm){
        int i,j;
        float *float_p;
        unsigned int IBMFloatBytes;
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                IBMFloatBytes=getendianchange(IBMFloatBytes);
                IBMFloatBytes=ibm2ieee(IBMFloatBytes);
                float_p=(float *)(&IBMFloatBytes);
                trace(i,j)=*float_p;
                if(isinf(abs(trace(i,j)))){
                    trace(i,j)=0;
                }
            }
        }
    }
    else if(!head.isibm && head.endian=='l'){
        int i,j;
        float *float_p;
        unsigned long IBMFloatBytes;
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                IBMFloatBytes=getendianchange(IBMFloatBytes);
                float_p=(float *)(&IBMFloatBytes);
                trace(i,j)=*float_p;
            }
        }
    }
    else if(!head.isibm && head.endian=='b'){
        int i,j;
        float *float_p;
        unsigned long IBMFloatBytes;
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                float_p=(float *)(&IBMFloatBytes);
                trace(i,j)=*float_p;
            }
        }
    }
}

void segyhead_open(segyhead & head)
{
    head.infile.open(head.filename,ios::binary);
    if(!head.infile) cout<<"file open error: "<<head.filename<<endl;
    head.infile.read((char *)(&head.head0), sizeof(head.head0));
    head.infile.read((char *)(&head.head1), sizeof(head.head1));
}

void segyhead_initialize(segyhead & head)
{
    head.filename[0] = '\0';
    head.endian = '0';
}

void segyhead_endianget(segyhead & head)
{
    if ((head.head1.format == 1) || (head.head1.format == 2) || \
        (head.head1.format == 3) || (head.head1.format == 4) || \
        (head.head1.format == 5) || (head.head1.format == 6) || \
        (head.head1.format == 7) || (head.head1.format == 8))
    {
        head.endian = 'b'; //大端(123456)，高位字节在前，低位字节在后
        cout<<"the format of the SEGY-data is Big-Endian: "<<\
            head.head1.format<<endl;
        if(head.head1.format == 1)
        {
            head.isibm=true;
            cout<<"the format of the SEGY-data is IBM: "<<\
                head.head1.format<<endl;
        }
        else
        {
            head.isibm=false;
            cout<<"the format of the SEGY-data is IEEE: "<<\
                head.head1.format<<endl;
        }
    }
    else if(1)
    {
        endianchange(head.head1.format);
        if ((head.head1.format == 1) || (head.head1.format == 2) || \
            (head.head1.format == 3) || (head.head1.format == 4) || \
            (head.head1.format == 5) || (head.head1.format == 6) || \
            (head.head1.format == 7) || (head.head1.format == 8))
        {
            head.endian = 'l'; //小端(654321)，与大端相反，需要调整为大端
            cout<<"the format of the SEGY-data is Little-Endian: "<<\
                head.head1.format<<endl;
            if(head.head1.format == 1)
            {
                head.isibm=true;
                cout<<"the format of the SEGY-data is IBM: "<<\
                    head.head1.format<<endl;
            }
            else
            {
                head.isibm=false;
                cout<<"the format of the SEGY-data is IEEE: "<<\
                    head.head1.format<<endl;
            }
        }
        endianchange(head.head1.format);
    } 
    else
    {
        head.endian = '0';
        cout<<"Error:the format of the SEGY-data is not know!!!"<<endl;
    }
}

template <typename T1> 
inline void endianchange(T1 & a)
{
    char *p,t;
    int i,bit_num;
    bit_num=sizeof(a);
    p=(char *)(&a);
    for(i=0;i<bit_num/2;i++)
    {
        t=*(p+i);
        *(p+i)=*(p+bit_num-i-1);
        *(p+bit_num-i-1)=t;
    }
}

void fmatendianchange(fmat & data,int n1,int n2)
{
    int i,j;
    float a;
    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            a=data(i,j);
            data(i,j)=floatendianchange(a);
        }
    }
}

inline float floatendianchange(float a)
{
    typedef union SWAP_UNION{
        float f;
        char  c[4];
        }SWAP_UNION;
    SWAP_UNION d1,d2;
    d1.f=a;
    d2.c[0]=d1.c[3];
    d2.c[1]=d1.c[2];
    d2.c[2]=d1.c[1];
    d2.c[3]=d1.c[0];
    return d2.f;
}

template <typename T1> 
inline float getendianchange(T1 a)
{
    char *p,t;
    int i,bit_num;
    bit_num=sizeof(a);
    p=(char *)(&a);
    for(i=0;i<bit_num/2;i++)
    {
        t=*(p+i);
        *(p+i)=*(p+bit_num-i-1);
        *(p+bit_num-i-1)=t;
    }
    float b;
    b=float(a);
    return b;
}

unsigned int ibm2ieee (const unsigned int num) 
{  
    
    unsigned int ibmCode(num);
    //unsigned int ieeeFloat(0);
    //unsigned int ieeeFloat(0);
if((ibmCode << 1) == 0)
{  
// If it is an IBM floating point number  
// Returns the corresponding IEEE binary encoding  
    return ibmCode;
}
//IBM浮点数： SEEEEEEE MMMMMMMM MMMMMMMM MMMMMMMM
//浮点数值 Value = (-1)^s * M * 16^(E-64)
   unsigned int  signCode(num>>31);//获取符号位, S=0或1， sign=00 00 00 0000000S
   /*
   int sign(1); //正数
   if (signCode==1)
      sign=-1;//负数
*/
   ibmCode= (ibmCode<<1); // 左移移出符号位, ibmCode= EEEEEEEM MMMMMMMM MMMMMMMM MMMMMMM0
 /*
   int exponent((int)(ibmCode >> 25));  // 获取阶数, exponent=00 00 00 0EEEEEEE
 
   unsigned long fraction(ibmCode << 7);  //移出符号位和阶数剩余的部分：尾数部分fraction=MMMMMMMM MMMMMMMM MMMMMMM 00000000
   fraction= (fraction >> 8);  //fraction=00000000 MMMMMMMM MMMMMMMM MMMMMMM
 
    double ratio(pow(2,24)); //使用静态局部变量，避免重复计算
   if ((exponent==0) && (fraction==0))  //00000000 00 00 00 或10000000 00 00 00 都表示0
      ieeeFloat=0;
      
   else //IBMfloat= (sign)*Fraction*16^(exponent-64);
   {
      double   P2(pow(16, exponent-64));
      double   P1(fraction/ratio);
      ieeeFloat=sign*P1*P2;
   }
    return   ieeeFloat;
*/

    unsigned int float_IBM(num);
/*
if((float_IBM << 1) == 0)
{  
// If it is an IBM floating point number  
// Returns the corresponding IEEE binary encoding  
    return float_IBM;
}*/
// Get the S, symbol part of the IBM floating point number  
unsigned int S_IBM_32(signCode<<31);  
// Get the E, exponential part of the IBM floating point number  
unsigned int exponent(ibmCode >> 25);
//int exponent2((int)(ibmCode >> 25));
unsigned int E_IBM_32(exponent << 24);  
// Get the F, decimal part of the IBM floating point number 
unsigned int fraction(ibmCode << 7);  //移出符号位和阶数剩余的部分：尾数部分fraction=MMMMMMMM MMMMMMMM MMMMMMM 00000000
   fraction= (fraction >> 8);  //fraction=00000000 MMMMMMMM MMMMMMMM MMMMMMM
unsigned int F_IBM_32(fraction);  
// Get the F, decimal part of the IBM floating point number  
unsigned int radix(0);  
unsigned int F_IEEE_32(F_IBM_32);  
while (radix <= 3 && F_IEEE_32 < 0x01000000) {  
    radix++;  
    F_IEEE_32 = F_IEEE_32 << 1;  }  
// Put it back in the appropriate position in the IEEE type F section  
F_IEEE_32 = (F_IEEE_32 - 0x01000000)>>1;  
// Get the E, exponential part of the IBM floating point number  
 // Start counting  
// Put it in the correct position  
/*
double ratio(pow(2,24));
double   P2(pow(16, exponent2-64));
double   P1(fraction/ratio);
ieeeFloat=sign*P1*P2;*/
unsigned int E_IEEE_32((((E_IBM_32>>22)-130)-(radix-1))<<23);  
// Whether overflow occurs  
if (E_IEEE_32 > 0x7F000000) {  
    //return S_IBM_32|0x7F000000;  }  
    return S_IBM_32;  }  
if (E_IEEE_32 < 0x10000000) {  
    return S_IBM_32;  } 
else  
    //return ieeeFloat;  
    return S_IBM_32 | E_IEEE_32 | F_IEEE_32;  
}  

////////////////////////////////////////////////
void Bubble_sort(float *data, int n)
{
    int i, j;
    float t;
    for (j = 1; j <= n; j++)
    {
        for (i = 0; i <= n - j; i++)
        {
            if (data[i] > data[i + 1])
            {
                t = data[i];
                data[i] = data[i + 1];
                data[i + 1] = t;
            }
        }
    }
}

void swap_float_4(float *tnf4)
{
    int *tni4 = (int *)tnf4;
    *tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
             ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

float efloat(float value)
{
    union
    {
        char bytes[4];
        float flt;
    } ieee;

    // swap bytes
    unsigned char temp;
    ieee.flt = value;

    temp = ieee.bytes[3];
    ieee.bytes[3] = ieee.bytes[0];
    ieee.bytes[0] = temp;

    temp = ieee.bytes[2];
    ieee.bytes[2] = ieee.bytes[1];
    ieee.bytes[1] = temp;

    return (ieee.flt);
}

void ibm_to_float(int from[], int to[], int n, int endian)
/***********************************************************************
ibm_to_float - convert between 32 bit IBM and IEEE floating numbers
************************************************************************
Input::
from		input vector
to		output vector, can be same as input vector
endian		byte order =0 little endian (DEC, PC's)
			    =1 other systems
*************************************************************************
Notes:
Up to 3 bits lost on IEEE -> IBM

Assumes sizeof(int) == 4

IBM -> IEEE may overflow or underflow, taken care of by
substituting large number or zero

Only integer shifting and masking are used.
*************************************************************************
Credits: CWP: Brian Sumner,  c.1985
*************************************************************************/
{
    register int fconv, fmant, i, t;

    for (i = 0; i < n; ++i)
    {
        fconv = from[i];
        /* if little endian, i.e. endian=0 do this */
        if (endian == 0)
            fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
                    ((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);
        if (fconv)
        {
            fmant = 0x00ffffff & fconv;
            /* The next two lines were added by Toralf Foerster */
            /* to trap non-IBM format data i.e. conv=0 data  */
            if (fmant == 0)
                printf("mantissa is zero data may not be in IBM FLOAT Format !");
            t = (int)((0x7f000000 & fconv) >> 22) - 130;
            while (!(fmant & 0x00800000))
            {
                --t;
                fmant <<= 1;
            }
            if (t > 254)
                fconv = (0x80000000 & fconv) | 0x7f7fffff;
            else if (t <= 0)
                fconv = 0;
            else
                fconv = (0x80000000 & fconv) | (t << 23) | (0x007fffff & fmant);
        }
        to[i] = fconv;
    }
    return;
}

void float_to_ibm(int from[], int to[], int n, int endian)
/**********************************************************************
 float_to_ibm - convert between 32 bit IBM and IEEE floating numbers
***********************************************************************
Input:
from	   input vector
n	   number of floats in vectors
endian	   =0 for little endian machine, =1 for big endian machines

Output:
to	   output vector, can be same as input vector

***********************************************************************
Notes:
Up to 3 bits lost on IEEE -> IBM

IBM -> IEEE may overflow or underflow, taken care of by
substituting large number or zero

Only integer shifting and masking are used.
***********************************************************************
Credits:     CWP: Brian Sumner
***********************************************************************/
{
    register int fconv, fmant, i, t;

    for (i = 0; i < n; ++i)
    {
        fconv = from[i];
        if (fconv)
        {
            fmant = (0x007fffff & fconv) | 0x00800000;
            t = (int)((0x7f800000 & fconv) >> 23) - 126;
            while (t & 0x3)
            {
                ++t;
                fmant >>= 1;
            }
            fconv = (0x80000000 & fconv) | (((t >> 2) + 64) << 24) | fmant;
        }
        if (endian == 0)
            fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
                    ((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);

        to[i] = fconv;
    }
    return;
}

#endif
