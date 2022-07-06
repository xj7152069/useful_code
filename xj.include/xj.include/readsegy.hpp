#ifndef READSEGY_HPP
#define READSEGY_HPP

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "../xjc.h"

using namespace std;

float ibm2ieee (float fromf);
void ibm_to_float(unsigned int from[], unsigned int to[], int n, int endian);
void segyhead_endianget(segyhead & head);
void segyhead_initialize(segyhead & head);
template <typename T1> inline void endianchange(T1 & a);
template <typename T1> inline float getendianchange(T1 a);
inline float floatendianchange(float a);
void segyhead_open(segyhead & head);
void fmatendianchange(fmat & data,int n1,int n2);
fmat segyhead_readonetrace_tofmat(segyhead & head, fmat & trace);
//void segyhead_readoneline_tofmat(segyhead & head, fmat & data);

///////////////////////////////////////////////////////////////////
fmat segyhead_readonetrace_tofmat(segyhead & head, fmat & trace)
{
    int nz;
    head.infile.read((char *)(&head.head2), sizeof(head.head2));
    if(head.endian=='b')
        nz=int(getendianchange(head.head1.hns));
    else
        nz=int(head.head1.hns);
    trace.zeros(nz,1);
    fmat rawtrace;
    rawtrace.copy_size(trace);

    if(head.isibm && head.endian=='b'){
        int i,j;
        float *float_p;
        unsigned int *int_p;
        float IBMFloatBytes;
        float_p=(float *)(&IBMFloatBytes);
        int_p=(unsigned int *)(&IBMFloatBytes);
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                rawtrace(i,j)=*float_p;
                IBMFloatBytes=getendianchange(IBMFloatBytes);
                ibm_to_float(int_p, int_p, 1, 1);
                trace(i,j)=*float_p;
                if(isinf(abs(trace(i,j)))){
                    trace(i,j)=0;
                }
            }
        }
            //cout<<i<<"ok"<<IBMFloatBytes<<endl;

    }
    else if(head.isibm && head.endian=='l'){
        int i,j;
        float *float_p;
        unsigned int *int_p;
        float IBMFloatBytes;
        float_p=(float *)(&IBMFloatBytes);
        int_p=(unsigned int *)(&IBMFloatBytes);
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                rawtrace(i,j)=*float_p;
                ibm_to_float(int_p, int_p, 1, 1);
                trace(i,j)=*float_p;
                if(isinf(abs(trace(i,j)))){
                    trace(i,j)=0;
                }
            }
        }
    }
    else if(!head.isibm && head.endian=='b'){
        int i,j;
        float *float_p;
        unsigned int IBMFloatBytes;
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                float_p=(float *)(&IBMFloatBytes);
                rawtrace(i,j)=*float_p;
                IBMFloatBytes=getendianchange(IBMFloatBytes);
                float_p=(float *)(&IBMFloatBytes);
                trace(i,j)=*float_p;
            }
        }
    }
    else if(!head.isibm && head.endian=='l'){
        int i,j;
        float *float_p;
        unsigned int IBMFloatBytes;
        for(i=0;i<trace.n_rows;i++){
            for(j=0;j<trace.n_cols;j++){
                head.infile.read((char *)&IBMFloatBytes, 4); 
                float_p=(float *)(&IBMFloatBytes);
                rawtrace(i,j)=*float_p;
                float_p=(float *)(&IBMFloatBytes);
                trace(i,j)=*float_p;
            }
        }
    }

    return rawtrace;
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
        head.endian = 'l'; //小端，高位字节在前，低位字节在后
        cout<<"the format of the SEGY-data is little-Endian: "<<\
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
            head.endian = 'b'; //大端(654321)，与小端相反，需要调整为小端
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
float ibm2ieee (float fromf) 
{
    float to;
    int *from = (int *)&fromf;
    // int fconv, fmant, i, t;
    int fconv, fmant, t;
    fconv = from[0];

    if (fconv)
    {
        fmant = 0x00ffffff & fconv;
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
    to = fconv;
    return to;
}
/*
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

   ibmCode= (ibmCode<<1); // 左移移出符号位, ibmCode= EEEEEEEM MMMMMMMM MMMMMMMM MMMMMMM0

    unsigned int float_IBM(num);

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
*/
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

void ibm_to_float(unsigned int from[], unsigned int to[], int n, int endian)
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
    unsigned int fconv, fmant, i, t;

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
            t = (unsigned int)((0x7f000000 & fconv) >> 22) - 130;
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
