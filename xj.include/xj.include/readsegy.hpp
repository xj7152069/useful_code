#ifndef READSEGY_HPP
#define READSEGY_HPP

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "../xjc.h"

using namespace std;

void segyhead_endianget(segyhead & head);
void segyhead_initialize(segyhead & head);
template <typename T1> inline void endianchange(T1 & a);
template <typename T1> inline float getendianchange(T1 a);
inline float floatendianchange(float a);
void segyhead_open(segyhead & head);
void fmatendianchange(fmat & data,int n1,int n2);
void segyhead_readonetrace_tofmat(segyhead & head, fmat & trace);
void segyhead_readoneline_tofmat(segyhead & head, fmat & data);

///////////////////////////////////////////////////////////////////
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

void segyhead_readonetrace_tofmat(segyhead & head, fmat & trace)
{
    int nz;
    head.infile.read((char *)(&head.head2), sizeof(head.head2));
    nz=int(getendianchange(head.head2.ns));
    trace.zeros(nz,1);
    trace=dataread(nz,1,head.infile);
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