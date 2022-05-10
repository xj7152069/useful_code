
#define ARMA_DONT_USE_BLAS //!!!!
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "./include/xjc.h"

using namespace std;
using namespace arma;
void radoninv(struct linerradon3d & par, bool recoverdata);

/*core function:
    beamforming_parset(nx,ny,nz,par);
    beamforming_parupdate(par);
    linerradon(par); 
    beamforming_cleardata(par);
    beamforminginv3d(par);
    rebuildsignal(par);
*/
int main()
{
    int i,j,k,usecopy(0),i1,j1,k1;
    //input 3D data
    string infile;
    int nz2(3501),nz(3501),nx(70),ny(1),nf(0);
    float dt;

    //cout<<"please input the file:";
    //cin>>infile;
    infile="./data/qj.data.1801.segy";
    char file[99],file2[99];
    file[0]='\0';
    for(k=0;k<infile.length();k++){
        file[k]=infile[k];
    }
    file[k]='\0';
    //strcat(file,"./data/qj.data.1801.segy");

    ofstream outf1,outf2,outf3,outf4;
    ifstream inf1,inf2,inf3;

    file2[0]='\0';
    strcat(file2,file);
    strcat(file2,"rebuild.bin");
    outf1.open(file2);
    
    file2[0]='\0';
    strcat(file2,file);
    strcat(file2,"err.bin");
    //outf3.open(file2);

    struct linerradon3d par;

    segyhead head;
    head.filename[0]='\0';
    strcat(head.filename,file);
    segyhead_open(head);
    segyhead_endianget(head);
    head.begnx=0;

    inf3.open("./data/rebuild.1801.bin");
    float offset,offset1,offset2,doffset(1000),js(0),lsjs(0),minnx(30);
    fmat alltp(nz2,161,fill::zeros);//???????
    alltp.fill(0);

while(!head.infile.eof()){
    nx=0;js=0;
    outf2.open("./data1.bin");
    outf4.open("./data1.offset.bin");

    //read segy data
    fmat a(1,1);
    while(!head.infile.eof())
    {
        float test;
        if(head.begnx==0){
            head.begnx=1;
            segyhead_readonetrace_tofmat(head,head.data);
            nz=head.data.n_rows;
            if(head.endian=='l'){
                offset1=getendianchange(head.head2.offset);
                offset=offset1;
            }else{
                offset1=head.head2.offset;
                offset=offset1;
            }
            nz=head.data.n_rows;
            //head.data=dataread(nz,1,inf3);

            datawrite(head.data,nz,1,outf2);
            a(0,0)=offset;
            datawrite(a,1,1,outf4);
            nx++;
            cout<<offset1<<"||"<<endl;
        }
        else if(head.begnx==1){
            segyhead_readonetrace_tofmat(head,head.data);
            nz=head.data.n_rows;
            //head.data=dataread(nz,1,inf3);

            if(head.endian=='l'){
                offset=getendianchange(head.head2.offset);
            }else{
                offset=(head.head2.offset);
            }
            js+=abs(offset-offset1);
            offset1=offset;
            if(js<doffset || nx<minnx){
            //if(js<doffset || nx<minnx || nx>=minnx){
                nz=head.data.n_rows;
                datawrite(head.data,nz,1,outf2);  
                a(0,0)=offset;
                datawrite(a,1,1,outf4);
                nx++;
            }
            else if(js>=doffset && nx>=minnx){
                head.begnx=2;
                break; //!!!!!!!!!!!
            }
        }
        else if(head.begnx==2){
            nz=head.data.n_rows;
            datawrite(head.data,nz,1,outf2);    
            a(0,0)=offset;
            datawrite(a,1,1,outf4);    
            head.begnx=1;
            nx++;
            cout<<offset1<<"||"<<endl;
        }

    head.infile.read((char *)(&test), sizeof(test));
    if(head.infile.eof())
    {break;}
    else
    {
        head.infile.seekg(-1*sizeof(test),ios::cur);
    }
    }
    outf2.close();
    outf4.close();
    if(head.endian=='l'){
        dt=getendianchange(head.head2.dt);
    }else{
        dt=(head.head2.dt);
    }
    dt/=1000000;

/////////////////////radon inv/////////////////////
//input data: (row,col,slice) of par.data is (nx,ny,nz
//output recover data: (row,col,slice) of par.realrebuildtx is (nx,ny,nz)
//output tau-p: (row,col,slice) of par.realdataTP is (npx,npy,nz)
//nx: Spatial sampling nummber
//ny: trace nummber
//nz: time sampling nummber
//npx: X Ray parameters sampling nummber
//npy: Y Ray parameters sampling nummber
    nz2=nz;
    if(1)
    {
        int nh(nx*ny);
        fmat basedata(nz2,nh),site(1,nh);
        inf1.open("data1.bin");  //local data
        dataread(basedata,nz2,nh,inf1);
        inf1.close();
        inf1.open("data1.offset.bin");  //offset
        dataread(site,1,nh,inf1);  //Seismic trace coordinates
        inf1.close();

//////////////////////////radon par-set////////////////////////////
        beamforming_parset(nx,ny,nz,par);
        par.dpx=0.000005;
        par.dpy=0.000005;
        par.dz=dt;

//The default px of central channel is zero
        par.npx=161;
        par.px_center=par.dpx*80;
        par.npy=1;
        par.py_center=0.0;   

        //ncpu
        par.numthread=5;

//Frequency calculation range (number)
        par.nf1=1;
        par.nf2=1000; 
//Low frequency constraint range (number)
        par.rulef1=1;
        par.rulef2=500;
//regularization parameter
        par.dig_n1=0*nx*ny;  //L1
        par.kerpar1=0.0;     //kernel function par
        par.dig_n2=nx*ny*1;  //L2, Tikhonov 
//Parameters updated
        beamforming_parupdate(par);
//Seismic trace coordinates

        for(k=0;k<par.nx;k++)
        {
            par.x_coord(k,0)=site(0,k)/2;
        }
//////////////////////////////////////////////////////////////////
    //data input, (row,col,slice) of par.data is (nx,ny,nz)
    //nx: Spatial sampling nummber
    //ny: trace nummber
    //nz: time sampling nummber
        fmat dataline(nz,nx,fill::zeros);
        dataline=basedata;

/******************************************/
        par.data.col(0)=dataline.t();
/******************************************/

//////////////////////////////////////////////////

        radoninv(par,true);

///////////////////output tau-p///////////////////
        fmat tpdata(nz,par.npx),tpdata_t(par.npx,nz);
    //output tau-p, (row,col,slice) of par.realdataTP is (npx,npy,nz)
    //npx: X Ray parameters sampling nummber
    //npy: Y Ray parameters sampling nummber
    //nz: time sampling nummber
        for(k=0;k<par.npy;k++){
/******************************************/
            tpdata_t=par.realdataTP.col(k);
/******************************************/
            tpdata=tpdata_t.t();
            for(i=0;i<nz2;i++){
                alltp.row(i)+=tpdata.row(i);
            }
        }
////////////////////output recover data/////////////////////
    //output recover data, (row,col,slice) of par.realrebuildtx is (nx,ny,nz)
    //npx: X Ray parameters sampling nummber
    //npy: Y Ray parameters sampling nummber
    //nz: time sampling nummber
        fmat filterdata_t(nx,nz);
        for(k=0;k<ny;k++)
        {
/******************************************/
            filterdata_t=par.realrebuildtx.col(k);
/******************************************/
            datawrite(dataline=filterdata_t.t(),nz2,par.nx,outf1);
        }
        cout<<"step - 2 is ok: "<<endl;
    }

}
    outf1.close();
    //outf3.close();
    inf3.close();
    datawrite(alltp,nz2,par.npx,"alltp.bin");

    if(1)
    {
        nx=800;
        int nh(nx*ny);
        fmat basedata(nz2,nh),site(1,nh);
        inf1.open("data1.bin");  //local data
        dataread(basedata,nz2,nh,inf1);
        inf1.close();
        inf1.open("data1.offset.bin");  //offset
        dataread(site,1,nh,inf1);  //Seismic trace coordinates
        inf1.close();

//////////////////////////radon par-set////////////////////////////
        beamforming_parset(nx,ny,nz,par);
        par.dpx=0.000005;
        par.dpy=0.000005;
        par.dz=dt;

//The default px of central channel is zero
        par.npx=161;
        par.px_center=par.dpx*80;
        par.npy=1;
        par.py_center=0.0;   

        //ncpu
        par.numthread=5;

//Frequency calculation range (number)
        par.nf1=1;
        par.nf2=1000; 
//Low frequency constraint range (number)
        par.rulef1=1;
        par.rulef2=500;
//regularization parameter
        par.dig_n1=0*nx*ny;  //L1
        par.kerpar1=0.0;     //kernel function par
        par.dig_n2=nx*ny*1;  //L2, Tikhonov 
//Parameters updated
        beamforming_parupdate(par);
//Seismic trace coordinates

        for(k=0;k<par.nx;k++)
        {
            par.x_coord(k,0)=5*k;
        }
//////////////////////////////////////////////////////////////////
    //data input, (row,col,slice) of par.data is (nx,ny,nz)
    //nx: Spatial sampling nummber
    //ny: trace nummber
    //nz: time sampling nummber
        fmat dataline(nz,nx,fill::zeros);
        dataline=basedata;

/******************************************/
        //par.data.col(0)=dataline.t();
/******************************************/

//////////////////////////////////////////////////


///////////////////output tau-p///////////////////
        fmat tpdata(nz,par.npx),tpdata_t(par.npx,nz);
    //output tau-p, (row,col,slice) of par.realdataTP is (npx,npy,nz)
    //npx: X Ray parameters sampling nummber
    //npy: Y Ray parameters sampling nummber
    //nz: time sampling nummber
        for(k=0;k<par.npy;k++){
/******************************************/
            tpdata=dataread(nz,par.npx,"alltp.bin");
            par.realdataTP.col(k)=tpdata.t();
        }
        rebuildsignal(par);
////////////////////output recover data/////////////////////
    //output recover data, (row,col,slice) of par.realrebuildtx is (nx,ny,nz)
    //npx: X Ray parameters sampling nummber
    //npy: Y Ray parameters sampling nummber
    //nz: time sampling nummber
        fmat filterdata_t(nx,nz);
        for(k=0;k<ny;k++)
        {
/******************************************/
            filterdata_t=par.realrebuildtx.col(k);
/******************************************/
            datawrite(dataline=filterdata_t.t(),nz2,par.nx,"allrebuild.bin");
        }
        cout<<"step - 2 is ok: "<<endl;
    }

    return 0;
}

void radoninv(struct linerradon3d & par, bool recoverdata=false)
{
    struct linerradon3d * ppar;
    ppar=&par;
    //slant stack firstly
    linerradon(ppar[0]); 

    //inv radon
    beamforming_cleardata(ppar[0]);
    beamforminginv3d(ppar[0],0);
    //beamformingCG3d(par); //not stable
    if(recoverdata){
        rebuildsignal(ppar[0]);
    }

}
