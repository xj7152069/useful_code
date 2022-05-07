#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <armadillo>
#include "./include/xjc.h"

using namespace arma;
using namespace std;

int main()
{
    float dpx(0.0000015);
    float pi(3.1415926);
    int i,j,k,usecopy(0),ncpu(5);
    int nz(5000),nx(900),ny(1),nf(5000);
    int nt(5000),nbs(8),nbsx(1),sx(450);
    float dx(5.0),dt(0.0003),f,df,dw,dkz,dkx,kx,kz,w,v(1500),ro(1000);
    fmat dataline(nz,nx,fill::zeros);
    fmat datal(nt,1);

    df=1.0/dt/nt/nbs;
    dkx=2.0*pi/dx/nx/nbsx;

    ofstream outf1,outf2;
    ifstream inf1,inf2;

    fmat rawdatap(nt,nx),rawdatavz(nt,nx);
    inf1.open("./data/datafrom701/obndata.vz.cr.bin500.0");
    rawdatavz=dataread(nt,nx,inf1);
    inf1.close();
    inf1.open("./data/datafrom701/obndata.p.cr.bin500.0");
    rawdatap=dataread(nt,nx,inf1);
    inf1.close();

if(0)
{
    struct linerradon3d par;
    beamforming_parset(nx,ny,nz,par);
    
    par.dx=5,par.dy=5,par.dz=0.0003;
    par.npx=601;
    par.npy=1;
    par.dpx=0.000003;
    par.dpy=0.000003;
    par.px_center=0.0;
    par.py_center=0.0;   
    par.x_center=dx*(nx/2-sx);
    par.y_center=0;
    par.numthread=5;
    par.dig_n=1;
    par.nf1=1;
    par.nf2=1000;
    par.rulef1=1;
    par.rulef2=1000;
    beamforming_parupdate(par);
    cout<<"par.nf = "<<par.nf<<endl;
    dpx=par.dpx;
    fmat tpdata(nz,par.npx),tpdata_t(par.npx,nz);
/////////////////////////////////////////////////
    k=0;
    for(i=0;i<ny;i++){
        for(j=0;j<nx;j++){
            dataline.col(j)=rawdatap.col(k);
            k++;
        }
        dataline=get_blackman_rightwin2d\
            (dataline,100);
        dataline=get_blackman_leftwin2d\
            (dataline,100);
        dataline=get_blackman_downwin2d\
            (dataline,1000);

        par.data.col(i)=dataline.t();
        dataline.fill(0.0);
    }
////////////////////////////////////////////////
    beamforming_cleardata(par);
    linerradon(par); 
    par.dig_n1=0*nx*ny;
    par.dig_n2=nx*ny/nx/ny;
    beamforming_cleardata(par);
    beamforminginv3d(par);

    outf1.open("./result/model3/data3d.p.tp.bin");
    for(k=0;k<par.npy;k++)
    {
        tpdata_t=par.realdataTP.col(k);
        datawrite(tpdata=tpdata_t.t(),nz,par.npx,outf1);
    }
    outf1.close();

/////////////////////////////////////////////////
    k=0;
    for(i=0;i<ny;i++){
        for(j=0;j<nx;j++){
            dataline.col(j)=rawdatavz.col(k);
            k++;
        }
        dataline=get_blackman_rightwin2d\
            (dataline,100);
        dataline=get_blackman_leftwin2d\
            (dataline,100);
        dataline=get_blackman_downwin2d\
            (dataline,100);

        par.data.col(i)=dataline.t();
        dataline.fill(0.0);
    }
////////////////////////////////////////////////
    beamforming_cleardata(par);
    linerradon(par); 
    par.dig_n1=0*nx*ny;
    par.dig_n2=nx*ny/nx/ny;
    //beamforming_cleardata(par);
    //beamforminginv3d(par);

    outf1.open("./result/model3/data3d.vz.tp.bin");
    for(k=0;k<par.npy;k++)
    {
        tpdata_t=par.realdataTP.col(k);
        datawrite(tpdata=tpdata_t.t(),nz,par.npx,outf1);
    }
    outf1.close();

    cout<<"step-1: Radon is ok."<<endl;
    nx=par.npx;
}
///////////////////////////////////////////////////////////////////////
    //par.df=1.0/par.dz/par.nf;
if(0)
{
    demultiple2d par;

    par.dx=dt,par.nx=nt,par.dy=dx,par.ny=nx;
    par.agcwx=100,par.agcwy=10;
    par.datapagc.zeros(nt,nx);
    par.datavzagc.zeros(nt,nx);
    par.dataup.zeros(nt,nx);
    par.datadown.zeros(nt,nx);
    par.datap.zeros(nt,nx);
    par.datavz.zeros(nt,nx);
    par.data2d.zeros(nt,nx);
    par.datadict.zeros(nt,nx);
    par.datad.zeros(nt,nx);
    par.datah.zeros(nt,nx);
    par.datahd.zeros(nt,nx);
    par.nwp=(20),par.nwph=(10),par.nwpd=(10),par.nwphd=(10),\
    par.n1w=(4800),par.d1w=(4750),par.lmd=0.001;

    //fmat datavz(nt,nx);
    fmat datap(nt,nx),datavz(nt,nx);
    //dataread(par.datap,nt,nx,"./result/data3d.code.tp.bin");
    //par.datap*=(-1);
    //dataread(par.datavz,nt,nx,"./result/obndata.p.sum.tp.bin");
    
    dataread(par.datavz,nt,nx,"./result/model3/data3d.vz.tp.bin");
    dataread(par.datap,nt,nx,"./result/model3/data3d.p.tp.bin");
    par.datap=par.datap/(par.datap.max()-par.datap.min());
    par.datavz=par.datavz/(par.datavz.max()-par.datavz.min());
    //par.datavz=par.datavz*ro*v;
    par.datavz=get_blackman_upwin2d\
        (par.datavz,100);
    par.datap=get_blackman_upwin2d\
        (par.datap,100);
    //par.datap=par.datap-datap;
    //par.datavz=par.datavz-datavz;
    /*
    fmat t(nx,1);
    for(i=0;i<nx;i++)
    {
        t(i,0)=round(300+sqrt(25*(i-250)*(i-250)+200*200)/1500/0.0003);
        for(j=0;j<t(i,0);j++)
        {
            par.datap(j,i)=0;
            par.datavz(j,i)=0;
            datap(j,i)=0;
            datavz(j,i)=0;
        }
    }*/

    fmat datap2(nt,nx);
    par.datap=par.datap/v/ro;
    
/////////////////////////////////////////////////////////////

    datap=par.datap;
    //getagc_demultiple2d(par,par.datap,par.datapagc);
    //datap=par.datap;
    //deagc_demultiple2d(par,datap,par.datapagc);
    datavz=par.datavz;
    //getagc_demultiple2d(par,par.datavz,par.datavzagc);
    //datavz=par.datavz;
    //deagc_demultiple2d(par,datavz,par.datavzagc);
    
    float a0(1.0/1500),q0;
    int med(nx/2);
    datap.col(med)=datap.col(med)*a0;
    for(i=1;i<(med-1);i++){
        q0=a0*a0-dpx*dpx*i*i;
        if(q0>=0)
            q0=sqrt(q0);
        else
            q0=0;
        datap.col(med-i)=datap.col(med-i)*q0;
        datap.col(med+i)=datap.col(med+i)*q0;
    }

    for(k=0;k<nx;k++){
        cout<<k<<"|"<<endl;
        par.datah.col(k)=hilbert1D(datal=datap.col(k),nt,dt);
        par.datad.col(k)=diff1D(datal=datap.col(k),nt,dt);
        par.datahd.col(k)=diff1D(datal=par.datah.col(k),nt,dt);
    }
    datawrite(par.datah,nt,nx,"./media/obndata.p2.h.bin");
    datawrite(par.datad,nt,nx,"./media/obndata.p2.d.bin");
    datawrite(par.datahd,nt,nx,"./media/obndata.p2.hd.bin");

    dataread(par.datah,nt,nx,"./media/obndata.p2.h.bin");
    dataread(par.datad,nt,nx,"./media/obndata.p2.d.bin");
    dataread(par.datahd,nt,nx,"./media/obndata.p2.hd.bin");

////////////////////////////////////////////////////////////////////
    par.data2d=par.datavz, par.datadict=par.datap;
    single_trace_dewave(par,ncpu);
    cout<<par.data2d.n_rows<<"|"<<par.data2d.n_cols<<endl;
    //deagc_demultiple2d(par,par.data2d,par.datavzagc);
    //datawrite(datap3=datap3/1,nt,nx,"./result/obndata.p.sum.tp.bin");
    datawrite(par.dataup,par.n1,par.n2,"./result/model3/obndata.p.upgoing.tp.bin");
    datawrite(par.datadown,par.n1,par.n2,"./result/model3/obndata.p.downgoing.tp.bin");
}
/////////////////////////////////////////////////////////////
if(0)
{
    fmat datapr(nt,nx),dataup(nt,nx),datadown(nt,nx);
    cx_fmat dataupfft(nt,nx),datadownfft(nt,nx),\
        datadownfft_t(nx,nt),dataprfft(nt,nx,fill::zeros);
    datapr(0,0)=0.001;

    dataread(dataup,nt,nx,"./result/model3/obndata.p.upgoing.tp.bin");
    dataread(datadown,nt,nx,"./result/model3/obndata.p.downgoing.tp.bin");

    for(k=0;k<nx;k++){
        dataupfft.col(k)=fft(dataup.col(k));
        datadownfft.col(k)=fft(datadown.col(k));   
    }
    datadownfft_t=datadownfft.t();
    for(k=0;k<nx;k++){
        for(i=0;i<nt/2;i++){
            dataprfft(i,k)=(dataupfft(i,k)*datadownfft_t(k,i))/\
                (datadownfft(i,k)*datadownfft_t(k,i)+datapr(0,0));
        }
    }
    for(k=0;k<nx;k++){
        datapr.col(k)=real(ifft(dataprfft.col(k)));
    }

    datawrite(datapr,nt,nx,"./result/model3/obndata.p.pr.tp.bin");
}
/////////////////////////////////////////////////////////////
if(0)
{
    struct linerradon3d par;
    nx=900;
    beamforming_parset(nx,ny,nz,par);
    
    par.dx=5,par.dy=5,par.dz=0.0003;
    par.npx=601;
    par.npy=1;
    par.dpx=0.000003;
    par.dpy=0.000003;
    par.px_center=0.0;
    par.py_center=0.0;  
    par.x_center=dx*(nx/2-sx);
    par.y_center=0; 
    par.numthread=5;
    par.dig_n=1;
    par.nf1=1;
    par.nf2=nz/2-1;
    par.rulef1=1;
    par.rulef2=nz/2-1;
    beamforming_parupdate(par);
    cout<<"par.nf = "<<par.nf<<endl;
    dpx=par.dpx;
/////////////////////////////////////////////////

    fmat datatp(nz,par.npx);
    datatp=dataread(nz,par.npx,"./result/model3/obndata.p.pr.tp.bin");
    float a0(1.0/1500),q0;
    int med(par.npx/2);
    datatp.col(med)=datatp.col(med)*a0;
    for(i=1;i<(med-1);i++){
        q0=a0*a0-dpx*dpx*i*i;
        if(q0>=0)
            q0=1;
        else
            q0=0;
        datatp.col(med-i)=datatp.col(med-i)*q0;
        datatp.col(med+i)=datatp.col(med+i)*q0;
    }
    for(i=0;i<5;i++){
        datatp.col(i).fill(0);
        datatp.col(par.npx-i-1).fill(0);
    }

    datawrite(datatp,nz,par.npx,"./result/model3/obndata.pr.tp.process.bin");
    par.realdataTP.col(0)=datatp.st();
////////////////////////////////////////////////
    rebuildsignal(par);

    fmat filterdata_t(nx,nz), db(1,1,fill::zeros);
    outf1.open("./result/model3/obndata.pr.rebuild.bin");
    for(k=0;k<ny;k++)
    {
        filterdata_t=par.realrebuildtx.col(k);
        datawrite(dataline=filterdata_t.t(),par.nz,par.nx,outf1);
        //filterdata_t=par.realrebuildtx.col(k);
    }
    outf1.close();
///////////////////////////////////////////////////////
    datatp=dataread(nz,par.npx,"./result/model3/obndata.p.upgoing.tp.bin");
    par.realdataTP.col(0)=datatp.st();
    rebuildsignal(par);
    outf1.open("./result/model3/obndata.upgoing.rebuild.bin");
    for(k=0;k<ny;k++)
    {
        filterdata_t=par.realrebuildtx.col(k);
        datawrite(dataline=filterdata_t.t(),par.nz,par.nx,outf1);
        //filterdata_t=par.realrebuildtx.col(k);
    }
    outf1.close();
////////////////////////////////////////////////////////////
    datatp=dataread(nz,par.npx,"./result/model3/obndata.p.downgoing.tp.bin");
    par.realdataTP.col(0)=datatp.st();
    rebuildsignal(par);
    outf1.open("./result/model3/obndata.downgoing.rebuild.bin");
    for(k=0;k<ny;k++)
    {
        filterdata_t=par.realrebuildtx.col(k);
        datawrite(dataline=filterdata_t.t(),par.nz,par.nx,outf1);
        //filterdata_t=par.realrebuildtx.col(k);
    }
    outf1.close();

}

/////////////////////////////////////////////////////////////
if(0)
{
    nx=900,nz=nt;
    int Z(nz),X(nx),T(nz);
    float dx(5),dz(5),dt(0.0003),f0(30),deepth(40*dz);
    int sx(500),sz(30+20/(1+exp(0.015*(sx-i)))),vocean(1500);
////////////////////////////////////////////////
float *obnz;
obnz=new float[nx];
    for(i=0;i<nx;i++)
    {
        float pi(3.1415926);
        float phase=(float(i)/nx)*pi*2.0;
        //obnz[i]=round(60+9*sin(phase));
        obnz[i]=30+20/(1+exp(0.015*(450-i)));
    }
////////////////////////////////////////////////
    int nt2(nz*2);
    nt=nz;
    fmat data1(nt2,nx,fill::zeros),data2(nt2,nx),data(nt,nx);
    fmat suf(T,nx),sufp(T,nx),sufs(T,nx),ureal(Z,X),tu2u1(nx,nx);
    float df(1.0/dt/nt2),d;

    //data=dataread(nt,nx,"./result/useful/data3d.demulti2.rebuild.bin");
    data=rawdatap;

    get_blackman_leftwin2d(data, 100);
    get_blackman_rightwin2d(data, 100);
    
    for(i=0;i<nx;i++){
        for(j=0;j<nx;j++){
            deepth=j;
            deepth=dz*(20+20/(1+exp(0.015*(450-deepth))));
            //deepth=dz*35;
            d=abs(i-j)*dx;
            d=sqrt(d*d+deepth*deepth);
            tu2u1(i,j)=-d/vocean-0.05;
        }
    }
    for(i=0;i<nt;i++){
        data1.row(i)=data.row(i);
    }
    multiple_code2d(data2,data1,tu2u1,df,0,4000,ncpu);

    datawrite(data2,T,nx,"./result/model3/obndata.cr.p.code.bin");
    //datawrite(data2,T,nx,"./result/obndata.demulti2.code.bin");
}
//////////////////////////////////////////////////////////////
if(1)
{
    demultiple2d par;

    par.dx=dt,par.nx=nt,par.dy=dx,par.ny=nx;
    par.agcwx=100,par.agcwy=10;
    par.datapagc.zeros(nt,nx);
    par.datavzagc.zeros(nt,nx);
    par.dataup.zeros(nt,nx);
    par.datadown.zeros(nt,nx);
    par.datap.zeros(nt,nx);
    par.datavz.zeros(nt,nx);
    par.data2d.zeros(nt,nx);
    par.datadict.zeros(nt,nx);
    par.datad.zeros(nt,nx);
    par.datah.zeros(nt,nx);
    par.datahd.zeros(nt,nx);
    par.nwp=(20),par.nwph=(10),par.nwpd=(10),par.nwphd=(10),\
    par.n1w=(4800),par.d1w=(4750),par.lmd=0.001;

    //fmat datavz(nt,nx);
    fmat datap(nt,nx),datavz(nt,nx);

    
    dataread(par.datavz,nt,nx,"./result/model3/obndata.downgoing.rebuild.bin");
    dataread(par.datap,nt,nx,"./result/model3/data3d.demulti2.rebuild.bin");
    par.datap=-par.datap/(par.datap.max()-par.datap.min());
    par.datavz=par.datavz/(par.datavz.max()-par.datavz.min());
    
/////////////////////////////////////////////////////////////

    datap=par.datap;
    datavz=par.datavz;

    for(k=0;k<nx;k++){
        cout<<k<<"|"<<endl;
        par.datah.col(k)=hilbert1D(datal=datap.col(k),nt,dt);
        par.datad.col(k)=diff1D(datal=datap.col(k),nt,dt);
        par.datahd.col(k)=diff1D(datal=par.datah.col(k),nt,dt);
    }
    datawrite(par.datah,nt,nx,"./media/obndata.p2.h.bin");
    datawrite(par.datad,nt,nx,"./media/obndata.p2.d.bin");
    datawrite(par.datahd,nt,nx,"./media/obndata.p2.hd.bin");

    dataread(par.datah,nt,nx,"./media/obndata.p2.h.bin");
    dataread(par.datad,nt,nx,"./media/obndata.p2.d.bin");
    dataread(par.datahd,nt,nx,"./media/obndata.p2.hd.bin");

/***************************/
    par.data2d=par.datavz, par.datadict=par.datap;
    single_trace_dewave(par,ncpu);
    cout<<par.data2d.n_rows<<"|"<<par.data2d.n_cols<<endl;

    datawrite(par.dataup,par.n1,par.n2,"./result/model3/obndata.downgoing.process.bin");
}
////////////////////////////////////////////////////////////

    cx_fmat test(1,1);
    test(0,0).real(1);
    test(0,0).imag(2);
    test=test/2;
    test.print();

    return 0;
}




    /*
    //dataread(datavz,nt,nx,"shotvz.bin");
    //dataread(datap,nt,nx,"shotp.bin");
    for(i=0;i<nt;i++){
        for(j=0;j<nx;j++){
            if(isnan(abs(datavz(i,j))))
            datavz(i,j)=0;
            if(isnan(abs(datap(i,j))))
            datap(i,j)=0;
        }
    }
    datawrite(datavz,nt,nx,"shotvz.bin");
    datawrite(datap,nt,nx,"shotp.bin");*/




