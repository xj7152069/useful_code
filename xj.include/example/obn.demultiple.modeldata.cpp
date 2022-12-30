#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "./include/xjc.h"
#include <omp.h>
#include <armadillo>
using namespace arma;
using namespace std;

int main(int argc, char * argv[])
{
    int i,j,k,usecopy(0),i1,j1,k1,nzdata(12000);
    float fmax(250),frule(50),factor_L2(0.1), factor_L1(0),\
        iterations_num(45), residual_ratio(0.01),ncpu(atoi(argv[1]));
    int nz(nzdata),nx(900),ny(1),nf(nzdata),nt(nz);
    float npx(501),npy(1),dt(0.0003),dx(5),dy(5),dz(5),dpx(0.000003),dpy(0.000003);
    float x0(-dx*50),y0(-dy*ny/2),px0(-dpx*npx/2),py0(-dpy*npy/2);
    float f,df,dw,dkz,dkx,kx,kz,w,velocity_water(1500),ro(1000);
    fcube data3d_p(nx,ny,nz),data3d_vz(nx,ny,nz),datatp3d(npx,npy,nz),\
        datarecover3d(nx,ny,nz),dataerr3d(nx,ny,nz);
    fmat data,datatp,datarecover,dataerr;
    fmat data_t,datatp_t,datarecover_t,dataerr_t;
    fmat coordx(nx,ny),coordy(nx,ny),seabase_depth(nx,ny);
    data.zeros(nz,nx);
    datatp.zeros(nz,nx);
    datarecover.zeros(nz,nx);
    dataerr.zeros(nz,nx);
    data_t.zeros(nx,nz);
    datatp_t.zeros(nx,nz);
    datarecover_t.zeros(nx,nz);
    dataerr_t.zeros(nx,nz);
    coordx.zeros(nx,ny);
    coordy.zeros(nx,ny);
    seabase_depth.zeros(nx,ny);
    data3d_p.zeros(nx,ny,nz);
    data3d_vz.zeros(nx,ny,nz),datatp3d.zeros(npx,npy,nz),\
    datarecover3d.zeros(nx,ny,nz),dataerr3d.zeros(nx,ny,nz);

    fmat datawin;
    datawin=get_datawin_tp3d(nt,npx,npx,dt,dpx,dpx,\
        px0,px0, 175,1500.0,20001,20001,0.1,0.1);
    datawin+=round(0.05/dt);
    fmat datawin2;
    datawin2=get_datawin_tp3d(nt,npx,npx,dt,dpx,dpx,\
        px0,px0, 670,1500.0,20001,20001,0.1,0.1);
    datawin2+=round(0.05/dt);

    ofstream outf1,outf2;
    ifstream inf1,inf2;
    char file1[99],file2[99];
    float sxbeg(300),sxend(300),dsx(100),sx;

for(sx=sxbeg;sx<=sxend;sx+=dsx){

    x0=(-dx*(sx-50));
    for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
        coordx(i,j)=dx*i+x0;
        coordy(i,j)=dy*j+y0;
        seabase_depth(i,j)=35*dz;
        //seabase_depth(i,j)=dz*(-10+40+20/(1+exp(0.1*(nx/2-i)))-0.04*i\
            +10/(1+exp(0.1*(3*nx/4-i)))+10/(1+exp(0.1*(nx/4-i))));
    }}
    file1[0]='\0';
    strcat(file1,"./data/data.cr.model4/obndata.vz.cr.bin");
    strcat(file1,numtostr(sx,5));
    inf1.open(file1);
    dataread3d_bycol_transpose(data3d_vz,nzdata,nx,inf1);
    inf1.close();
    data_t=data3d_vz.col(0);
    data_t=get_blackman_downwin2d(data_t,100);
    data_t=get_blackman_upwin2d(data_t,100);
    data_t=get_blackman_rightwin2d(data_t,2000);
    data3d_vz.col(0)=data_t;

    file1[0]='\0';
    strcat(file1,"./data/data.cr.model4/obndata.p.cr.bin");
    strcat(file1,numtostr(sx,5));
    inf1.open(file1);
    dataread3d_bycol_transpose(data3d_p,nzdata,nx,inf1);
    data3d_p=data3d_p/data3d_p.max()/0.01;
    inf1.close();
    data_t=data3d_p.col(0);
    data_t=get_blackman_downwin2d(data_t,100);
    data_t=get_blackman_upwin2d(data_t,100);
    data_t=get_blackman_rightwin2d(data_t,2000);
    data3d_p.col(0)=data_t;

if(0)
{
    k=Beamforming_CG_3D(datatp3d,datarecover3d,dataerr3d,\
        data3d_vz,coordx,coordy, nz, nx,ny,dt,\
        npx,px0, dpx,npy,py0, dpy,fmax,frule,\
        ncpu, factor_L2, factor_L1, iterations_num, residual_ratio);
    datawrite3d_bycol_transpose(datatp3d,nz,npx,"./result/model4/tp.vz.bin");

    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,fmax,frule,ncpu);

    k=Beamforming_CG_3D(datatp3d,datarecover3d,dataerr3d,\
        data3d_p,coordx,coordy, nz, nx,ny,dt,\
        npx,px0, dpx,npy,py0, dpy,fmax,frule,\
        ncpu, factor_L2, factor_L1, iterations_num, residual_ratio);
    datawrite3d_bycol_transpose(datatp3d,nz,npx,"./result/model4/tp.p.bin");

}
///////////////////////////////////////////////////////////////////////
//PZ sum
if(0)
{
    demultiple2d par;
    par.dt=dt,par.nt=nt,par.dx=dpx,par.nx=npx;
    par.data_remove_wave.zeros(par.nt,par.nx);
    par.data_antiremove_wave.copy_size(par.data_remove_wave);
    par.data2d.copy_size(par.data_remove_wave);
    par.datadict.copy_size(par.data_remove_wave);
    par.datad.copy_size(par.data_remove_wave);
    par.datah.copy_size(par.data_remove_wave);
    par.datahd.copy_size(par.data_remove_wave);
    par.nwp=(40),par.nwph=(20),par.nwpd=(20),par.nwphd=(20),\
    par.n1w=(800),par.d1w=(690),par.lmd=0.00001;
    par.datawinbeg.zeros(par.nx,1);
    par.datawinend.copy_size(par.datawinbeg);

    par.datawinend.col(0)=datawin.col(par.nx/2);
    //datawrite(par.datawinend,npx,1,"datawin.bin");

    fmat datap,datavz;
    datap.copy_size(par.data_remove_wave);
    datavz.copy_size(par.data_remove_wave);

    dataread(datavz,"./result/model4/tp.vz.bin");
    dataread(datap,"./result/model4/tp.p.bin");
    datap=datap/(datap.max()-datap.min());
    datavz=datavz/(datavz.max()-datavz.min());
    datavz=get_blackman_upwin2d\
        (datavz,100);
    datap=get_blackman_upwin2d\
        (datap,100);

/////////////////////////////////////////////////////////////

    float a0(1.0/velocity_water),q0;
    int med(par.nx/2);
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

omp_set_num_threads(ncpu);
#pragma omp parallel for
    for(k=0;k<par.nx;k++){
        fmat datal(nt,1);
        //cout<<k<<"|";
        par.datah.col(k)=hilbert1D(datal=datap.col(k),nt,dt);
        par.datad.col(k)=diff1D(datal=datap.col(k),nt,dt);
        par.datahd.col(k)=diff1D(datal=par.datah.col(k),nt,dt);
    }

    datawrite(par.datah,"./media/obndata.p2.h.bin");
    datawrite(par.datad,"./media/obndata.p2.d.bin");
    datawrite(par.datahd,"./media/obndata.p2.hd.bin");

    dataread(par.datah,"./media/obndata.p2.h.bin");
    dataread(par.datad,"./media/obndata.p2.d.bin");
    dataread(par.datahd,"./media/obndata.p2.hd.bin");

////////////////////////////////////////////////////////////////////
    par.data2d=datavz, par.datadict=datap;
    //single_trace_dewave(par,ncpu);
    par.datawinbeg.col(0)=datawin.col(npx/2)-200;
    par.datawinend.col(0)=datawin.col(npx/2)+200;
    single_trace_dewave_withdatawin(par,ncpu);
    cout<<par.data2d.n_rows<<"|"<<par.data2d.n_cols<<endl;
    datawrite(par.data_remove_wave,"./result/model4/tp.p.upgoing.bin");
    //datavz=datavz/(datavz.max()-datavz.min());
    //datawrite(par.data_remove_wave-=datavz,"./result/model4/tp.p.upgoing.err.bin");
    datatp3d.col(0)=par.data_remove_wave.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,fmax,frule,ncpu);
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.vz.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
////////////////////////////////////////////////////////////////////
    par.datawinbeg.col(0)=datawin2.col(npx/2)-150;
    par.datawinend.col(0)=datawin2.col(npx/2)+150;
    single_trace_dewave_withdatawin(par,ncpu);
    cout<<par.data2d.n_rows<<"|"<<par.data2d.n_cols<<endl;
    datawrite(par.data_remove_wave,nt,npx,"./result/model4/tp.p.downgoing.bin");
    //datawrite(par.data_remove_wave-=datavz,nt,npx,"./result/model4/tp.p.downgoing.err.bin");
    datatp3d.col(0)=par.data_remove_wave.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,fmax,frule,ncpu);
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.downgoing.vz.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
}
/////////////////////////////////////////////////////////////
//data code
if(0)
{
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.vz.bin");
    strcat(file1,numtostr(sx,5));
    dataread3d_bycol_transpose(datarecover3d,nz,nx,file1);
    
    cx_fcube datafx3d,datafx_code3d;
    datafx3d.copy_size(datarecover3d);
    datafx_code3d.copy_size(datarecover3d);
    for(i=0;i<nx;i++){
        if((i%2)!=0){
            //datarecover3d.row(i).fill(0);
        }
    }

    df=1.0/dt/nz;
    tx2fx_3d_thread(datafx3d,datarecover3d, ncpu);
    multiple_code3d(datafx_code3d, datafx3d, seabase_depth,\
        coordx, coordy,  velocity_water,\
        dx, dy, df, 1, fmax, ncpu);
    fx2tx_3d_thread(datarecover3d,datafx_code3d, ncpu);
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.code.vz.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    k=Beamforming_CG_3D(datatp3d,datarecover3d,dataerr3d,\
        datarecover3d,coordx,coordy, nz, nx,ny,dt,\
        npx,px0, dpx,npy,py0, dpy,fmax,frule,\
        ncpu, factor_L2, factor_L1, iterations_num, residual_ratio);
    datawrite3d_bycol_transpose(datatp3d,nz,npx,\
        "./result/model4/tp.upgoing.code.vz.bin");
}
/////////////////////////////////////////////////////////////
//demultiple 1
if(0)
{
    fmat subfmat1,subfmat2,datacode,datatp1,datatp2,coordx_local,filter;
    fmat subfmat1_t,subfmat2_t,datatpcode_t(npx,nz),datatpmulti_t(npx,nz),\
        datatpdemulti_t(npx,nz),datatpfilter_t(npx,nz);
    //datarecover_t=datarecover3d.col(0);
    //datacode=datarecover_t.t();
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.code.vz.bin");
    strcat(file1,numtostr(sx,5));
    datacode=dataread(nt,nx,file1);
    int ntrace_local(100),ny_local(1),trace1,trace2;
    float x0_local(-dx*ntrace_local/2);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.vz.bin");
    strcat(file1,numtostr(sx,5));
    dataread(data,nt,nx,file1);
    data=fmatsmooth(data,nt,nx,9);


    for(trace1=0;trace1<=nx-ntrace_local;trace1+=ntrace_local){
        trace2=min(nx,trace1+ntrace_local);
        subfmat1=getsubfmat(data,0,nt-1,trace1,trace2-1);
        subfmat2=getsubfmat(datacode,0,nt-1,trace1,trace2-1);
        subfmat1_t=subfmat1.t();
        subfmat2_t=subfmat2.t();
        cout<<"now is process trace:"<<trace2<<endl;

        int nx_local(trace2-trace1);
        float x0_local(-dx*nx_local/2);
        coordx_local.zeros(nx_local,ny_local);
        for(i=0;i<nx_local;i++){
        for(j=0;j<ny_local;j++){
            coordx_local(i,j)=dx*i+x0_local;
        }}

        k= Beamforming_CG_2D(datatpmulti_t,datarecover_t,dataerr_t,\
            subfmat1_t, coordx_local, nt, nx_local, dt,npx,px0, dpx,\
            fmax,frule,ncpu, factor_L2,factor_L1,\
            iterations_num, residual_ratio);
        k= Beamforming_CG_2D(datatpcode_t,datarecover_t,dataerr_t,\
            subfmat2_t, coordx_local, nt, nx_local, dt,npx,px0, dpx,\
            fmax,frule,ncpu, factor_L2,factor_L1,\
            iterations_num, residual_ratio);
        datatpdemulti_t.fill(0);
        get_taup2d_CCA(datatpdemulti_t, datatpfilter_t, datatpmulti_t, datatpcode_t,\
            ncpu,20, 20,1, 2, 99);

        k= Beamforming_recoverdata_2D(subfmat1_t,datatpdemulti_t,\
            coordx_local, nz, nx_local, dt, npx,px0, dpx, fmax, frule, ncpu);
        subfmat1=subfmat1_t.t();
        copydata2subfmat(data,subfmat1,0,nt-1,trace1,trace2-1);
    }
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.demulti.1.vz.bin");
    strcat(file1,numtostr(sx,5));
    datawrite(data,nt,nx,file1);
    //datawrite(datatp1=datatpdemulti_t.t(),nt,npx,"./result/model5/tp.upgoing.demulti.vz.bin");
}
/////////////////////////////////////////////////////////////
//demultiple 2
if(0)
{
    demultiple2d par;
    par.dt=dt,par.nt=nt,par.dx=dpx,par.nx=npx;
    par.data_remove_wave.zeros(par.nt,par.nx);
    par.data_antiremove_wave.copy_size(par.data_remove_wave);
    par.data2d.copy_size(par.data_remove_wave);
    par.datadict.copy_size(par.data_remove_wave);
    par.datad.copy_size(par.data_remove_wave);
    par.datah.copy_size(par.data_remove_wave);
    par.datahd.copy_size(par.data_remove_wave);
    par.nwp=(5),par.nwph=(5),par.nwpd=(5),par.nwphd=(5),\
    par.n1w=(100),par.d1w=(90),par.lmd=1;
    par.datawinbeg.zeros(par.nx,1);
    par.datawinend.copy_size(par.datawinbeg);

    par.datawinend.col(0)=datawin.col(par.nx/2);
    //datawrite(par.datawinend,npx,1,"datawin.bin");

    fmat datap,datavz;
    datap.copy_size(par.data_remove_wave);
    datavz.copy_size(par.data_remove_wave);

    dataread(datavz,"./result/model4/tp.p.upgoing.bin");
    dataread(datap,"./result/model4/tp.upgoing.code.vz.bin");
    datap=datap/(datap.max()-datap.min());
    datavz=datavz/(datavz.max()-datavz.min());
    datavz=get_blackman_upwin2d\
        (datavz,100);
    datap=get_blackman_upwin2d\
        (datap,100);

/////////////////////////////////////////////////////////////

omp_set_num_threads(ncpu);
#pragma omp parallel for
    for(k=0;k<par.nx;k++){
        fmat datal(nt,1);
        //cout<<k<<"|";
        par.datah.col(k)=hilbert1D(datal=datap.col(k),nt,dt);
        par.datad.col(k)=diff1D(datal=datap.col(k),nt,dt);
        par.datahd.col(k)=diff1D(datal=par.datah.col(k),nt,dt);
    }

    datawrite(par.datah,"./media/obndata.p2.h.bin");
    datawrite(par.datad,"./media/obndata.p2.d.bin");
    datawrite(par.datahd,"./media/obndata.p2.hd.bin");

    dataread(par.datah,"./media/obndata.p2.h.bin");
    dataread(par.datad,"./media/obndata.p2.d.bin");
    dataread(par.datahd,"./media/obndata.p2.hd.bin");

////////////////////////////////////////////////////////////////////
    par.data2d=datavz, par.datadict=datap;
    //single_trace_dewave(par,ncpu);
    cout<<par.data2d.n_rows<<"|"<<par.data2d.n_cols<<endl;
    datawrite(par.data_remove_wave,"./result/model4/tp.upgoing.demulti.vz.bin");
    datatp3d.col(0)=par.data_remove_wave.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,fmax,frule,ncpu);
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.upgoing.demulti.vz.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
/////////////////////////////////////////////////////////////
}
/////////////////////////////////////////////////////////////
//UDD
fmat dataw;
if(1)
{
    fmat datapr(nt,npx),dataup(nt,npx),datadown(nt,npx);
    cx_fmat dataupfft(nt,npx),datadownfft(nt,npx),\
        datadownfft_t(npx,nt),dataprfft(nt,npx,fill::zeros);

    dataread(dataup,nt,npx,"./result/model4/tp.p.upgoing.bin");
    dataread(datadown,nt,npx,"./result/model4/tp.p.downgoing.bin");
    datapr(0,0)=dataup.max()*0.1;
    for(k=0;k<npx;k++){
        dataupfft.col(k)=fft(dataup.col(k));
        datadownfft.col(k)=fft(datadown.col(k));   
    }
    datadownfft_t=datadownfft.t();
    for(k=0;k<npx;k++){
        for(i=0;i<nt/2;i++){
            dataprfft(i,k)=(dataupfft(i,k)*datadownfft_t(k,i))/\
                (datadownfft(i,k)*datadownfft_t(k,i)+datapr(0,0));
        }
    }
    for(k=0;k<npx;k++){
        datapr.col(k)=real(ifft(dataprfft.col(k)));
    }
    datapr=get_blackman_upwin2d(datapr,500);
    for(k=8000;k<nz;k++){
        datapr.row(k).fill(0);
    }
    
    file1[0]='\0';
    strcat(file1,"./result/model4/tp.pr.UDD.bin");
    strcat(file1,numtostr(sx,5));
    datawrite(datapr,nt,npx,file1);
    dataw=datapr;

    df=1.0/dt/nz;
    datatp3d.col(0)=datapr.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,(nz/2-3)*df,frule,ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.rebuild.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    cx_fcube datafx3d,datafx_code3d;
    datafx3d.copy_size(datarecover3d);
    datafx_code3d.copy_size(datarecover3d);
    datarecover3d.col(0)=get_blackman_leftwin2d(datarecover3d.col(0), 200);
    datarecover3d.col(0)=get_blackman_rightwin2d(datarecover3d.col(0), 200);

    float wavelet_delay(0.05);
    fmat tu2u1(nx,nx);
    for(i=0;i<nx;i++){
        for(j=0;j<nx;j++){
            float deepth=j;
            //deepth=dz*(20+20/(1+exp(0.015*(450-deepth))));
            deepth=dz*35;
            float d=abs(i-j)*dx;
            d=sqrt(d*d+deepth*deepth);
            tu2u1(i,j)=d/velocity_water+wavelet_delay;
        }
    }
    datapr.zeros(nx,nz);
    datapr=datarecover3d.col(0);
    datapr=multiple_code2d(datapr.st(),tu2u1,df,10*df,(nz/2-3)*df,ncpu);
    datarecover3d.col(0)=datapr.st();
    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.rebuild.code.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    datapr.zeros(nx,nz);
    datapr=datarecover3d.col(0);
    fmat datapr_pro=fmatsmooth(datapr,nx,nz,25);
    for(k=8000;k<nz;k++){
        datapr_pro.col(k).fill(0);
        datapr.col(k).fill(0);
    }
    datapr_pro=datapr-datapr_pro;

    fmat wavelet(nz,1,fill::zeros);
    cx_fmat waveletfft(nz,1,fill::zeros);
    for(k=0;k<nz;k++){
        wavelet(k,0)=10*wavelet01(k,dt,30.0);
    }
    waveletfft.col(0)=fft(wavelet.col(0));
    wavelet.fill(0);
    //waveletfft.set_imag(wavelet);
    dataprfft.zeros(nx,nz);
    for(k=0;k<nx;k++){
        dataprfft.row(k)=(fft(datapr_pro.row(k)));
        for(i=0;i<nz/2;i++){
            dataprfft(k,i)=dataprfft(k,i)*waveletfft(i,0);
        }
        datapr_pro.row(k)=real(ifft(dataprfft.row(k)));
    }
    datapr_pro(span(0,nx-1),span(0,8000))\
        =datapr_pro(span(0,nx-1),span(0+167,8000+167));
    for(k=0;k<10;k++){
        datapr_pro.row(k).fill(0);
        datapr_pro.row(nx-k-1).fill(0);
    }
    datapr_pro=get_blackman_upwin2d(datapr_pro,50);
    datapr_pro=get_blackman_downwin2d(datapr_pro,50);

    datarecover3d.col(0)=datapr_pro;

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.rebuild.pro.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
}
/////////////////////////////////////////////////////////////
//UDD timespace inv
if(0)
{
    fmat datapr(nt,npx),dataup(1*nt,npx),datadown(1*nt,npx);

    dataread(dataup,nt,npx,"./result/model4/tp.p.upgoing.bin");
    dataread(datadown,nt,npx,"./result/model4/tp.p.downgoing.bin");
    /*
    for(k=0;k<nt;k++){
        dataup.row(nt+k)=dataup.row(k);
        dataup.row(2*nt+k)=dataup.row(k);
        datadown.row(nt+k)=datadown.row(k);
        datadown.row(2*nt+k)=datadown.row(k);
    }*/
    demultiple2d par;
    par.dt=dt,par.nt=1*nt,par.dx=dpx,par.nx=npx;
    par.data_remove_wave.zeros(par.nt,par.nx);
    par.data_antiremove_wave.copy_size(par.data_remove_wave);
    par.data2d.copy_size(par.data_remove_wave);
    par.datadict.copy_size(par.data_remove_wave);
    par.datad.copy_size(par.data_remove_wave);
    par.datah.copy_size(par.data_remove_wave);
    par.datahd.copy_size(par.data_remove_wave);
    par.nwp=(2000),par.nwph=(0),par.nwpd=(0),par.nwphd=(0),\
    par.n1w=(1*nt-5),par.d1w=(nt),par.lmd=10;
    par.datawinbeg.zeros(par.nx,1);
    par.datawinend.copy_size(par.datawinbeg);

    par.data2d=dataup, par.datadict=datadown;
    //dataw*=0;
    dataw=single_trace_dewave(par,dataw,ncpu);

    datapr.fill(0);
    for(k=0;k<par.nwp;k++){
        datapr.row(k)=dataw.row(k);
    }

    datapr=get_blackman_upwin2d(datapr,500);
    file1[0]='\0';
    strcat(file1,"./result/model4/tp.pr.UDD.time.bin");
    strcat(file1,numtostr(sx,5));
    datawrite(datapr,nt,npx,file1);

    df=1.0/dt/nz;
    datatp3d.col(0)=datapr.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,(nz/2-3)*df,frule,ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.time.rebuild.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    cx_fcube datafx3d,datafx_code3d;
    datafx3d.copy_size(datarecover3d);
    datafx_code3d.copy_size(datarecover3d);
    datarecover3d.col(0)=get_blackman_upwin2d(datarecover3d.col(0), 100);
    datarecover3d.col(0)=get_blackman_downwin2d(datarecover3d.col(0), 100);

    float wavelet_delay(0.0);
    tx2fx_3d_thread(datafx3d,datarecover3d, ncpu);
    multiple_code3d(datafx_code3d, datafx3d, seabase_depth,\
        coordx, coordy,  velocity_water,dx, dy, df, 1*df, \
        (nz/2-3)*df, ncpu, wavelet_delay);
    fx2tx_3d_thread(datarecover3d,datafx_code3d, ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.time.rebuild.code.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
}
/////////////////////////////////////////////////////////////
//DGD
if(1)
{
    fmat datapr(nt,npx,fill::zeros),datadown_s(nt,npx),datadown(nt,npx);
    cx_fmat datadown_sfft(nt,npx),datadownfft(nt,npx),\
        datadownfft_t(npx,nt),dataprfft(nt,npx,fill::zeros);

    dataread(datadown,nt,npx,"./result/model4/tp.p.downgoing.bin");
    datadown_s=datadown;
    datapr(0,0)=datadown.max()*0.1;
    datapr(1,1)=-1;

    for(i=0;i<npx;i++){
        for(j=datawin(i,npx/2)-150;j<=datawin(i,npx/2)+150;j++){
            datadown_s(j,i)=0;
        }
    }
    datawrite(datadown_s,nt,npx,"./result/model4/tp.without.s.downgoing.bin");
    datadown*=-1;

    for(k=0;k<npx;k++){
        datadown_sfft.col(k)=fft(datadown_s.col(k));
        datadownfft.col(k)=fft(datadown.col(k));   
    }
    datadownfft_t=datadownfft.t();
    for(k=0;k<npx;k++){
        for(i=0;i<nt/2;i++){
            dataprfft(i,k)=(datadown_sfft(i,k)*datadownfft_t(k,i))/\
                ((datadownfft(i,k))*datadownfft_t(k,i)+datapr(0,0));
        }
    }
    for(k=0;k<npx;k++){
        datapr.col(k)=real(ifft(dataprfft.col(k)));
    }
    datapr=get_blackman_upwin2d(datapr,500);
    for(k=8000;k<nz;k++){
        datapr.row(k).fill(0);
    }

    file1[0]='\0';
    strcat(file1,"./result/model4/tp.pr.DGR.bin");
    strcat(file1,numtostr(sx,5));
    datawrite(datapr,nt,npx,file1);

    df=1.0/dt/nz;
    datatp3d.col(0)=datapr.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,(nz/2-3)*df,frule,ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.DGR.rebuild.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    cx_fcube datafx3d,datafx_code3d;
    datafx3d.copy_size(datarecover3d);
    datafx_code3d.copy_size(datarecover3d);
    datarecover3d.col(0)=get_blackman_leftwin2d(datarecover3d.col(0), 200);
    datarecover3d.col(0)=get_blackman_rightwin2d(datarecover3d.col(0), 200);

    float wavelet_delay(0.05);
    fmat tu2u1(nx,nx);
    for(i=0;i<nx;i++){
        for(j=0;j<nx;j++){
            float deepth=j;
            //deepth=dz*(20+20/(1+exp(0.015*(450-deepth))));
            deepth=dz*35;
            float d=abs(i-j)*dx;
            d=sqrt(d*d+deepth*deepth);
            tu2u1(i,j)=d/velocity_water+wavelet_delay;
        }
    }
    datapr.zeros(nx,nz);
    datapr=datarecover3d.col(0);
    datapr=multiple_code2d(datapr.st(),tu2u1,df,1*df,(nz/2-3)*df,ncpu);
    datarecover3d.col(0)=datapr.st();

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.DGR.rebuild.code.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
}
/////////////////////////////////////////////////////////////
//UGD
if(0)
{
    fmat datapr(nt,npx,fill::zeros),dataup_adds(nt,npx,fill::zeros),\
        datadown(nt,npx),dataup(nt,npx);
    cx_fmat dataup_sfft(nt,npx),dataupfft(nt,npx),\
        dataupfft_adds_t(npx,nt),dataprfft(nt,npx,fill::zeros);

    dataread(datadown,nt,npx,"./result/model4/tp.p.downgoing.bin");
    dataread(dataup,nt,npx,"./result/model4/tp.p.upgoing.bin");
    dataup_adds=(-1)*dataup;
    datapr(0,0)=datadown.max()*0.1;
    datapr(1,1)=-1;

    for(i=0;i<npx;i++){
        for(j=datawin(i,npx/2)-150;j<=datawin(i,npx/2)+150;j++){
            dataup_adds(j,i)+=datadown(j,i);
        }
    }
    datawrite(dataup_adds,nt,npx,"./result/model4/tp.add.s.upgoing.bin");

    for(k=0;k<npx;k++){
        dataup_sfft.col(k)=fft(dataup_adds.col(k));
        dataupfft.col(k)=fft(dataup.col(k));   
    }
    dataupfft_adds_t=dataup_sfft.t();
    for(k=0;k<npx;k++){
        for(i=0;i<nt/2;i++){
            dataprfft(i,k)=(dataupfft(i,k)*dataupfft_adds_t(k,i))/\
                ((dataup_sfft(i,k))*dataupfft_adds_t(k,i)+datapr(0,0));
        }
    }
    for(k=0;k<npx;k++){
        datapr.col(k)=real(ifft(dataprfft.col(k)));
    }
    datapr=get_blackman_upwin2d(datapr,500);
    for(k=8000;k<nz;k++){
        datapr.row(k).fill(0);
    }

    file1[0]='\0';
    strcat(file1,"./result/model4/tp.pr.UGR.bin");
    strcat(file1,numtostr(sx,5));
    datawrite(datapr,nt,npx,file1);

    df=1.0/dt/nz;
    datatp3d.col(0)=datapr.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,(nz/2-3)*df,frule,ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UGR.rebuild.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    cx_fcube datafx3d,datafx_code3d;
    datafx3d.copy_size(datarecover3d);
    datafx_code3d.copy_size(datarecover3d);
    datarecover3d.col(0)=get_blackman_leftwin2d(datarecover3d.col(0), 50);
    datarecover3d.col(0)=get_blackman_rightwin2d(datarecover3d.col(0), 50);

    float wavelet_delay(0.05);
    fmat tu2u1(nx,nx);
    for(i=0;i<nx;i++){
        for(j=0;j<nx;j++){
            float deepth=j;
            //deepth=dz*(20+20/(1+exp(0.015*(450-deepth))));
            deepth=dz*35;
            float d=abs(i-j)*dx;
            d=sqrt(d*d+deepth*deepth);
            tu2u1(i,j)=d/velocity_water+wavelet_delay;
        }
    }
    datapr.zeros(nx,nz);
    datapr=datarecover3d.col(0);
    datapr=multiple_code2d(datapr.st(),tu2u1,df,1*df,(nz/2-3)*df,ncpu);
    datarecover3d.col(0)=datapr.st();

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UGR.rebuild.code.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
}
/////////////////////////////////////////////////////////////
//DGD timespace inv
if(0)
{
    fmat datapr(nt,npx),dataup(3*nt,npx),datadown(3*nt,npx),dataw;

    dataread(datadown,nt,npx,"./result/model4/tp.p.downgoing.bin");
    dataup=datadown;
    for(i=0;i<npx;i++){
        for(j=datawin(i,npx/2)-100;j<=datawin(i,npx/2)+100;j++){
            dataup(j,i)=0;
        }
    }
    for(k=0;k<nt;k++){
        dataup.row(nt+k)=dataup.row(k);
        dataup.row(2*nt+k)=dataup.row(k);
        datadown.row(nt+k)=datadown.row(k);
        datadown.row(2*nt+k)=datadown.row(k);
    }
    demultiple2d par;
    par.dt=dt,par.nt=3*nt,par.dx=dpx,par.nx=npx;
    par.data_remove_wave.zeros(par.nt,par.nx);
    par.data_antiremove_wave.copy_size(par.data_remove_wave);
    par.data2d.copy_size(par.data_remove_wave);
    par.datadict.copy_size(par.data_remove_wave);
    par.datad.copy_size(par.data_remove_wave);
    par.datah.copy_size(par.data_remove_wave);
    par.datahd.copy_size(par.data_remove_wave);
    par.nwp=(5000),par.nwph=(0),par.nwpd=(0),par.nwphd=(0),\
    par.n1w=(3*nt-5),par.d1w=(nt),par.lmd=500;
    par.datawinbeg.zeros(par.nx,1);
    par.datawinend.copy_size(par.datawinbeg);

    par.data2d=dataup, par.datadict=datadown;
    //dataw=single_trace_dewave(par,ncpu);

    datapr.fill(0);
    for(k=0;k<par.nwp;k++){
        datapr.row(k)=dataw.row(k);
    }

    datapr=get_blackman_upwin2d(datapr,500);
    file1[0]='\0';
    strcat(file1,"./result/model4/tp.pr.UDD.bin");
    strcat(file1,numtostr(sx,5));
    datawrite(datapr,nt,npx,file1);

    df=1.0/dt/nz;
    datatp3d.col(0)=datapr.t();
    k= Beamforming_recoverdata_3D(datarecover3d,datatp3d,\
        coordx,coordy, nz, nx,ny,dt,npx,px0, dpx,npy,py0,\
        dpy,(nz/2-3)*df,frule,ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.rebuild.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);

    cx_fcube datafx3d,datafx_code3d;
    datafx3d.copy_size(datarecover3d);
    datafx_code3d.copy_size(datarecover3d);
    datarecover3d.col(0)=get_blackman_upwin2d(datarecover3d.col(0), 100);
    datarecover3d.col(0)=get_blackman_downwin2d(datarecover3d.col(0), 100);

    float wavelet_delay(0.0);
    tx2fx_3d_thread(datafx3d,datarecover3d, ncpu);
    multiple_code3d(datafx_code3d, datafx3d, seabase_depth,\
        coordx, coordy,  velocity_water,dx, dy, df, 1*df, \
        (nz/2-3)*df, ncpu, wavelet_delay);
    fx2tx_3d_thread(datarecover3d,datafx_code3d, ncpu);

    file1[0]='\0';
    strcat(file1,"./result/model4/tx.pr.UDD.rebuild.code.bin");
    strcat(file1,numtostr(sx,5));
    datawrite3d_bycol_transpose(datarecover3d,nz,nx,file1);
}
/////////////////////////////////////////////////////////////
}
    cx_fmat test(1,1);
    test(0,0).real(1);
    test(0,0).imag(2);
    test.print();
    test.fill(0);
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




