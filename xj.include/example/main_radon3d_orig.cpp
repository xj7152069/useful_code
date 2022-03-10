
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
void debadtrace(fmat & data);
int main()
{
    int i,j,k,usecopy(0);
    int nz(1024),nx(200),ny(200),nf(1024);
    //int nz(1024),nx(200),ny(1),nf(1024);
    //int nz(2048),nx(200),ny(1),nf(2048);
    //int nz(3000),nx(101),ny(1),nf(nz);
    ofstream outf1,outf2;
    ifstream inf1,inf2;
    //inf1.open("./data/agcdata.3000.974.bin");
    //inf1.open("./data/txz.three.1024.200.bin");
    //inf1.open("./data/txz9.bin");
    //inf1.open("./data/txz.two.2048.200.bin");
    inf1.open("./data/data3d.two.1024.200.200.bin");
//注意处理实际数据的时候，可能一个测线包含的道数不相等;
//进行三维处理时需要重新调整;
    int nh(nx*ny);
    fmat basedata(nz,nh);
//读取数据并加入随机噪声，实际数据可以不用加
    dataread(basedata,nz,nh,inf1);
    inf1.close();
    //basedata=10*basedata/(max(max(basedata)));
    if(usecopy==1)
    {
        fmat basedataold(nz,nh);
        basedataold=basedata;
        nz=nz*2;nf=nz;
        basedata.resize(nz,nh);
        for(i=0;i<nz/2;i++)
        {
            basedata.row(i)=basedataold.row(i);
            basedata.row(nz-i-1)=basedataold.row(i);
        }
    }

    struct linerradon3d par;
    //struct beamformingCG3d parcg;
    beamforming_parset(nx,ny,nz,par);
    cout<<"npx="<<par.npx<<endl;
    cout<<"dpx="<<par.dpx<<endl;
    cout<<"npy="<<par.npy<<endl;
    cout<<"dpy="<<par.dpy<<endl;
    
    //cin>>k;
    par.npx=51;
    par.npy=51;
    //par.dpx*=6;
    par.dpx=0.00001;
    par.dpy=0.00001;
    par.px_center=0.00000*20;
    par.py_center=0.0;   
    par.numthread=6;
    par.dig_n=1;
    beamforming_parupdate(par);
    par.nf1=1;
    par.nf2=90;
    par.rulef1=1;
    par.rulef2=90;
/////////////////////////////////////////////////
    cout<<"par.px0 = "<<par.px_center<<endl;
    cout<<"par.npx = "<<par.npx<<endl;
    cout<<"par.dpx = "<<par.dpx<<endl;
    cout<<"par.py0 = "<<par.py_center<<endl;
    cout<<"par.npy = "<<par.npy<<endl;
    cout<<"par.dpy = "<<par.dpy<<endl;
    cout<<"par.nf2 = "<<par.nf2<<endl;
    cout<<"dig="<<par.dig_n<<endl;
    //cin>>k;

    fcube backupdata;
    fmat dataline(nz,nx,fill::zeros);
    backupdata.copy_size(par.data);
    k=0;
    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
        {
            dataline.col(j)=basedata.col(k);
            k++;
        }
        par.data.col(i)=dataline.t();
        backupdata.col(i)=dataline.t();
        //par.data.col(i)=get_blackman_rightwin2d\
            (par.data.col(i),par.nz/10.0);
        datawrite(dataline,nz,nx,"data.bin");
        dataline.fill(0.0);
    }
    cout<<"step - 0 is ok"<<endl;
////////////////////////////////////////////////
    linerradon(par); 
    par.dig_n1=0;
    par.dig_n2=400;
    //beamforming_cleardata(par);
    beamforminginv3d(par);

    //beamformingCG3d(par,15);
    cout<<"step - 1 is ok"<<endl;
////////////////////////////////////////////////

    par.p_power.fill(0);
    for(i=0;i<par.nz;i++)
    {
    par.p_power+=abs(par.realdataTP.slice(i));
    }
    datawrite(par.p_power,par.npx,par.npy,\
        "./result/data3d.p_power.bin");
///////////////////////////////////////////////

    fmat tpdata(nz,par.npx),tpdata_t(par.npx,nz);
    outf1.open("./result/data3d.tp.bin");
    //datawrite(par.p_power,par.npx,par.npy,\
        "./result/data3d.p_power.bin");
    for(k=0;k<par.npy;k++)
    {
        tpdata_t=par.realdataTP.col(k);
        tpdata_t.col(0)=(1.0/50)*par.p_power.col(k);
        datawrite(tpdata=tpdata_t.t(),nz,par.npx,outf1);
    }
    outf1.close();
/////////////////////////////////////////////////////////
    //cout<<"step - beamformingCG3d : "<<endl;

    rebuildsignal(par);

    for(k=1023;k<nz;k++)
    {
        //par.realrebuildtx.slice(k).fill(0);
        backupdata.slice(k).fill(0);
    }
    outf1.open("./result/data3d.rebuild.bin");
    outf2.open("./result/data3d.err.bin");
    fmat filterdata_t(nx,nz), db(1,1,fill::zeros);
    for(k=0;k<0;k++)
    {
        par.realrebuildtx.row(k).fill(0.001);
        par.realrebuildtx.row(nx-k-1).fill(0.001);
        backupdata.row(k).fill(0.001000);
        backupdata.row(nx-k-1).fill(0.001000);
    }
    for(k=0;k<0;k++)
    {
        par.realrebuildtx.col(k).fill(0.001);
        par.realrebuildtx.col(ny-1-k).fill(0.001);
        backupdata.col(k).fill(0.001000);
        backupdata.col(ny-1-k).fill(0.001000);
    }

    for(k=0;k<ny;k++)
    {
        filterdata_t=par.realrebuildtx.col(k);
        datawrite(dataline=filterdata_t.t(),par.nz,par.nx,outf1);
        filterdata_t=backupdata.col(k);
        datawrite(dataline=dataline-filterdata_t.t(),\
            par.nz,par.nx,outf2);
        filterdata_t=par.realrebuildtx.col(k);
        db=db+(sum(sum(abs(dataline)))/sum(sum(abs(filterdata_t))));
    }
    outf1.close();
    outf2.close();
    db=db/ny;
    db=20.0*log10(db);
    cout<<"step - 2 is ok: "<<db(0,0)<<" DB"<<endl;

    //int nn(par.npx*par.npy);
    //fmat hess1(nn,nn),hess2(nn,nn);
    //hess1=dataread(nn,nn,"hessmat.all.bin");
    //hess2=dataread(nn,nn,"hessmat.simple.bin");
    //datawrite(hess1=(hess1-hess2),nn,nn,"hessmat.err.bin");

    return 0;
}

void gyh(fmat & data)
{
    float maxpower(-1),minpower(1),abspower;
    int i,j;
    maxpower=max(max(data));  
    minpower=min(min(data));  
    abspower=maxpower-minpower;
    cout<<maxpower<<"||"<<minpower<<"||"<<abspower<<"||"<<endl;
    for(i=0;i<data.n_rows;i++)
    {
        for(j=0;j<data.n_cols;j++)
        {
            data(i,j)=data(i,j)/abspower;
        }
    }
}

void getrectangledata(fmat & data, int row1, int row2,\
    int col1, int col2)
{
    int i,j;
    for(i=0;i<data.n_rows;i++)
    {
    if(i<row1 || i>row2)
    {
        data.row(i).fill(0);
    }
    }
    for(i=0;i<data.n_cols;i++)
    {
    if(i<col1 || i>col2)
    {
        data.col(i).fill(0);
    }
    }
}

fmat getrectangledata(fmat & data, float perc, int n3, int n4)
{
    int i,j,i1,j1,n1,n2;
    n1=data.n_rows,n2=data.n_cols;
    fmat data2(n1,n2),win(n1,n2);
    data2=matmul(data,data,n1,n2);
    win.fill(0.0);
    float maxpower;
    maxpower=max(max(data2));

    for(i=0;i<n1;i++)
    {
        for(j=0;j<n2;j++)
        {
            if(data2(i,j)>=maxpower*perc)
            {win(i,j)=1;}
            else if(i>n3 && i<n1-n3 && j<n2-n4 && j>n4)
            {
                for(i1=i-n3;i1<i+n3;i1++)
                {
                    for(j1=j-n4;j1<j+n4;j1++)
                    {
                        if(data2(i1,j1)>=maxpower*perc)
                        {
                            win(i,j)=1;
                        }
                    }
                }
            }
        }
    }
    return win;
}





