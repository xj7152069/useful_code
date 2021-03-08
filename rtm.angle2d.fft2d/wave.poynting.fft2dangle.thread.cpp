
#include <iostream>
#include <thread>
#include <future>
#include "xjc.h"
using namespace std; 

void test(int par1, int par2, int par3);
int main()
{
    int par_num(6),par_beg(50),par_gep(50),par_gep_in(10);
    int i,j;
    void *p=operator new(par_num*sizeof(thread));
    thread *t=static_cast<thread *>(p);
    for(i=0;i<par_num;i++)
    {
        new(t+i) thread(test,par_beg+i*par_gep,par_beg+(i+1)*par_gep-par_gep_in,par_gep_in);
    }
    for(i=0;i<par_num;i++)
    {
        (t+i)->join();
    }
    delete []p;
    /*/����һ���߳� 
    std::thread t1(test,100);
    std::thread t2(test,200);
    std::thread t3(test,300);
    t1.join();
    t2.join(); 
    t3.join(); */
    std::cout << "All thread finished\n";
    
    return 0;
}

void test(int par1, int par2, int par3)
{
    int source_beg, source_end, source_gep;
    source_beg=par1;
    source_end=par2;
    source_gep=par3;

    int T,Z,Zfft,X,Xfft,sx,ds,s_id,s_num,zf;
    T=2000;
    Z=200;
    X=401;
    Zfft=256;
    Xfft=512;
    zf=35;
    sx=250;

    int i,j,k;
    float *s, *sh;
    fmat suf(T,X),suf_v(T,X),b(Z,X),d(Z,X);
    fmat sufh(T,X),sufh_v(T,X),ns(T,X),nsh(T,X);
    char file1[99];
    ofstream outfs, outfsh, outfup, outfdown;
    wave2D S(Z,X), Sh(Z,X);
    fmat cp(Z,X),cpfft(Zfft,Xfft),cpfft2(Zfft,Xfft);
    cx_fmat comwave(Z,X), wfft(Zfft,Xfft), wfft2(Zfft,Xfft), wifft(Zfft,Xfft);

    float **vx,**vz,**swave,**rwave,**imag,a_beg(-200.0),a_gep(2.0),na(200),a_err(15.0),*waveangle;
    int nz,nx,zfft,xfft,nwin(37);
    nz=Z; nx=X; zfft=Zfft; xfft=Xfft;
    Poynting_2D pytS(3,nz,nx);
    Poynting_2D pytR(3,nz,nx);
    fft2d_anglegather2d test(na, nz, nx, Zfft, Xfft, a_beg, a_gep);
    swave=newfmat(zfft,xfft);
    rwave=newfmat(zfft,xfft);
    imag=newfmat(nz,nx);
    vx=newfmat(nz,nx);
    vz=newfmat(nz,nx);
    ifstream infs,infr,infrs,infrsh;
    ofstream outffacig,outfdacig,outfimag,outfrecvacig;
    ofstream outstheta,outrtheta;

    fmat w0(Xfft,Zfft);
    w0=windowsZ(Xfft,Zfft,5.0);
    fmat w(Zfft,Xfft);
    w=w0.t();
    //w=fmatsmooth(w,Zfft,Xfft,5);
    datawrite(w,Zfft,Xfft,"window2d.bin");
    waveangle=new float[nwin];
    for(k=0;k<nwin;k++)
    {
        waveangle[k]=-90+k*5;
        //cout<<"{"<<waveangle[k]<<"} "<<endl;
    }
    test.createanglewin(nwin,a_err,waveangle);
    test.creatRecvAngle(nwin,-90.0,5.0);

    s=new float[T];
    wavelet01(s, T, S.dt);
    sh=hilbert1D(s, T, S.dt);

    sx=source_beg-source_gep;
    while(sx<source_end)
    {
        sx+=source_gep;

        S.cleardata();
        Sh.cleardata();
        dataread(S.p2,Z,X,"./data/model.dat");  
        dataread(Sh.p2,Z,X,"./data/model.dat");  
        //outfs.open("./data/zy.movie.s.bin");
        //outfsh.open("./data/zy.movie.sh.bin");
        //outfup.open("./data/zy.movie.up.bin");
        file1[0]='\0';
        strcat(file1,"./data/zy.movie.down.bin");
        strcat(file1,numtostr(source_beg,5));
        outfdown.open(file1);
        for(k=0;k<T;k++)
            {
            S.s2[zf][sx]+=s[k];
            S.timeslicecal();
            S.timeslicecopy();
            Sh.s2[zf][sx]+=sh[k];
            Sh.timeslicecal();
            Sh.timeslicecopy();
            
            //datawrite(S.s3, Z, X, outfs);
            //datawrite(Sh.s3, Z, X, outfsh);

            cp=matcopy(S.s3,Z,X);
            comwave.set_real(cp);
            cp=matcopy(Sh.s3,Z,X);
            comwave.set_imag(cp);
            wfft=fft2(comwave,Zfft,Xfft);

            cpfft=real(wfft);
            cpfft2=matmul(cpfft,w,Zfft,Xfft);
            wfft.set_real(cpfft2);
            cpfft=arma::imag(wfft);
            cpfft2=matmul(cpfft,w,Zfft,Xfft);
            wfft.set_imag(cpfft2);

            cpfft.fill(0.0);
            wfft2.set_real(cpfft);
            wfft2.set_imag(cpfft);
            for(i=Zfft/2;i<Zfft;i++)
            {
                wfft2.row(i)=wfft.row(i);
            }
            wifft=ifft2(wfft2,Zfft,Xfft);
            cpfft=real(wifft);
            datawrite(cpfft,Zfft,Xfft,outfdown);

            b=matcopy(S.s3,Z,X);
            suf.row(k)=b.row(zf);
            b=matcopy(Sh.s3,Z,X);
            sufh.row(k)=b.row(zf);
            if(k%500==0)
                cout<<k<<endl;
            }
        //outfs.close();
        //outfsh.close();
        //outfup.close();
        outfdown.close();

        S.cleardata();
        S.setvelocity(3000.0);
        Sh.cleardata();
        //wave2D C(Z,X);
        Sh.setvelocity(3000.0);
        for(k=0;k<T;k++)
        {
            S.s2[zf][sx]+=s[k];
            S.timeslicecal();
            S.timeslicecopy();
            Sh.s2[zf][sx]+=sh[k];
            Sh.timeslicecal();
            Sh.timeslicecopy();

            b=matcopy(S.s3,Z,X);
            suf_v.row(k)=b.row(zf);
            b=matcopy(Sh.s3,Z,X);
            sufh_v.row(k)=b.row(zf);
            if(k%500==0)
                cout<<k<<endl;
        }
        S.cleardata();
        Sh.cleardata();
        dataread(S.p2,Z,X,"./data/model.dat");  
        dataread(Sh.p2,Z,X,"./data/model.dat");  
        matsmooth(S.p2,S.p2,Z,X,25);
        matsmooth(Sh.p2,Sh.p2,Z,X,25);
        datawrite(S.p2,Z,X,"./data/model.smooth.dat");

        ns=suf-suf_v;
        //datawrite(ns,T,X,"./data/ns_suf.bin");
        nsh=sufh-sufh_v;
        //datawrite(nsh,T,X,"./data/nsh_suf.bin");

        file1[0]='\0';
        strcat(file1,"./data/ns.movie.s.bin");
        strcat(file1,numtostr(source_beg,5));
        outfs.open(file1);
        file1[0]='\0';
        strcat(file1,"./data/ns.movie.sh.bin");
        strcat(file1,numtostr(source_beg,5));
        outfsh.open(file1);
        //outfdown.open("./data/ns.movie.down.bin");
        for(k=0;k<T;k++)
            {
            d=matcopy(S.s2,Z,X);
            d.row(zf)=d.row(zf)+ns.row(T-1-k);
            matcopy(S.s2,d,Z,X);
            S.timeslicecal();
            S.timeslicecopy();
            d=matcopy(Sh.s2,Z,X);
            d.row(zf)=d.row(zf)+nsh.row(T-1-k);
            matcopy(Sh.s2,d,Z,X);
            Sh.timeslicecal();
            Sh.timeslicecopy();
            datawrite(S.s3, Z, X, outfs);
            datawrite(Sh.s3, Z, X, outfsh);

            if(k%500==0)
                cout<<k<<endl;
            }
        outfs.close();
        outfsh.close();
        //outfup.close();
        //outfdown.close();

    /////////////////////////////////angle calculation//////////////////////////////////

        //matcopy(p,0.0,nz,nx);
        matcopy(imag,0.0,nz,nx);
        matcopy(vx,0.0,nz,nx);
        matcopy(vz,0.0,nz,nx);

        file1[0]='\0';
        strcat(file1,"./data/zy.movie.down.bin");
        strcat(file1,numtostr(source_beg,5));
        infs.open(file1);
        cout<<file1<<endl;
        file1[0]='\0';
        strcat(file1,"./data/ns.movie.s.bin");
        strcat(file1,numtostr(source_beg,5));
        infrs.open(file1);
        file1[0]='\0';
        strcat(file1,"./data/ns.movie.sh.bin");
        strcat(file1,numtostr(source_beg,5));
        infrsh.open(file1);
        //outstheta.open("./data/thetaS.down.bin");
        //outrtheta.open("./data/thetaR.up.bin");
        //out3.open("./data/wave.ns.pyt");
        pytS.dataClear();
        test.cleardata();
        
        for(k=100;k<T-500;k++)
        {
            infrsh.seekg(Z*X*(T-k-1)*sizeof(float), ios::beg);
            infrs.seekg(Z*X*(T-k-1)*sizeof(float), ios::beg);
            cp=dataread(Z,X,infrs);
            comwave.set_real(cp);
            //datawrite(cp,Z,X,"waveformr.bin");
            cp=dataread(Z,X,infrsh);
            comwave.set_imag(cp);
            
            infs.seekg(zfft*xfft*(k)*sizeof(float), ios::beg);
            dataread(swave,zfft,xfft,infs);
            //datawrite(swave,zfft,xfft,"waveforms.bin");
            pytS.addtimeslicecal(swave);     
            pytS.velocityCalculate();

            matsmooth(vz,pytS.vz,nz,nx,1);
            matsmooth(vx,pytS.vx,nz,nx,1);
            matmul(vz,-1.0,nz,nx);
            matmul(vx,-1.0,nz,nx);
            test.theatcal_S(vz, vx);
            if(k==600)
            {datawrite(test.ps,Z,X,"pyt.theta.bin");}
            
            //datawrite(test.pr,nz,nx,outrtheta);
            //datawrite(test.ps,nz,nx,outstheta);
            
            test.addFA_DA(swave, comwave);
            
            if(k%100==0)
                cout<<k<<endl;
        }
        infrs.close();
        infrsh.close();
        infs.close();
        //outstheta.close();

        file1[0]='\0';
        strcat(file1,"./data/fa.angle.cig.gather");
        strcat(file1,numtostr(float(sx),5));
        outffacig.open(file1);
        file1[0]='\0';
        strcat(file1,"./data/da.angle.cig.gather");
        strcat(file1,numtostr(float(sx),5));
        outfdacig.open(file1);
        file1[0]='\0';
        strcat(file1,"./data/recv.angle.cig.gather");
        strcat(file1,numtostr(float(sx),5));
        outfrecvacig.open(file1);
        for(i=0;i<nx;i++)
        {
            datawrite(test.FA[i],nz,na,outffacig);
            datawrite(test.DA[i],nz,na,outfdacig);
            datawrite(test.recvA[i],nz,nwin,outfrecvacig);
            //cout<<"output the CIG"<<endl;
            for(j=0;j<nz;j++)
            {
                for(k=0;k<na;k++)
                {
                    imag[j][i]+=test.FA[i][j][k];
                }
            }
        }
        file1[0]='\0';
        strcat(file1,"./data/fa.angle.imag");
        strcat(file1,numtostr(float(sx),5));
        datawrite(imag,nz,nx,file1);
        matcopy(imag,0.0,nz,nx);
        outfdacig.close();
        outffacig.close();
        outfrecvacig.close();

        cout<<"finshed --> "<<sx<<endl;
    }
        //return 0;
}

