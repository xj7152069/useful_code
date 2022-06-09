
clear
% float v,dx,dy,dpx,dpy,dt,zr,xr,yr;
% int nt,nx,ny,npx,npy,i,j,k,kx,ky,kpx,kpy;

v=1500;

nx=20001;
ny=20001;
dx=0.1;
dy=0.1;
dt=0.0003;
npx=501;
npy=10;
dpx=0.000003;
dpy=0.000003;
px0=-dpx*npx/2;
py0=-dpy*npy/2;

zr=200;
xr=0;
yr=0;

datatime0=zeros(npx,npy);
for kx=1:nx
    xr=kx*dx;
    for ky=1:ny
        yr=ky*dy;
        t=(1/v)*sqrt(zr*zr+xr*xr+yr*yr);
        px=(1/v)*(xr/sqrt(zr*zr+xr*xr+yr*yr));
        py=(1/v)*(yr/sqrt(zr*zr+xr*xr+yr*yr));
        t0=t-px*xr-py*yr;
        kpx=round((px-px0)/dpx);
        kpy=round((py-py0)/dpy);
        
        if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0)
        datatime0(kpx,kpy)=t0;
        end
        
        kpx=round((-px-px0)/dpx);
        kpy=round((py-py0)/dpy);
        if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0)
        datatime0(kpx,kpy)=t0;
        end
        
        kpx=round((px-px0)/dpx);
        kpy=round((-py-py0)/dpy);
        if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0)
        datatime0(kpx,kpy)=t0;
        end
        
        kpx=round((-px-px0)/dpx);
        kpy=round((-py-py0)/dpy);
        if(kpx<npx && kpy<npy && kpx>0 && kpy>0 && t0>0)
        datatime0(kpx,kpy)=t0;
        end
    end
end

datatime0x=datatime0;
datatime0y=datatime0;
for kpy=1:npy
        t0=0;
    for kpx=1:npx
        if(datatime0y(kpx,kpy)>0)
            t0=datatime0y(kpx,kpy);
        else
            datatime0y(kpx,kpy)=t0;
        end
    end
    for kpx=npx:-1:1
        if(datatime0y(kpx,kpy)>0)
            t0=datatime0y(kpx,kpy);
        else
            datatime0y(kpx,kpy)=t0;
        end
    end
end
for kpx=1:npx
        t0=0;
    for kpy=1:npy
        if(datatime0x(kpx,kpy)>0)
            t0=datatime0x(kpx,kpy);
        else
            datatime0x(kpx,kpy)=t0;
        end
    end
    for kpy=npy:-1:1
        if(datatime0x(kpx,kpy)>0)
            t0=datatime0x(kpx,kpy);
        else
            datatime0x(kpx,kpy)=t0;
        end
    end
end
for kpx=1:npx
    for kpy=1:npy
        datatime0(kpx,kpy)=max(datatime0y(kpx,kpy),datatime0x(kpx,kpy));
        datatime0(kpx,kpy)=(datatime0(kpx,kpy)+0.05)/dt;
    end
end
x=px0:dpx:-px0-dpx;
y=py0:dpy:-py0-dpy;
[x3,y3]=meshgrid(x,y);

figure()
mesh(x3,y3,datatime0')
set(gca,'ZDir','reverse');
set(gca,'ZLim',[0 6000]);
xlabel('Ray parameters X (m/s)')
% ylabel('Ray parameters Y (m/s)')
zlabel('Time sampling')
%调整色标范围以及倒转Z轴
view(-0,1);

