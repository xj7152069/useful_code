% 
% outfhead.write((char *)(&gx), sizeof(gx));
% outfhead.write((char *)(&gy), sizeof(gy));
% outfhead.write((char *)(&sx), sizeof(sx));
% outfhead.write((char *)(&sy), sizeof(sy));
% outfhead.write((char *)(&offset), sizeof(offset));
% outfhead.write((char *)(&selev), sizeof(selev));
% outfhead.write((char *)(&gelev), sizeof(gelev));
% outfhead.write((char *)(&sdepth), sizeof(sdepth));
% outfhead.write((char *)(&swdep), sizeof(swdep));

clc
clear
% close all
%% 设置参数
    name = {'.\data\2685.head','.\data\2701.head','.\data\2717.head',...
        '.\data\2733.head','.\data\2749.head','.\data\2765.head','.\data\2781.head',...
        '.\data\2797.head','.\data\2813.head','.\data\2821.head','.\data\2829.head',...
        '.\data\2837.head','.\data\2845.head','.\data\2853.head','.\data\2861.head',...
        '.\data\2877.head','.\data\2885.head','.\data\2893.head','.\data\2925.head','.\data\3149.head'};
    obnline = {'RL2685','RL2701','RL2717',...
        'RL2733','RL2749','RL2765','RL2781',...
        'RL2797','RL2813','RL2821','RL2829',...
        'RL2837','RL2845','RL2853','RL2861',...
        'RL2877','RL2885','RL2893','RL2925','RL3149'};
    signalntrace=[20672, 20672, 20672, 20672,...
        20672, 20672, 20672, 20672,...
        20672, 51680, 20672, 51680,...
        20672, 51680, 20672, 20672,...
        51680, 20672, 20672, 20672];
    figure()
        %% 读取数据
    npar=9;
    data2d=zeros(1,1);
    data2dswdepth=zeros(npar,1);
    kk1=1;
    kk2=20;
for kk=kk1:kk2
kk
    %三维数据采样点数

    ntrace=signalntrace(kk)*4;
    ncomp=4;
    n1=ntrace;
    n2=npar;

    data2dunit=zeros(npar,ntrace);
    %打开二进制文件路径
    fid=fopen(name{kk},'r');
    if fid==-1
        disp('can not find file!');
    end
    js=0;
    while ~feof(fid)
        js=js+1;
        [data2dunit,count1]=fread(fid,[npar,ntrace],'float');
        if js==1
            data2d=data2dunit;
             data2dswdepth=[data2dswdepth data2dunit];
        elseif size(data2dunit,2)>0
            data2d=[data2d data2dunit(:,1)];
        end
    end  
    fclose(fid);%关闭文件
%     [datasort,site]=sort(data2d,2);
%     for s=1:size(site,2)
%         datasort(:,s)=data2d(:,site(10,s));
%     end
%     data2d=datasort;
% figure()
%     sou(kk)=scatter(data2d(3,1:ntrace),data2d(4,1:ntrace),5);
%     hold on
%     plot(data2d(3,1:ntrace),data2d(4,1:ntrace));
%     hold on
    obn(kk)=scatter(data2d(1,ntrace-1:end),data2d(2,ntrace-1:end),20,'filled');
    hold on
    legend(name{kk})
    hold on
    xlabel('SX')
    ylabel('SY')
    axis([1470000 1650000 20565000 20660000])
    hold on
end
%     legend(sou(kk1:kk2),obnline{kk1:kk2})
    legend(obn(kk1:kk2),obnline{kk1:kk2})
    hold on
title('DF11')
% 
[xq, yq] = meshgrid(min(data2dswdepth(3,2:end)):50:max(data2dswdepth(3,2:end)),...
    min(data2dswdepth(4,2:end)):50:max(data2dswdepth(4,2:end)));
zq = griddata(data2dswdepth(3,2:end),data2dswdepth(4,2:end),data2dswdepth(9,2:end),xq,yq);
% [XI,YI,ZI] = griddata(.......,method)
% 用指定的算法method 计算：
% ‘linear’：基于三角形的线性插值（缺省算法）；
% ‘cubic’： 基于三角形的三次插值；
% ‘nearest’：最邻近插值法；
% ‘v4’：MATLAB 4 中的griddata 算法。
zq=zq/10.0;
xq=xq/10.0;
yq=yq/10.0;
figure()
% plot3(data2dsx(2:end,1),data2dsy(2:end,1),data2dseabase)
surf(xq, yq, zq)
xlabel('SX (m)')
ylabel('SY (m)')
zlabel('Water Depth (m)')
hold on
title('Seabed Depth')
set(gca,'ZDir','reverse');
axis square
shading interp
zlim([0 100])

zq2=zq(22:1810,3:3222);
xq2=xq(22:1810,3:3222);
yq2=yq(22:1810,3:3222);
    fid=fopen('.\data\df11.seafloor.depth.bin','w');
    count=fwrite(fid,zq2,'float');
    fclose(fid);%关闭文件
    fid=fopen('.\data\df11.seafloor.sx.bin','w');
    count=fwrite(fid,xq2,'float');
    fclose(fid);%关闭文件
    fid=fopen('.\data\df11.seafloor.sy.bin','w');
    count=fwrite(fid,yq2,'float');
    fclose(fid);%关闭文件
% figure()
% nx=1600;
% plot(yq(:,nx)/10, zq(:,nx)/10,'LineWidth',2)
% xq(1,nx)
% hold on
% nx=850;
% plot(yq(:,nx)/10, zq(:,nx)/10,'LineWidth',2)
% xq(1,nx)
% hold on
% nx=100;
% plot(yq(:,nx)/10, zq(:,nx)/10,'LineWidth',2)
% xq(1,nx)
% hold on
% xlabel('SY (m)')
% ylabel('Water Depth (m)')
% ylim([00 100])
% set(gca,'YDir','reverse');
% set(gcf, 'position',[100, 100, 700, 300]);
% 
% figure()
% ny=900;
% plot(xq(ny,:)/10, zq(ny,:)/10,'LineWidth',2)
% yq(ny,1)
% hold on
% ny=500;
% plot(xq(ny,:)/10, zq(ny,:)/10,'LineWidth',2)
% yq(ny,1)
% hold on
% ny=100;
% plot(xq(ny,:)/10, zq(ny,:)/10,'LineWidth',2)
% yq(ny,1)
% hold on
% xlabel('SX (m)')
% ylabel('Water Depth (m)')
% ylim([0 100])
% set(gca,'YDir','reverse');
