
% outfhead.write((char *)(&gx), sizeof(gx));
% outfhead.write((char *)(&gy), sizeof(gy));
% outfhead.write((char *)(&sx), sizeof(sx));
% outfhead.write((char *)(&sy), sizeof(sy));
% outfhead.write((char *)(&offset), sizeof(offset));
% outfhead.write((char *)(&selev), sizeof(selev));
% outfhead.write((char *)(&gelev), sizeof(gelev));
% outfhead.write((char *)(&sdepth), sizeof(sdepth));
% outfhead.write((char *)(&swdep), sizeof(swdep));
% outfhead.write((char *)(&tracf), sizeof(tracf));

clc
clear
close all
% close all
%% 设置参数
%道头的参数数量必须明确，在不清楚具体道数的情况下可以尽可能取大，以读取所有的数据
ntrace=3994999;
npar=10;
dcdpx=25;
dcdpy=25;

fid=fopen('..\data\3149.cdp.headAllHead','rb');
%读取文件，在到达文件末尾时会自动关闭，读取的data尺寸会取决于读取的数据量
[data,count2]=fread(fid,[npar,ntrace],'float');
fclose(fid);%关闭文件

n=size(data);
n1=npar;
n2=n(2);
datahead2d=single(zeros(n1,n2));
datahead2d=data;
clear data
%计算CMP位置
dataheadadd=single(zeros(2,n2));
dataheadadd(1,:)=(datahead2d(3,:)+datahead2d(1,:))/2;
dataheadadd(2,:)=(datahead2d(4,:)+datahead2d(2,:))/2;
%计算炮点分布范围
minsx=min(datahead2d(3,:));
minsy=min(datahead2d(4,:));
maxsx=max(datahead2d(3,:));
maxsy=max(datahead2d(4,:));
u_1=single(minsx:dcdpx:maxsx);
u_2=single(minsy:dcdpy:maxsy);
[U_1,U_2]=(meshgrid(u_1,u_2));
datacdpnum=single(zeros(size(U_1)));
s=size(dataheadadd);
%统计CMP覆盖次数
for k=1:s(2)
    k2=round((dataheadadd(1,k)-minsx)/dcdpx);
    k1=round((dataheadadd(2,k)-minsy)/dcdpy);
    datacdpnum(k1,k2)=datacdpnum(k1,k2)+1;
end
%绘图，确定图的宽高和位置，各个图层参数被保存在数组 draw 中；
figure('Units','centimeter','Position',[5 5 30 15])
draw(1)=scatter(datahead2d(3,:),datahead2d(4,:),15,'filled','b');
hold on
draw(2)=scatter(datahead2d(1,:),datahead2d(2,:),20,'filled','r');
hold on
%为了设置不透明的图例单独绘制了一个图层
draw(3)=scatter(dataheadadd(1,1),dataheadadd(2,1),0.1,'filled','g');
draw(3).MarkerFaceAlpha = 0.5; %设置散点的透明度，该属性也会表现在图例中
hold on
draw(4)=scatter(dataheadadd(1,:),dataheadadd(2,:),15,'filled','g');
draw(4).MarkerFaceAlpha = 0.02; %设置散点的透明度
hold on
%设置图例
legend(draw(1:3),'Source','OBN','CMP')
hold on
xlabel('X (m)')
ylabel('Y (m)')
%设置X轴和Y轴范围
axis([minsx-10*dcdpx maxsx+10*dcdpx minsy-10*dcdpy maxsy+10*dcdpy])
%倒转Z轴
% set(gca,'ZDir','reverse');
hold on
%显示网格线
grid on
%显示坐标轴
axis on
%显示画框
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','centimeter','Position',[5 5 30 15])
pcolor(U_1,U_2,datacdpnum)
hold on
shading interp;
axis([minsx-10*dcdpx maxsx+10*dcdpx minsy-10*dcdpy maxsy+10*dcdpy])
view(0,90)
xlabel('X (m)')
ylabel('Y (m)')
hold on
grid on
axis on
box on
