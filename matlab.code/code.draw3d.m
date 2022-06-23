
clc
clear
close all
%% 设置参数
%三维数据采样点数
nline=75;
ntrace=75;
nt=2048;
n1=nline;
n2=ntrace;
n3=nt;

%三维数据采样间隔
dline=0.00001;
dtrace=0.00001;
dt=0.00010;
d1=dline;
d2=dtrace;
d3=dt;

%三维数据采样坐标范围
line0=-0.00037;
trace0=-0.00037;
t0=0;
n1min=line0;
n2min=trace0;
n3min=t0;
n1max=n1min+d1*(n1-1);
n2max=n2min+d2*(n2-1);
n3max=n3min+d3*(n3-1);

%设置缺失数据块角点的三维坐标点，即凹陷角坐标
%此三维坐标点对应采样点数而非采样间隔，表示各个维度上缺失的数据点数
pline=37;
ptrace=60;
pt=1500;
p1=pline;
p2=ptrace;
p3=nt-pt;

%% 读取数据
data3d=zeros(n2,n1,n3);
data2d=zeros(n3,n2);
data2d1=zeros(n2,n3);
data2d2=zeros(n2,n3);
%打开二进制文件路径
fid=fopen('datatp.obn.bin','rb');
    if(fid>0)
        for k=1:n1
            %读取数据的维度，列方向维度为nt，行方向维度为道数ntrace
            %一个二维数据为一条line
            [data2d,count1]=fread(fid,[n3,n2],'float');
            data2d1=data2d';
            for k2=1:n3
                %倒转数据Z轴，使得时间/深度从小到大显示
                data2d2(:,n3+1-k2)=data2d1(:,k2);
            end
            data3d(:,k,:)=data2d2;
        end
    end   
fclose(fid);%关闭文件
data2d=zeros(1,1);
data2d1=zeros(1,1);
data2d2=zeros(1,1);

%% 定义坐标网格
%分别定义三个维度上的规则网格坐标
u_1=n1min:d1:n1max;
u_2=n2min:d2:n2max;
u_3=n3max:-d3:n3min;
%定义三维的规则网格坐标，每个矩阵表示一个维度的坐标值

%% 拷贝数据
Z=data3d;
data3d=zeros(1,1,1);

%% 绘图：缺角图
%绘图时需要缺失一个角数据块，因而分三块绘制数据
%获取数据块1及其坐标
d1=Z(p2:end,:,:);
u_11=u_1;
u_21=u_2(p2:end);
u_31=u_3;
[U_11,U_21,U_31]=meshgrid(u_11,u_21,u_31);

%获取数据块2及其坐标
d2=Z(1:p2,:,1:p3);
u_12=u_1;
u_22=u_2(1:p2);
u_32=u_3(1:p3);
[U_12,U_22,U_32]=meshgrid(u_12,u_22,u_32);

%获取数据块3及其坐标
d3=Z(1:p2,p1:end,p3:end);
u_13=u_1(p1:end);
u_23=u_2(1:p2);
u_33=u_3(p3:end);
[U_13,U_23,U_33]=meshgrid(u_13,u_23,u_33);

figure(2)
hold on
%绘制数据块1
slice(U_11,U_21,U_31,d1,u_11,u_21,u_31);
U_11=zeros(1,1,1);
U_21=zeros(1,1,1);
U_31=zeros(1,1,1);
hold on
%绘制数据块1的边框线；需要框线的三维坐标
plot3(ones(1,length(u_21))*u_11(1),u_21,ones(1,length(u_21))*u_31(end),'k-')
hold on
plot3(u_11(1:p1),ones(1,length(u_11(1:p1)))*u_21(1),ones(1,length(u_11(1:p1)))*u_31(end),'k-')
hold on
plot3(ones(1,length(u_31(p3:end)))*u_11(1),ones(1,length(u_31(p3:end)))*u_21(1),u_31(p3:end),'k-')

%绘制数据块2
hold on
slice(U_12,U_22,U_32,d2,u_12,u_22,u_32);
U_12=zeros(1,1,1);
U_22=zeros(1,1,1);
U_32=zeros(1,1,1);
hold on
%绘制数据块2的边框线；需要框线的三维坐标
plot3(ones(1,length(u_22))*u_12(1),u_22,ones(1,length(u_22))*u_32(end),'k-')
hold on
plot3(u_12(1:p1),ones(1,length(u_12(1:p1)))*u_22(1),ones(1,length(u_12(1:p1)))*u_32(end),'k-')
hold on
plot3(ones(1,length(u_32))*u_12(1),ones(1,length(u_32))*u_22(1),u_32,'k-')

%绘制数据块3
hold on
slice(U_13,U_23,U_33,d3,u_13,u_23,u_33);
U_13=zeros(1,1,1);
U_23=zeros(1,1,1);
U_33=zeros(1,1,1);
hold on
%绘制数据块3的边框线；需要框线的三维坐标
plot3(ones(1,length(u_23))*u_13(1),u_23,ones(1,length(u_23))*u_33(end),'k-')
hold on
plot3(u_13,ones(1,length(u_13))*u_23(1),ones(1,length(u_13))*u_33(end),'k-')
hold on
plot3(ones(1,length(u_33))*u_13(1),ones(1,length(u_33))*u_23(1),u_33,'k-')

hold on
%绘制缺失数据块的边框线；需要框线的三维坐标
plot3(ones(1,length(u_2))*u_1(p1),u_2,ones(1,length(u_2))*u_3(p3),'k--')
hold on
plot3(u_1,ones(1,length(u_1))*u_2(p2),ones(1,length(u_1))*u_3(p3),'k--')
hold on
plot3(ones(1,length(u_3))*u_1(p1),ones(1,length(u_3))*u_2(p2),u_3,'k--')

xlabel('u_1')
ylabel('u_2')
zlabel('u_3')
%调整色标范围以及倒转Z轴
caxis([-0.00005 0.00005])
set(gca,'ZDir','reverse');
axis([n1min n1max n2min n2max n3min n3max])

%设置三维观察视角
view(-60,20);
shading flat
colorbar;

%% 绘图：三维图
% Z=sin(U_1)+cos(U_2)+exp(U_3);
[U_1,U_2,U_3]=meshgrid(u_1,u_2,u_3);
figure(1)
slice(U_1,U_2,U_3,Z,u_1,u_2,u_3);
U_1=zeros(1,1,1);
U_2=zeros(1,1,1);
U_3=zeros(1,1,1);
hold on
plot3(ones(1,length(u_2))*u_1(1),u_2,ones(1,length(u_2))*u_3(end),'k-')
hold on
plot3(u_1,ones(1,length(u_1))*u_2(1),ones(1,length(u_1))*u_3(end),'k-')
hold on
plot3(ones(1,length(u_3))*u_1(1),ones(1,length(u_3))*u_2(1),u_3,'k-')
xlabel('u_1')
ylabel('u_2')
zlabel('u_3')
set(gca,'ZDir','reverse');
% axis([0 0.3 0 0.4 0 0.5])
shading flat
colorbar;






