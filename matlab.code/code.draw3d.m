
clc
clear
close all
%% ���ò���
%��ά���ݲ�������
nline=75;
ntrace=75;
nt=2048;
n1=nline;
n2=ntrace;
n3=nt;

%��ά���ݲ������
dline=0.00001;
dtrace=0.00001;
dt=0.00010;
d1=dline;
d2=dtrace;
d3=dt;

%��ά���ݲ������귶Χ
line0=-0.00037;
trace0=-0.00037;
t0=0;
n1min=line0;
n2min=trace0;
n3min=t0;
n1max=n1min+d1*(n1-1);
n2max=n2min+d2*(n2-1);
n3max=n3min+d3*(n3-1);

%����ȱʧ���ݿ�ǵ����ά����㣬�����ݽ�����
%����ά������Ӧ�����������ǲ����������ʾ����ά����ȱʧ�����ݵ���
pline=37;
ptrace=60;
pt=1500;
p1=pline;
p2=ptrace;
p3=nt-pt;

%% ��ȡ����
data3d=zeros(n2,n1,n3);
data2d=zeros(n3,n2);
data2d1=zeros(n2,n3);
data2d2=zeros(n2,n3);
%�򿪶������ļ�·��
fid=fopen('datatp.obn.bin','rb');
    if(fid>0)
        for k=1:n1
            %��ȡ���ݵ�ά�ȣ��з���ά��Ϊnt���з���ά��Ϊ����ntrace
            %һ����ά����Ϊһ��line
            [data2d,count1]=fread(fid,[n3,n2],'float');
            data2d1=data2d';
            for k2=1:n3
                %��ת����Z�ᣬʹ��ʱ��/��ȴ�С������ʾ
                data2d2(:,n3+1-k2)=data2d1(:,k2);
            end
            data3d(:,k,:)=data2d2;
        end
    end   
fclose(fid);%�ر��ļ�
data2d=zeros(1,1);
data2d1=zeros(1,1);
data2d2=zeros(1,1);

%% ������������
%�ֱ�������ά���ϵĹ�����������
u_1=n1min:d1:n1max;
u_2=n2min:d2:n2max;
u_3=n3max:-d3:n3min;
%������ά�Ĺ����������꣬ÿ�������ʾһ��ά�ȵ�����ֵ

%% ��������
Z=data3d;
data3d=zeros(1,1,1);

%% ��ͼ��ȱ��ͼ
%��ͼʱ��Ҫȱʧһ�������ݿ飬����������������
%��ȡ���ݿ�1��������
d1=Z(p2:end,:,:);
u_11=u_1;
u_21=u_2(p2:end);
u_31=u_3;
[U_11,U_21,U_31]=meshgrid(u_11,u_21,u_31);

%��ȡ���ݿ�2��������
d2=Z(1:p2,:,1:p3);
u_12=u_1;
u_22=u_2(1:p2);
u_32=u_3(1:p3);
[U_12,U_22,U_32]=meshgrid(u_12,u_22,u_32);

%��ȡ���ݿ�3��������
d3=Z(1:p2,p1:end,p3:end);
u_13=u_1(p1:end);
u_23=u_2(1:p2);
u_33=u_3(p3:end);
[U_13,U_23,U_33]=meshgrid(u_13,u_23,u_33);

figure(2)
hold on
%�������ݿ�1
slice(U_11,U_21,U_31,d1,u_11,u_21,u_31);
U_11=zeros(1,1,1);
U_21=zeros(1,1,1);
U_31=zeros(1,1,1);
hold on
%�������ݿ�1�ı߿��ߣ���Ҫ���ߵ���ά����
plot3(ones(1,length(u_21))*u_11(1),u_21,ones(1,length(u_21))*u_31(end),'k-')
hold on
plot3(u_11(1:p1),ones(1,length(u_11(1:p1)))*u_21(1),ones(1,length(u_11(1:p1)))*u_31(end),'k-')
hold on
plot3(ones(1,length(u_31(p3:end)))*u_11(1),ones(1,length(u_31(p3:end)))*u_21(1),u_31(p3:end),'k-')

%�������ݿ�2
hold on
slice(U_12,U_22,U_32,d2,u_12,u_22,u_32);
U_12=zeros(1,1,1);
U_22=zeros(1,1,1);
U_32=zeros(1,1,1);
hold on
%�������ݿ�2�ı߿��ߣ���Ҫ���ߵ���ά����
plot3(ones(1,length(u_22))*u_12(1),u_22,ones(1,length(u_22))*u_32(end),'k-')
hold on
plot3(u_12(1:p1),ones(1,length(u_12(1:p1)))*u_22(1),ones(1,length(u_12(1:p1)))*u_32(end),'k-')
hold on
plot3(ones(1,length(u_32))*u_12(1),ones(1,length(u_32))*u_22(1),u_32,'k-')

%�������ݿ�3
hold on
slice(U_13,U_23,U_33,d3,u_13,u_23,u_33);
U_13=zeros(1,1,1);
U_23=zeros(1,1,1);
U_33=zeros(1,1,1);
hold on
%�������ݿ�3�ı߿��ߣ���Ҫ���ߵ���ά����
plot3(ones(1,length(u_23))*u_13(1),u_23,ones(1,length(u_23))*u_33(end),'k-')
hold on
plot3(u_13,ones(1,length(u_13))*u_23(1),ones(1,length(u_13))*u_33(end),'k-')
hold on
plot3(ones(1,length(u_33))*u_13(1),ones(1,length(u_33))*u_23(1),u_33,'k-')

hold on
%����ȱʧ���ݿ�ı߿��ߣ���Ҫ���ߵ���ά����
plot3(ones(1,length(u_2))*u_1(p1),u_2,ones(1,length(u_2))*u_3(p3),'k--')
hold on
plot3(u_1,ones(1,length(u_1))*u_2(p2),ones(1,length(u_1))*u_3(p3),'k--')
hold on
plot3(ones(1,length(u_3))*u_1(p1),ones(1,length(u_3))*u_2(p2),u_3,'k--')

xlabel('u_1')
ylabel('u_2')
zlabel('u_3')
%����ɫ�귶Χ�Լ���תZ��
caxis([-0.00005 0.00005])
set(gca,'ZDir','reverse');
axis([n1min n1max n2min n2max n3min n3max])

%������ά�۲��ӽ�
view(-60,20);
shading flat
colorbar;

%% ��ͼ����άͼ
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






