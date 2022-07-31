
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
%% ���ò���
%��ͷ�Ĳ�������������ȷ���ڲ�����������������¿��Ծ�����ȡ���Զ�ȡ���е�����
ntrace=3994999;
npar=10;
dcdpx=25;
dcdpy=25;

fid=fopen('..\data\3149.cdp.headAllHead','rb');
%��ȡ�ļ����ڵ����ļ�ĩβʱ���Զ��رգ���ȡ��data�ߴ��ȡ���ڶ�ȡ��������
[data,count2]=fread(fid,[npar,ntrace],'float');
fclose(fid);%�ر��ļ�

n=size(data);
n1=npar;
n2=n(2);
datahead2d=single(zeros(n1,n2));
datahead2d=data;
clear data
%����CMPλ��
dataheadadd=single(zeros(2,n2));
dataheadadd(1,:)=(datahead2d(3,:)+datahead2d(1,:))/2;
dataheadadd(2,:)=(datahead2d(4,:)+datahead2d(2,:))/2;
%�����ڵ�ֲ���Χ
minsx=min(datahead2d(3,:));
minsy=min(datahead2d(4,:));
maxsx=max(datahead2d(3,:));
maxsy=max(datahead2d(4,:));
u_1=single(minsx:dcdpx:maxsx);
u_2=single(minsy:dcdpy:maxsy);
[U_1,U_2]=(meshgrid(u_1,u_2));
datacdpnum=single(zeros(size(U_1)));
s=size(dataheadadd);
%ͳ��CMP���Ǵ���
for k=1:s(2)
    k2=round((dataheadadd(1,k)-minsx)/dcdpx);
    k1=round((dataheadadd(2,k)-minsy)/dcdpy);
    datacdpnum(k1,k2)=datacdpnum(k1,k2)+1;
end
%��ͼ��ȷ��ͼ�Ŀ�ߺ�λ�ã�����ͼ����������������� draw �У�
figure('Units','centimeter','Position',[5 5 30 15])
draw(1)=scatter(datahead2d(3,:),datahead2d(4,:),15,'filled','b');
hold on
draw(2)=scatter(datahead2d(1,:),datahead2d(2,:),20,'filled','r');
hold on
%Ϊ�����ò�͸����ͼ������������һ��ͼ��
draw(3)=scatter(dataheadadd(1,1),dataheadadd(2,1),0.1,'filled','g');
draw(3).MarkerFaceAlpha = 0.5; %����ɢ���͸���ȣ�������Ҳ�������ͼ����
hold on
draw(4)=scatter(dataheadadd(1,:),dataheadadd(2,:),15,'filled','g');
draw(4).MarkerFaceAlpha = 0.02; %����ɢ���͸����
hold on
%����ͼ��
legend(draw(1:3),'Source','OBN','CMP')
hold on
xlabel('X (m)')
ylabel('Y (m)')
%����X���Y�᷶Χ
axis([minsx-10*dcdpx maxsx+10*dcdpx minsy-10*dcdpy maxsy+10*dcdpy])
%��תZ��
% set(gca,'ZDir','reverse');
hold on
%��ʾ������
grid on
%��ʾ������
axis on
%��ʾ����
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
