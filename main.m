clear;clc;close all;
tic
% 就三步相移吧 这里是bmp，就单通道
img0(:,:,1)=double(imread('1.bmp'));
img0(:,:,2)=double(imread('2.bmp'));
img0(:,:,3)=double(imread('3.bmp'));
n0=size(img0);
bs=2^10;%均衡化精度

img1=img0(round(0.05*n0(1)):round(0.95*n0(1)),round(0.05*n0(2)):round(0.95*n0(2)),:); %xjb裁剪掉不要的地方
%n=size(img1);
p1=a2ph(img1); 
%转为截断相位

p1_1=p1';
p1_1=p1_1(1:5120);

p1_1=p1_1(:);
clear p1,clear n;

p1_1_c=pheq(p1_1,bs,min(p1_1(:)),max(p1_1(:)));
% 可以排序然后滤波降噪，但还是算了
p1_1_c=p1_1_c-p1_1;

p0=a2ph(img0);
%下面提取有效条纹区域和背景振幅之类的 xjb写 能用就行
% img0_mean3=mean(img0,3);
% img0_mean3_highmean=mean(img0_mean3(img0_mean3>0.5*max(img0_mean3)));
% img0_mean3_lowmean=mean(img0_mean3(img0_mean3<0.5*max(img0_mean3)));
% img0_mean3_t=((img0_mean3_highmean*2-max(img0_mean3(:)))+(img0_mean3_lowmean*2-min(img0_mean3(:))))/2;
% img0_u=img0_mean3>img0_mean3_t;
% u=double(img0_u);
% u(u==0)=NaN;
% clear img0_mean3;clear img0_mean3_t;clear img0_mean3_highmean;clear img0_mean3_lowmean;clear img0_u;
u=1; % 不提取就取1吧
%相位对比转换
t0=(p1_1+pi)/(2*pi)*(bs-1);
t1=0:(bs-1);
[~,x]=min(abs(t1-t0));
tr=p1_1_c(x);
clear x;clear t0;clear t1;
%表建立完成
p2=tr(round((p0+pi)/(2*pi)*(bs-1)+1));
p2=(p2+p0).*u;

p0=p0.*u;
% figure(1);
% plot(1:numel(p1_1),sort(p1_1),1:numel(p1_1_c),sort(p1_1_c));
% figure(2);
% subplot(221);imshow(p0,[]);title('原相位');
% subplot(222);imshow(p2,[]);title('校正后相位');
% subplot(223);h0=histogram(p0,128);title('校正前直方图');
% subplot(224);h1=histogram(p2,128);title('校正后直方图');
% figure(3);
% plot(p2(10,:));
% std(h0.Values/sum(h0.Values))*100
% std(h1.Values/sum(h1.Values))
toc

temprms

%函数
function I=a2ph(img)
%转相位
n=size(img);
a=zeros(n(1),n(2));
b=a;
for k=0:n(3)-1
    a=a+img(:,:,k+1)*sin(2*pi*k/n(3));
    b=b+img(:,:,k+1)*cos(2*pi*k/n(3));
end
I=-atan2(a,b);
end
function img1=pheq(img0,n,r0,r1)
%其实只是均衡化,感觉不如用matalb自带的
%但至少输入输出的范围没那么麻烦
%对NaN也适用改造
img0=double(img0);
img1=zeros(size(img0));
img0=(img0-min(img0(:)))./(max(img0(:))-min(img0(:)));
img0=round(img0*(n-1));
num=sum(~isnan(img0(:)));
for k=0:n-1
    img1=img1+(img0==k)*sum(img0(:)<=k)/num;
end
img1(isnan(img0))=NaN;
img1=(img1-min(img1(:)))./(max(img1(:))-min(img1(:)));
img1=img1*(r1-r0)+r0;
end


