%from2022041401
%看看18步
clear;clc;close all;
tic
img0(:,:,1)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\1.bmp'));
img0(:,:,2)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\2.bmp'));
img0(:,:,3)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\3.bmp'));
img0(:,:,4)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\4.bmp'));
img0(:,:,5)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\5.bmp'));
img0(:,:,6)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\6.bmp'));
img0(:,:,7)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\7.bmp'));
img0(:,:,8)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\8.bmp'));
img0(:,:,9)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\9.bmp'));
img0(:,:,10)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\10.bmp'));
img0(:,:,11)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\11.bmp'));
img0(:,:,12)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\12.bmp'));
img0(:,:,13)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\13.bmp'));
img0(:,:,14)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\14.bmp'));
img0(:,:,15)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\15.bmp'));
img0(:,:,16)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\16.bmp'));
img0(:,:,17)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\17.bmp'));
img0(:,:,18)=double(imread('D:\CLASS\毕业论文\mail\平面 18步相移+3步相移图像\32\18.bmp'));
n0=size(img0);
bs=1024;%均衡化精度

img1=img0;
n=size(img1);

p1=a2ph(img1);  
p1_1=p1(1:10,:);
p1_1=p1_1(:);
clear p1;clear n;

p1_1=sortrows(p1_1);
p1_1_c=pheq(p1_1,bs,min(p1_1(:)),max(p1_1(:)));
%p1_1_c=linspace(min(p1_1),max(p1_1),numel(p1_1));p1_1_c=p1_1_c';%目前看起来还是排序后拉线好使
%*************滤波一下*********************(可跳过此段)
% enl=p1_1-p1_1_c;
% temp=enl(round(numel(enl)/2):end-1);
% enl=cat(1,temp,enl);
% [b,a]=butter(3,2*20/numel(enl),'low');
% enl=filter(b,a,enl(end:-1:1));
% enl=filter(b,a,enl(end:-1:1));
% enl=enl(numel(temp)+1:end);
% p1_1=enl+p1_1_c;
% clear enl;clear temp;clear a;clear b;
%************滤波完成*********************
p1_1=[p1_1;pi;-pi];p1_1_c=[p1_1_c;pi;-pi];
p1_1=sortrows(p1_1);p1_1_c=sortrows(p1_1_c);
%p1_1_c(end)=-pi;

p0=a2ph(img0);
%下面提取有效条纹区域和背景振幅之类的
img0_mean3=mean(img0,3);
img0_mean3_highmean=mean(img0_mean3(img0_mean3>0.5*max(img0_mean3)));
img0_mean3_lowmean=mean(img0_mean3(img0_mean3<0.5*max(img0_mean3)));
img0_mean3_t=((img0_mean3_highmean*2-max(img0_mean3(:)))+(img0_mean3_lowmean*2-min(img0_mean3(:))))/2;
img0_u=img0_mean3>img0_mean3_t;
u=double(img0_u);
u(u==0)=NaN;
clear img0_mean3;clear img0_mean3_t;clear img0_mean3_highmean;clear img0_mean3_lowmean;clear img0_u;
clear img0;
%相位对比转换
t0=(p1_1+pi)/(2*pi)*(bs-1);
t1=0:(bs-1);
[~,x]=min(abs(t1-t0));
tr=p1_1_c(x);
clear x;clear t0;clear t1;
%表建立完成
p2=tr(round((p0+pi)/(2*pi)*(bs-1)+1));
p2=p2.*u;

figure(1);
plot(1:numel(p1_1),sort(p1_1),1:numel(p1_1_c),sort(p1_1_c));
figure(2);
subplot(221);imshow(p0.*u,[]);title('原相位');
subplot(222);imshow(p2,[]);title('校正后相位');
subplot(223);h0=histogram(p0.*u,64);title('校正前直方图');
subplot(224);h1=histogram(p2,64);title('校正后直方图');
figure(3);
plot(p2(10,:));
std(h0.Values/sum(h0.Values))
std(h1.Values/sum(h1.Values))
toc
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


