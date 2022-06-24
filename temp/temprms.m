%求均方根误差的脚本
% close all;
load('matlab1.mat');
temp=abs(p000-p2);
temp=abs(temp-(temp>pi)*2*pi);
histogram(temp,64)
p20=temp;
temp=temp.^2;
temp=mean(temp(:));
temp=sqrt(temp)
% close all;