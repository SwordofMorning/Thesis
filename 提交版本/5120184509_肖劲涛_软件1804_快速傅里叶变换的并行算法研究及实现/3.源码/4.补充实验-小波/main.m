clear;clc;

irPic = imread('./3.jpg');
visPic = imread('./2.jpg');

irPic = rgb2gray(irPic);
visPic = rgb2gray(visPic);

% figure(1);
% contour3(irPic, 200); 
% figure(2);
% contour3(visPic, 200);

%% 小波变换

[c, s] = wavedec2(irPic, 2, 'coif3');
%设置尺度向量
n = [1,2];
%设置阈值向量p
p = [10.12,23.28];
%对三个方向高频系数进行阈值处理
nc = wthcoef2('h',c,s,n,p,'s');
nc = wthcoef2('v',nc,s,n,p,'s');
nc = wthcoef2('d',nc,s,n,p,'s');
%对新的小波分解结构[c,s]进行重构
irPic1 = waverec2(nc,s,'coif3');
irPic1 = uint8(irPic1);

%% 傅里叶变换
F = fft2(irPic);
myangle = angle(F);             %相位谱(没有进行移位的)
FS = abs(fftshift(F));          % 移位，使低频成分集中到图像中心，并得到幅度谱

S = log(1+abs(FS));

% 对幅度图进行操作，去除外围的高频成分的幅度值，也就是将高频成分能量去除了
%（此处只是简单地保留了图像中心 200 X 200 的图像块）
[m,n] = size(FS);
FS(1:m/2-100,:) = 0;
FS(m/2+100:m,:) = 0;
FS(m/2-100:m/2+100,1:n/2-100) = 0;
FS(m/2-100:m/2+100,n/2+100:n) = 0;

SS = log(1+abs(FS));

aaa = ifftshift(FS);          % 将处理后的幅度图反移位，恢复到正常状态
bbb = aaa.*cos(myangle) + aaa.*sin(myangle).*1i;      % 幅度值和相位值重新进行结合，得到复数
fr = abs(ifft2(bbb));               % 进行傅里叶反变换，得到处理后的时域图像
irPic2 = im2uint8(mat2gray(fr));       

%% imshow
figure(3);
subplot(1, 3, 1);
imshow(irPic);
title("原始图像");
subplot(1, 3, 2);
imshow(irPic1);
title("小波");
subplot(1, 3, 3);
imshow(irPic2)
title("傅里叶");