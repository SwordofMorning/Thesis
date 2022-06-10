clear;clc;

irPic = imread('./3.jpg');
visPic = imread('./2.jpg');

irPic = rgb2gray(irPic);
visPic = rgb2gray(visPic);

% figure(1);
% contour3(irPic, 200); 
% figure(2);
% contour3(visPic, 200);

%% С���任

[c, s] = wavedec2(irPic, 2, 'coif3');
%���ó߶�����
n = [1,2];
%������ֵ����p
p = [10.12,23.28];
%�����������Ƶϵ��������ֵ����
nc = wthcoef2('h',c,s,n,p,'s');
nc = wthcoef2('v',nc,s,n,p,'s');
nc = wthcoef2('d',nc,s,n,p,'s');
%���µ�С���ֽ�ṹ[c,s]�����ع�
irPic1 = waverec2(nc,s,'coif3');
irPic1 = uint8(irPic1);

%% ����Ҷ�任
F = fft2(irPic);
myangle = angle(F);             %��λ��(û�н�����λ��)
FS = abs(fftshift(F));          % ��λ��ʹ��Ƶ�ɷּ��е�ͼ�����ģ����õ�������

S = log(1+abs(FS));

% �Է���ͼ���в�����ȥ����Χ�ĸ�Ƶ�ɷֵķ���ֵ��Ҳ���ǽ���Ƶ�ɷ�����ȥ����
%���˴�ֻ�Ǽ򵥵ر�����ͼ������ 200 X 200 ��ͼ��飩
[m,n] = size(FS);
FS(1:m/2-100,:) = 0;
FS(m/2+100:m,:) = 0;
FS(m/2-100:m/2+100,1:n/2-100) = 0;
FS(m/2-100:m/2+100,n/2+100:n) = 0;

SS = log(1+abs(FS));

aaa = ifftshift(FS);          % �������ķ���ͼ����λ���ָ�������״̬
bbb = aaa.*cos(myangle) + aaa.*sin(myangle).*1i;      % ����ֵ����λֵ���½��н�ϣ��õ�����
fr = abs(ifft2(bbb));               % ���и���Ҷ���任���õ�������ʱ��ͼ��
irPic2 = im2uint8(mat2gray(fr));       

%% imshow
figure(3);
subplot(1, 3, 1);
imshow(irPic);
title("ԭʼͼ��");
subplot(1, 3, 2);
imshow(irPic1);
title("С��");
subplot(1, 3, 3);
imshow(irPic2)
title("����Ҷ");