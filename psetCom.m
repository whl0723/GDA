clc
clear
close all

x=linspace(-1.5,1.5,300);
y=linspace(-1.5,1.5,300);                                                  %像素300*300 每个像素是10um*10um
[x,y]=meshgrid(x,y);
A=100*exp(-0.05*(x.^2+y.^2));  %背景项
B=80*exp(-0.01*(x.^2+y.^2));   %调制项
% A=100*ones(size(x));  %背景项
% B=100*ones(size(x));   %调制项

p0=peaks(300) ;          %预设相位 
% p0=p0-min(p0(:));
% p=p0./max(p0(:))*0.6*pi;
% p=p-p(1,1);
% p0=255-double(rgb2gray(imread('cell.jpg')));
% p0(p0>115)=115;
p0=p0./max(p0(:))*pi*0.8;
p = imresize(p0,size(A),'bilinear');

Ib(:,:,1)=A+B.*cos(p);
Ib(:,:,2)=A+B.*cos(p+2);
Ib(:,:,3)=A+B.*cos(p+3);

noise=40;
Ib(:,:,1)=awgn(Ib(:,:,1),noise,'measured');
Ib(:,:,2)=awgn(Ib(:,:,2),noise,'measured');
Ib(:,:,3)=awgn(Ib(:,:,3),noise,'measured');

I1=Ib(:,:,1)-Ib(:,:,2);
I2=Ib(:,:,1)-Ib(:,:,3);
pw1=GS1(I1,I2);
pw2=rand(size(I1));
pw3=I1;
[pw11,RR1,alph]=GD_GS(I1,I2,pw1);
[pw21,RR2,alph]=GD_GS(I1,I2,pw2);
[pw31,RR3,alph]=GD_GS(I1,I2,pw3);

p1=unwrap2(pw11);
p2=unwrap2(pw21);
p3=unwrap2(pw31);

e1=p1-p;
e2=p2-p;
e3=p3-p;

RMSE_1=std2(e1)
RMSE_2=std2(e2)
RMSE_3=std2(e3)

figure 
subplot 121
hold on
plot(RR1(2:end),'r*')
plot(RR2(2:end),'b--')
plot(RR3(2:end),'ko')
legend('GS','Random','Interferogram')
xlabel('Iterations number')
title('Phase range = 0.8 \pi')

p0=peaks(300) ;          %预设相位 
% p0=p0-min(p0(:));
% p=p0./max(p0(:))*0.6*pi;
% p=p-p(1,1);
% p0=255-double(rgb2gray(imread('cell.jpg')));
% p0(p0>115)=115;
p0=p0./max(p0(:))*pi*4;
p = imresize(p0,size(A),'bilinear');

Ib(:,:,1)=A+B.*cos(p);
Ib(:,:,2)=A+B.*cos(p+2);
Ib(:,:,3)=A+B.*cos(p+3);

noise=4000;
Ib(:,:,1)=awgn(Ib(:,:,1),noise,'measured');
Ib(:,:,2)=awgn(Ib(:,:,2),noise,'measured');
Ib(:,:,3)=awgn(Ib(:,:,3),noise,'measured');

I1=Ib(:,:,1)-Ib(:,:,2);
I2=Ib(:,:,1)-Ib(:,:,3);
pw1=GS1(I1,I2);
pw2=rand(size(I1));
pw3=I1;
[pw11,RR1,alph]=GD_GS(I1,I2,pw1);
[pw21,RR2,alph]=GD_GS(I1,I2,pw2);
[pw31,RR3,alph]=GD_GS(I1,I2,pw3);

p1=unwrap2(pw11);
p2=unwrap2(pw21);
p3=unwrap2(pw31);

e1=p1-p;
e2=p2-p;
e3=p3-p;

RMSE_1=std2(e1)
RMSE_2=std2(e2)
RMSE_3=std2(e3)


subplot 122
hold on
plot(RR1(2:end),'r*')
plot(RR2(2:end),'b--')
plot(RR3(2:end),'ko')
legend('GS','Random','Interferogram')
xlabel('Iterations number ')
title('Phase range = 4 \pi')


% figure
% subplot 121
% hold on
% plot(RR1(2:end),'r*')
% plot(RR2(2:end),'b--')
% legend('GS','Random')
% xlabel('Number of iterations')
% title('Phase range = 0.8 \pi')
% 
% subplot 122
% hold on
% plot(RR1(2:end),'r*')
% plot(RR2(2:end),'b--')
% legend('GS','Random')
% xlabel('Number of iterations ')
% title('Phase range = 4 \pi')