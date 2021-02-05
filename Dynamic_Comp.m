clc
clear
close all

x=linspace(-1.5,1.5,300);
y=linspace(-1.5,1.5,300);                                                  %像素300*300 每个像素是10um*10um
[x,y]=meshgrid(x,y);
A=100*exp(-0.05*(x.^2+y.^2));  %背景项
B=80*exp(-0.01*(x.^2+y.^2));   %调制项

% p0=peaks(300) ;          %预设相位 
% p0=p0-min(p0(:));
% p=p0./max(p0(:))*0.4*pi;
% p=p-p(1,1);
p0=255-double(rgb2gray(imread('cell.jpg')));
p0(p0>115)=115;
p0=p0./max(p0(:))*pi*0.8;
p = imresize(p0,size(A),'bilinear');

Ib(:,:,1)=A+B.*cos(p);
Ib(:,:,2)=A+B.*cos(p+2);
pss=2.1:0.2:4;
% psr=(pss-1)./2;
noise=40;
Ib(:,:,1)=awgn(Ib(:,:,1),noise,'measured');
Ib(:,:,2)=awgn(Ib(:,:,2),noise,'measured');
ps=linspace(0,4,3);MaxIter=500;
fl=1.5;
fh=60;
for i=1:length(pss)
Ib(:,:,3)=A+B.*cos(p+pss(i));
Ib(:,:,3)=awgn(Ib(:,:,3),noise,'measured');

%%%%%%%%%%%%
% figure
% subplot 131
% imshow(Ib(:,:,1),[])
% subplot 132
% imshow(Ib(:,:,2),[])
% subplot 133
% imshow(Ib(:,:,3),[])

%%%%%%%%%%%%%%%% 高斯滤波
% sigma = 3;
% gausFilter = fspecial('gaussian', [5,5], sigma);
% Ib(:,:,1)= imfilter(Ib(:,:,1), gausFilter, 'replicate');
% Ib(:,:,2)= imfilter(Ib(:,:,2), gausFilter, 'replicate');
% Ib(:,:,3)= imfilter(Ib(:,:,3), gausFilter, 'replicate');

%%%%%%%%%%%%% IGS GS
I1=Ib(:,:,1)-Ib(:,:,2);
I2=Ib(:,:,1)-Ib(:,:,3);
pwGS=GS1(I1,I2);
% % tic
[pwIGS,RR,alph]=GD_GS(I1,I2,pwGS);
% tigs=toc

%%%%%%%%%%%%%%%% AIA

% tic
[pwAIA,delta,iter] = S_AIA(Ib,ps,MaxIter);
% taia=toc

%%%%%%%%%%%%%%%%%%% MSSM

% tic
pwMSSM=wy3(Ib(:,:,1),Ib(:,:,2),Ib(:,:,3),fl,fh);
% tmssm=toc

%%%%%%%%%%%%%%%%%%%
pwLSR=PCA_HEFS(Ib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 对比
% pIGS=unwrap2(pwIGS);
% pGS=unwrap2(pwGS);
pAIA=unwrap2(pwAIA);
pMSSM=unwrap2(pwMSSM);
pLSR=unwrap2(pwLSR);

pIGS=(pwIGS);
pGS=(pwGS);
% pAIA=(pwAIA);
% pMSSM=-(pwMSSM);

pIGS=pIGS-pIGS(1,1);
pGS=pGS-pGS(1,1);
pAIA=pAIA-pAIA(1,1);
pMSSM=pMSSM-pMSSM(1,1);
pLSR=pLSR-pLSR(1,1);

eIGS=pIGS-p;
eGS=pGS-p;
eAIA=pAIA-p;
eMSSM=pMSSM-p;
eLSR=pLSR-p;

RMSE_IGS(i)=std2(eIGS);
RMSE_GS(i)=std2(eGS);
RMSE_AIA(i)=std2(eAIA);
RMSE_MSSM(i)=std2(eMSSM);
RMSE_LSR(i)=std2(eLSR);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 显示
% figure
% subplot 131
% imshow(Ib(:,:,1),[])
% subplot 132
% imshow(Ib(:,:,2),[])
% subplot 133
% imshow(Ib(:,:,3),[])

figure
plot(pss,RMSE_IGS,'bo-');
hold on
% plot(pss,RMSE_GS,'r*-');
plot(pss,RMSE_AIA,'gx-');
plot(pss,RMSE_MSSM,'k.-');
plot(pss,RMSE_LSR,'r*-');
legend('SONIA','AIA','MSSM','HEFS')
xlabel('phase shift /rad ')
ylabel('RMSE /rad ')


% figure
% plot(psr,RMSE_IGS,'bo-');
% hold on
% plot(psr,RMSE_GS,'r*-');
% plot(psr,RMSE_AIA,'gx-');
% plot(psr,RMSE_MSSM,'k.-');
% legend('IGS','GS','AIA','MSSM')


% figure
% subplot 331
% imagesc(p);title('phase of SET');colorbar;xlabel('(a)')
% subplot 332
% imagesc(pIGS);title('phase of IGS');colorbar;xlabel('(b)')
% subplot 333
% imagesc(pGS);title('phase of GS');colorbar;xlabel('(c)')
% subplot 334
% imagesc(pAIA);title('phase of AIA');colorbar;xlabel('(d)')
% subplot 335
% imagesc(pMSSM);title('phase of MSSM');colorbar;xlabel('(e)')
% subplot 336
% imagesc(eIGS);title('error of IGS');colorbar;xlabel('(f)')
% subplot 337
% imagesc(eGS);title('error of GS');colorbar;xlabel('(g)')
% subplot 338
% imagesc(eAIA);title('error of AIA');colorbar;xlabel('(h)')
% subplot 339
% imagesc(eMSSM);title('error of MSSM');colorbar;xlabel('(i)')

