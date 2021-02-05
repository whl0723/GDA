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

% p0=peaks(300) ;          %预设相位 
% p0=p0-min(p0(:));
% p=p0./max(p0(:))*8*pi;
% p=p-p(1,1);
p0=255-double(rgb2gray(imread('cell.jpg')));
p0(p0>115)=115;
p0=p0./max(p0(:))*pi*0.8;
p = imresize(p0,size(A),'bilinear');

Ib(:,:,1)=A+B.*cos(p);
Ib(:,:,2)=A+B.*cos(p+2);
Ib(:,:,3)=A+B.*cos(p+3.3);

noise=40;
Ib(:,:,1)=awgn(Ib(:,:,1),noise,'measured');
Ib(:,:,2)=awgn(Ib(:,:,2),noise,'measured');
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
tic
% [pwIGS,ii,ee,pwGS]=IGS2(I1,I2);
[pwIGS,RR,alph]=GD_GS(I1,I2,pwGS);
tigs=toc
% figure
% plot(RR)



%%%%%%%%%%%%%%%% AIA
ps=linspace(0,3,3);MaxIter=500;
tic
[pwAIA,delta,iter] = S_AIA(Ib,ps,MaxIter);
taia=toc

%%%%%%%%%%%%%%%%%%% MSSM
fl=1.5;
fh=60;
tic
pwMSSM=wy3(Ib(:,:,1),Ib(:,:,2),Ib(:,:,3),fl,fh);
tmssm=toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
pwLSR=PCA_HEFS(Ib);
% [Dt,Nt,pwLSR]=PCA_LEF(Ib(:,:,1),Ib(:,:,2))
tLSR=toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 对比
pIGS=unwrap2(pwIGS);
pGS=unwrap2(pwGS);
pAIA=unwrap2(pwAIA);
pMSSM=unwrap2(pwMSSM);
pLSR=unwrap2(pwLSR);

pIGS=pIGS-min(pIGS(:));
pGS=pGS-min(pGS(:));
pAIA=pAIA-min(pAIA(:));
pMSSM=pMSSM-min(pMSSM(:));
pLSR=pLSR-min(pLSR(:));

figure
imagesc(pLSR)

eIGS=pIGS-p;
eGS=pGS-p;
eAIA=pAIA-p;
eMSSM=pMSSM-p;
eLSR=pLSR-p;

eIGS=eIGS-mean2(eIGS);
eGS=eGS-mean2(eGS);
eAIA=eAIA-mean2(eAIA);
eMSSM=eMSSM-mean2(eMSSM);
eLSR=eLSR-mean2(eLSR);


RMSE_IGS=std2(eIGS)
RMSE_GS=std2(eGS)
RMSE_AIA=std2(eAIA)
RMSE_MSSM=std2(eMSSM)
RMSE_LSR=std2(eLSR)

% PVRMSE_IGS=std2(eIGS)./(max(pIGS(:))-min(pIGS(:)))
% PVRMSE_GS=std2(eGS)./(max(pGS(:))-min(pGS(:)))
% PVRMSE_AIA=std2(eAIA)./(max(pAIA(:))-min(pAIA(:)))
% PVRMSE_MSSM=std2(eMSSM)./(max(pMSSM(:))-min(pMSSM(:)))
% PVRMSE_LSR=std2(eLSR)./(max(pLSR(:))-min(pLSR(:)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 显示
% figure
% subplot 131
% imshow(Ib(:,:,1),[])
% subplot 132
% imshow(Ib(:,:,2),[])
% subplot 133
% imshow(Ib(:,:,3),[])


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

%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot 421
% imagesc(p);title({'Phase of';'SET'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(a)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 422
% imagesc(Ib(:,:,1));title({'One of';'Interferngrams'},'FontSize',9);colorbar;
% text(25,25,'(b)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 423
% imagesc(pAIA);title({'Phase of';'AIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(c)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 424
% imagesc(eAIA);title({'Phase Deviation of';'AIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(d)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 425
% imagesc(pLSR);title({'Phase of';'HEFS'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(e)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 426
% imagesc(eLSR);title({'Phase Deviation of';'HEFS'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(f)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 427
% imagesc(pIGS);title({'Phase of';'SONIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(g)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 428
% imagesc(eIGS);title({'Phase Deviation of';'SONIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(h)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')


figure
subplot 231
imagesc(p);title({'Phase of';'SET'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(25,25,'(a)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 232
imagesc(Ib(:,:,1));title({'One of';'Interferngrams'},'FontSize',9);colorbar;
text(25,25,'(b)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 233
imagesc(eIGS);title({'Phase Deviation of';'SCOI'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(25,25,'(c)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 234
imagesc(eAIA);title({'Phase Deviation of';'AIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(25,25,'(d)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 235
imagesc(eMSSM);title({'Phase Deviation of';'MSSM'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(25,25,'(e)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 236
imagesc(eLSR);title({'Phase Deviation of';'HEFS'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(25,25,'(f)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')

