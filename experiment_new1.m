clc
clear
close all
%%%%%%%%%%%%%%% 参考相位
load('C:\Users\whl\Pictures\试验数据\SONIA\0619\3\REFm.mat');
p=AIA_phase-min(AIA_phase(:));

A=double(imread('C:\Users\whl\Pictures\试验数据\SONIA\0619\3\A\CCD1\768403.bmp'));
B=double(imread('C:\Users\whl\Pictures\试验数据\SONIA\0619\3\B\CCD1\770142.bmp'));
AB=sqrt(A.*B);

%%%%%%%%%%%%% 计算带物体的干涉图的相位分布
index=[1,30,70];

path1=strcat('C:\Users\whl\Pictures\试验数据\SONIA\0619\3\g\CCD1\',num2str(index(1)),'.bmp');
path2=strcat('C:\Users\whl\Pictures\试验数据\SONIA\0619\3\g\CCD1\',num2str(index(2)),'.bmp');
path3=strcat('C:\Users\whl\Pictures\试验数据\SONIA\0619\3\g\CCD1\',num2str(index(3)),'.bmp');

Is(:,:,1)=double(imread(path1));
Is(:,:,2)=double(imread(path2));
Is(:,:,3)=double(imread(path3));

% Ib(:,:,1)=(Is(1:end,1:end,1)-A-B)./AB;
% Ib(:,:,2)=(Is(1:end,1:end,2)-A-B)./AB;
% Ib(:,:,3)=(Is(1:end,1:end,3)-A-B)./AB;

    Ib(:,:,1)=Is(1:end,1:end,1);
    Ib(:,:,2)=Is(1:end,1:end,2);
    Ib(:,:,3)=Is(1:end,1:end,3);

% Ic(:,:,1)=Is(1:end,1:end,1);
% Ic(:,:,2)=Is(1:end,1:end,2);
% Ic(:,:,3)=Is(1:end,1:end,3);
% 
% % x1=231;
% % x2=414;
% % y1=127;
% % y2=442;
% 
% x1=1;
% x2=450;
% y1=1;
% y2=450;
% 
% Ib(:,:,1)=Ic(x1:x2,y1:y2,1);
% Ib(:,:,2)=Ic(x1:x2,y1:y2,2);
% Ib(:,:,3)=Ic(x1:x2,y1:y2,3);

%%%%%%%%%%%%% IGS GS
I1=Ib(:,:,1)-Ib(:,:,2);
I2=Ib(:,:,1)-Ib(:,:,3);
pwGS=GS1(I1,I2);
tic
% [pwIGS,ii,ee,pwGS]=IGS2(I1,I2);
[pwIGS,RR,alph]=GD_GS(I1,I2,I1);
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

% %%%%%%%%%%%%%%%%%%%%%%%%%% 计算干涉图背景信息
% path1b=strcat('H:\whl\0618\50\CCD1\',num2str(index(1)),'.bmp');
% path2b=strcat('H:\whl\0618\50\CCD1\',num2str(index(2)),'.bmp');
% path3b=strcat('H:\whl\0618\50\CCD1\',num2str(index(3)),'.bmp');
% 
% Isb(:,:,1)=double(imread(path1b));
% Isb(:,:,2)=double(imread(path2b));
% Isb(:,:,3)=double(imread(path3b));
% 
% % Ib(:,:,1)=Is(10:end,10:end,1);
% % Ib(:,:,2)=Is(10:end,10:end,2);
% % Ib(:,:,3)=Is(10:end,10:end,3);
% 
% Icb(:,:,1)=Isb(10:end,10:end,1);
% Icb(:,:,2)=Isb(10:end,10:end,2);
% Icb(:,:,3)=Isb(10:end,10:end,3);
% 
% Ibb(:,:,1)=Icb(x1:x2,y1:y2,1);
% Ibb(:,:,2)=Icb(x1:x2,y1:y2,2);
% Ibb(:,:,3)=Icb(x1:x2,y1:y2,3);
% 
% %%%%%%%%%%%%% IGS GS --- 背景
% I1b=Ibb(:,:,1)-Ibb(:,:,2);
% I2b=Ibb(:,:,1)-Ibb(:,:,3);
% pwGSb=GS1(I1b,I2b);
% tic
% [pwIGSb,ii,ee,pwGSb]=IGS2(I1b,I2b);
% % [pwIGS,RR,alph]=GD_GS(I1,I2,pwGS);
% tigs=toc;
% % figure
% % plot(RR)
% %%%%%%%%%%%%%%%% AIA --- 背景
% ps=linspace(0,3,3);MaxIter=500;
% tic
% [pwAIAb,deltab,iter] = S_AIA(Ibb,ps,MaxIter);
% taia=toc;
% 
% %%%%%%%%%%%%%%%%%%% MSSM
% fl=1.5;
% fh=60;
% tic
% pwMSSMb=wy3(Ibb(:,:,1),Ibb(:,:,2),Ibb(:,:,3),fl,fh);
% tmssm=toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 对比
pIGS=unwrap2(pwIGS);
pGS=unwrap2(pwGS);
pAIA=unwrap2(pwAIA);
pMSSM=-unwrap2(pwMSSM);
pLSR=unwrap2(pwLSR);

% pIGS=unwrap2(pwIGS)-unwrap2(pwIGSb);
% pGS=unwrap2(pwGS)-unwrap2(pwIGSb);
% pAIA=unwrap2(pwAIA)-unwrap2(pwIGSb);
% pMSSM=unwrap2(pwMSSM)-unwrap2(pwIGSb);

% pIGS=unwrap2(pwIGS)-pb;
% pGS=unwrap2(pwGS)-pb;
% pAIA=unwrap2(pwAIA)-pb;
% pMSSM=unwrap2(pwMSSM)-pb;

pIGS=pIGS-min(pIGS(:));
pGS=pGS-min(pGS(:));
pAIA=pAIA-min(pAIA(:));
pMSSM=pMSSM-min(pMSSM(:));
pLSR=pLSR-min(pLSR(:));

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

RMSE_IGS=std2(eIGS);
RMSE_GS=std2(eGS);
RMSE_AIA=std2(eAIA);
RMSE_MSSM=std2(eMSSM);
RMSE_LSR=std2(eLSR);

% RMSE_IGS(iii)=std2(eIGS);
% RMSE_GS(iii)=std2(eGS);
% RMSE_AIA(iii)=std2(eAIA);
% RMSE_MSSM(iii)=std2(eMSSM);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 显示
% figure
% subplot 131
% imshow(Ib(:,:,1),[])
% subplot 132
% imshow(Ib(:,:,2),[])
% subplot 133
% imshow(Ib(:,:,3),[])

% figure
% plot(RMSE_IGS)
% hold on
% plot(RMSE_AIA)
% plot(RMSE_MSSM)
% % plot(RMSE_IGS)
% legend('IGS','AIA','MSSM')
% figure
% plot(RMSE_MSSM)

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
% imagesc(pMSSM);title({'Phase of';'MSSM'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(e)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 426
% imagesc(eMSSM);title({'Phase Deviation of';'MSSM'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(f)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 427
% imagesc(pIGS);title({'Phase of';'SONIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(g)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
% subplot 428
% imagesc(eIGS);title({'Phase Deviation of';'SONIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
% text(25,25,'(h)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
aaa=45;
figure
subplot 231
imagesc(p);title({'REF'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(aaa,aaa,'(a)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 232
imagesc(Ib(:,:,1));title({'One of';'Interferngrams'},'FontSize',9);colorbar;
text(aaa,aaa,'(b)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 233
imagesc(eIGS);title({'Phase Deviation of';'SCOI'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(aaa,aaa,'(c)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 234
imagesc(eAIA);title({'Phase Deviation of';'AIA'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(aaa,aaa,'(d)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 235
imagesc(eMSSM);title({'Phase Deviation of';'MSSM'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(aaa,aaa,'(e)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')
subplot 236
imagesc(eLSR);title({'Phase Deviation of';'HEFS'},'FontSize',9);colorbar;hcb=colorbar;title(hcb,'rad')
text(aaa,aaa,'(f)','horiz','center','BackgroundColor','w','color','b','fontsize',9,'FontWeight','bold')