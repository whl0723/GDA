% 本程序生成随机相位曲面，随机相移量，随机相位波动范围，高斯噪声随机35-45dB


clc
clear
close all
j=1;
for i=1:1000
    p=rand(30,3);  %
    x=p(:,1);y=p(:,2);z=p(:,3);
    [xi,yi]=meshgrid(linspace(min(x),max(x),256),linspace(min(y),max(y),256));
    zi=griddata(x,y,z,xi,yi,'v4');
    dxzi=zi(1:end-1,:)-zi(2:end,:);
    dyzi=zi(:,1:end-1)-zi(:,2:end);
    index_dxzi=dxzi>pi;
    index_dyzi=dyzi>pi;
    
    if sum(index_dxzi(:))<1&&sum(index_dyzi(:))<1
        
        % h=scatter3(x,y,z,30,5*ones(size(x)),'r','filled');%,'filled'??????
        p0=zi-min(zi(:));
        A=100*exp(-0.05*(xi.^2+yi.^2));  %背景项
        B=80*exp(-0.01*(xi.^2+yi.^2));   %调制项
        
        range=rand(1)*5+0.2;
        p=p0./max(p0(:))*pi*range;
        pw=mod(p,2*pi)-pi;
        theta1=0.5+rand(1)*2;
        theta2=0.7+rand(1)*2;
        
        Ib(:,:,1)=A+B.*cos(p);
        Ib(:,:,2)=A+B.*cos(p+theta1);
        Ib(:,:,3)=A+B.*cos(p+theta1+theta2);
        
        noise=35+10*rand(1);
        Ib(:,:,1)=awgn(Ib(:,:,1),noise,'measured');
        Ib(:,:,2)=awgn(Ib(:,:,2),noise,'measured');
        Ib(:,:,3)=awgn(Ib(:,:,3),noise,'measured');
        
        I1=Ib(:,:,1)-Ib(:,:,2);
        I2=Ib(:,:,1)-Ib(:,:,3);
        % pwGS=GS1(I1,I2);
        [pwIGS,ii,ee,pwGS]=IGS2(I1,I2);
        %%%%%%%%%%%%%%%% AIA
        ps=linspace(0,3,3);MaxIter=100;
        %     tic
        [pwAIA,delta,iter] = S_AIA(Ib,ps,MaxIter);
        %     taia=toc;
        
        %%%%%%%%%%%%%%%%%%% MSSM
        fl=1.5;
        fh=60;
        %     tic
%         pwMSSM=wy3(Ib(:,:,1),Ib(:,:,2),Ib(:,:,3),fl,fh);
        [pwMSSM,~]=UV_factorization2(I1,I2);
        %     tmssm=toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     tic
        pwLSR=PCA_HEFS(Ib);
        % [Dt,Nt,pwLSR]=PCA_LEF(Ib(:,:,1),Ib(:,:,2))
        %     tLSR=toc;
        %%%%%%%%%%%%%%%%%%%%%
        
        pIGS=unwrap2(pwIGS);
        pGS=unwrap2(pwGS);
        pAIA=unwrap2(pwAIA);
        pMSSM=unwrap2(pwMSSM);
        pLSR=unwrap2(pwLSR);
        
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
        
        
        R_noise(i)=noise;
        db(i)=RMSE_IGS(i)./R_noise(i)*100;
        if RMSE_IGS(i)>0.06
            err(:,:,j)=pIGS-p;
            pp(:,:,j)=p;
            pp_out(:,:,j)=pIGS;
            R_theta1(j)=theta1;
            R_theta2(j)=theta2;
            R_interfergram(:,:,j)=Ib(:,:,1);
            j=j+1;
        end
    end
    
end




% plot(RMSE_LSR,'g--')
% legend('GDA','AIA','MSSM','HEFS')
% xlabel('Times of Calculations')
% ylabel('RMSE /rad')
% 
% figure
% plot(R_noise)


figure
subplot 221
plot(RMSE_IGS,'r*-');title('GDA')
xlabel('Times of Calculations')
ylabel('RMSE /rad')
subplot 222
plot(RMSE_AIA,'b.-');title('AIA')
xlabel('Times of Calculations')
ylabel('RMSE /rad')
subplot 223
plot(RMSE_MSSM,'ko-');title('MSSM')
xlabel('Times of Calculations')
ylabel('RMSE /rad')
subplot 224
plot(RMSE_LSR,'g--');title('HEFS')
xlabel('Times of Calculations')
ylabel('RMSE /rad')


figure
subplot 221
plot(RMSE_IGS,'r*-');title('GDA')
subplot 222
plot(RMSE_AIA,'b.-');title('AIA')
subplot 223
plot(RMSE_MSSM,'ko-');title('UV2')
subplot 224
plot(RMSE_LSR,'g--');title('HEFS')

if j>1
    figure
    imagesc(pp(:,:,1))
    figure
    imagesc(pp_out(:,:,1))
    figure
    imagesc(err(:,:,1))
end

Good_IGS=sum(RMSE_IGS<0.08)
Bad_IGS=sum(RMSE_IGS<0.2)-Good_IGS
Wrong_IGS=sum(RMSE_IGS>=0.2)

Good_AIA=sum(RMSE_AIA<0.08)
Bad_AIA=sum(RMSE_AIA<0.2)-Good_AIA
Wrong_AIA=sum(RMSE_AIA>=0.2)

% Good_MSSM=sum(RMSE_MSSM<0.08)
% Bad_MSSM=sum(RMSE_MSSM<0.2)-Good_MSSM
% Wrong_MSSM=sum(RMSE_MSSM>=0.2)

Good_UV2=sum(RMSE_MSSM<0.08)
Bad_UV2=sum(RMSE_MSSM<0.2)-Good_UV2
Wrong_UV2=sum(RMSE_MSSM>=0.2)

Good_LSR=sum(RMSE_LSR<0.08)
Bad_LSR=sum(RMSE_LSR<0.2)-Good_LSR
Wrong_LSR=sum(RMSE_LSR>=0.2)


% figure
% imagesc(p)
% figure
% imagesc(pIGS-p)
