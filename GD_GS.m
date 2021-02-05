function [pw2,RR1,alph]=GD_GS(Im1,Im2,pw1)
% pw1=GS1(Im1,Im2);

%%
norm1=sqrt(sum(sum(Im1.*Im1)));
Im11=Im1./norm1;                                    %对Im1进行归一化
proj=sum(sum((Im2.*Im11)))/sum(sum(Im11.*Im11));
Im22=Im2-proj*Im11;                                  %对Im2进行正交化
norm2=sqrt(sum(sum(Im22.*Im22)));
Im22=Im22./norm2;
S=-Im22./Im11;
% bb=1;
%%
[M,N]=size(Im1);
K=M*N;
i=2;
ee(1)=1;
RR=ones(100,1);
alph(1)=0.01;
while i<40  && ee>1e-5
    m1=sqrt(sum(cos(pw1(:)).^2,'omitnan'));
    m2=sum(sin(pw1(:)).*cos(pw1(:)),'omitnan')./(m1.^2);
    m3=sqrt(K-m1*m1-m1*m1*m2*m2);
    
    dm1=-cos(pw1).*sin(pw1)./m1;
    dm2=2*m2/m1/m1*cos(pw1).*sin(pw1)+1/m1/m1*(cos(pw1).^2-sin(pw1).^2);
    dm3=(1-m2^2)/m3*cos(pw1).*sin(pw1)-m2/m3*(cos(pw1).^2-sin(pw1).^2);
    
    dm3m1=dm3.*m1-m3*dm1;
    pw2=atan2(-Im22*m3+Im11*m1*m2,real(Im11*m1)+eps);
    f=pw2-pw1;
    RR(i)=std(rmmissing(f(:)));
    ee=abs(RR(i)-RR(i-1));
    df=(S.*dm3m1+dm2)./(1+(S.*m3/m1+m2).^2)-1;
    bb=(sum(df(:),'omitnan')+eps);
%     bb=1+i*5;
    alph=alph./bb;
    pw1=pw2-alph*df;
    i=i+1;
end
RR1=RR(2:i-1);