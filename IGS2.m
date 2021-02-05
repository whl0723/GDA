function [pw2,ii,ee,pw0]=IGS2(I1,I2)
I1=double(I1);
I2=double(I2);
Im1=I1;
Im2=I2;
% tic
pw1=GS1(I1,I2);
% tgs=toc
pw0=pw1;
% p1=unwrap2(pw1);
% p1=p1-p1(1,1);
ii=1;
ee=1;
norm1=sqrt(sum(sum(Im1.*Im1)));
Im1=Im1./norm1;                                    %对Im1进行归一化
proj=sum(sum((Im2.*Im1)))/sum(sum(Im1.*Im1));
Im2=Im2-proj*Im1;                                  %对Im2进行正交化
norm2=sqrt(sum(sum(Im2.*Im2)));
Im2=Im2./norm2;
while ee>1e-16  && ii<1000
    e1=(sin(pw1)).^2;e=sum(e1(:));
    b1=(cos(pw1)).^2;b=sum(b1(:));
    c1=(sin(pw1)).*(cos(pw1));c=sum(c1(:));
    m1=(e-(c*c)/b)^0.5;
    m2=b^0.5;
    m3=c/b;
    
    pw2=atan2(real((-Im2)*m1+m3*m2*Im1),real(Im1*m2));
    ee=std2(pw2-pw1);
    pw1=pw2;
    ii=ii+1;
end
% tigs=toc