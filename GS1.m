function pw=GS1(Im1,Im2)
%% GS方法，这里Im1和Im2都是去除背景项A后的两幅调制强度分布图；pw是得到的包裹相位
norm1=sqrt(sum(sum(Im1.*Im1)));
Im1=Im1./norm1;                                    %对Im1进行归一化
proj=sum(sum((Im2.*Im1)))/sum(sum(Im1.*Im1));
Im2=Im2-proj*Im1;                                  %对Im2进行正交化
norm2=sqrt(sum(sum(Im2.*Im2))); 
Im2=Im2./norm2;                                    %对Im2进行归一化
pw=atan2(-Im2,Im1);                                %求得的包裹相位                             

