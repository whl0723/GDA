function pw=GS1(Im1,Im2)
%% GS����������Im1��Im2����ȥ��������A�����������ǿ�ȷֲ�ͼ��pw�ǵõ��İ�����λ
norm1=sqrt(sum(sum(Im1.*Im1)));
Im1=Im1./norm1;                                    %��Im1���й�һ��
proj=sum(sum((Im2.*Im1)))/sum(sum(Im1.*Im1));
Im2=Im2-proj*Im1;                                  %��Im2����������
norm2=sqrt(sum(sum(Im2.*Im2))); 
Im2=Im2./norm2;                                    %��Im2���й�һ��
pw=atan2(-Im2,Im1);                                %��õİ�����λ                             

