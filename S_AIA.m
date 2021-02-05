function [pw,delta,iter] = S_AIA(Is,ps,MaxIter)
%% ˵��    ��������������У��Ľ�����С���˵����㷨����������ͼ
% ���룺
% Is��Is(:,:,i)��ʾÿ������ͼ�Ĺ�ǿ�ֲ�����K��
% ps�Ǽ��裨��õ�������������һ��������,����������������ŷ��ˣ���ô��õ���λҲ�Ƿ���
% MaxIter���趨������������
% �����
% pw����õİ�����λ��  delta����õ���������  iter�Ǵﵽ������ﵽ��������������ĵ�������
%% 
epsilon = 1e-3;        % �趨������������ֵ
err = 1;               % �������ε����������仯����ʾ
iter = 0;              % ����������ʾ
dim = size(Is);  
K = dim(1)*dim(2);       %������
N = dim(3);              %ͼ��
Ip = zeros(N,K);
for  i=1:N
    Ip(i,:) = reshape(Is(:,:,i),1,K);%��ÿ������ͼ�ֱ��ʾ���������������ϳɾ���
end  
Ip1=Ip';
phi = zeros(1,K);        %������λ���ڼ���ı�ʾ
% Mod = zeros(1,K);
E=eye(3);
%%
%%
while  ((err > epsilon) && (iter < MaxIter))         
%Fisrt setp:        
%Ϊ��ʹAΪ���������Ҫ����������ͼ����һ����������֪
        ps1 = ps;
        A = [ N         ,           sum(cos(ps)),             sum(sin(ps))
            sum(cos(ps)),           sum(cos(ps).^2),          sum(cos(ps).*sin(ps))
            sum(sin(ps)),           sum(sin(ps).*cos(ps)),    sum(sin(ps).^2)];
        A1=A\E;
%����ͼ��λֵ���� 
        s1=sum(Ip);          % size: [1,K]
        s2=cos(ps)*Ip;
        s3=sin(ps)*Ip;    
        s=[s1;s2;s3];
        X=A1*s;
        phi=atan2(-X(3,:),X(2,:));
         %  Mod(i) = sqrt(X(3,:).*X(3,:)+X(2,:).*X(2,:));
        
%Second  step:
%��һ����λ��֪
        A = [  K         ,         sum(cos(phi)),               sum(sin(phi))
            sum(cos(phi)),         sum(cos(phi).^2 ),           sum(cos(phi).*sin(phi))
            sum(sin(phi)),         sum(sin(phi).*cos(phi)),     sum(sin(phi).^2)];
        A1=A\E;
%����������   
        s1=sum(Ip1);          % size: [1,N]
        s2=cos(phi)*Ip1;
        s3=sin(phi)*Ip1;  
        s=[s1;s2;s3];
        X=A1*s;
        ps=atan2(-X(3,:),X(2,:));
        
        err = max(abs(ps-ps1));
        iter = iter+1;
end
%%
%%
pw = reshape(phi,dim(1),dim(2));
% pwMod = reshape(Mod,dim(1),dim(2));   
delta=ps;
for j=2:N
        delta(j)=delta(j)-2*pi*round((delta(j)-delta(j-1))/(2*pi));            %�������������,Ҫ����������ͼ����ʵ������֮��ô���pi����Ȼ�޷����      
end
delta = delta-delta(1);
if  delta(N)<0
   delta = -delta;
   pw=-pw;
end