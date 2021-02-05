function [pw,delta,iter] = S_AIA(Is,ps,MaxIter)
%% 说明    单波长干涉测量中，改进的最小二乘迭代算法，至少三幅图
% 输入：
% Is即Is(:,:,i)表示每幅干涉图的光强分布，共K幅
% ps是假设（或得到）的相移量，一个行向量,若假设的相移量符号反了，那么求得的相位也是反的
% MaxIter是设定的最大迭代次数
% 输出：
% pw是求得的包裹相位；  delta是求得的相移量；  iter是达到收敛或达到最大迭代次数所需的迭代次数
%% 
epsilon = 1e-3;        % 设定相移量精度阈值
err = 1;               % 相邻两次迭代相移量变化量表示
iter = 0;              % 迭代次数表示
dim = size(Is);  
K = dim(1)*dim(2);       %像素数
N = dim(3);              %图数
Ip = zeros(N,K);
for  i=1:N
    Ip(i,:) = reshape(Is(:,:,i),1,K);%将每幅干涉图分别表示成行向量，多幅组合成矩阵
end  
Ip1=Ip';
phi = zeros(1,K);        %设置相位用于计算的表示
% Mod = zeros(1,K);
E=eye(3);
%%
%%
while  ((err > epsilon) && (iter < MaxIter))         
%Fisrt setp:        
%为了使A为非奇异矩阵，要求至少三幅图，这一步相移量已知
        ps1 = ps;
        A = [ N         ,           sum(cos(ps)),             sum(sin(ps))
            sum(cos(ps)),           sum(cos(ps).^2),          sum(cos(ps).*sin(ps))
            sum(sin(ps)),           sum(sin(ps).*cos(ps)),    sum(sin(ps).^2)];
        A1=A\E;
%干涉图相位值计算 
        s1=sum(Ip);          % size: [1,K]
        s2=cos(ps)*Ip;
        s3=sin(ps)*Ip;    
        s=[s1;s2;s3];
        X=A1*s;
        phi=atan2(-X(3,:),X(2,:));
         %  Mod(i) = sqrt(X(3,:).*X(3,:)+X(2,:).*X(2,:));
        
%Second  step:
%这一步相位已知
        A = [  K         ,         sum(cos(phi)),               sum(sin(phi))
            sum(cos(phi)),         sum(cos(phi).^2 ),           sum(cos(phi).*sin(phi))
            sum(sin(phi)),         sum(sin(phi).*cos(phi)),     sum(sin(phi).^2)];
        A1=A\E;
%相移量计算   
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
        delta(j)=delta(j)-2*pi*round((delta(j)-delta(j-1))/(2*pi));            %用于相移量解包,要求相邻两幅图的真实相移量之差不得大于pi，不然无法解包      
end
delta = delta-delta(1);
if  delta(N)<0
   delta = -delta;
   pw=-pw;
end