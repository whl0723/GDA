% ���룺���ͼI(:,:,1),I(:,:,2)...

function phi=PCA_HEFS(I)

[h,l,n]=size(I); k=h*l;
X=zeros(n,k);
for  i=1:n
    X(i,:)=reshape(I(:,:,i),1,[]);     % X��ÿ�б�ʾһ������ͼ��ÿ�б�ʾ��ͬ���ص㴦��ͬ����ͼ��Ӧ��ǿ��ֵ
end
X = X-ones(size(X,1),1)*mean(X);       %ȥ������ǿ�ȣ� mean(X)��ʾ�Ծ���X��ÿ����ƽ��ֵ
covarianceMatrix = X*X';               %�����Э�������
[U,~,~] = svd(covarianceMatrix);  
Up=zeros(2,n);
Up(1,:)=U(:,1)';
Up(2,:)=U(:,2)';
Y = Up*X;

% ������й�һ��,�����������Ǵ��¾��ȷֲ���[0��2*pi]��ʱ������������̿�����߾���
srss1=sqrt(sum(Y(1,:).*Y(1,:)));
V(:,1)=Y(1,:)./srss1;
srss2=sqrt(sum(Y(2,:).*Y(2,:)));
V(:,2)=Y(2,:)./srss2;
[alpha,beta]=semiHyperLeastSquareFitInsection3C(V);
pp=ellipseFitPhaseRetrievalInSection3B(V,alpha,beta);

phi=reshape(pp,size(I(:,:,1)));
end

function [alpha,b]=semiHyperLeastSquareFitInsection3C(V)
N=size(V,1);
b=max(abs(V(:)));
x=V(:,1);
y=V(:,2);
Chi=[x.^2,2*x.*y,y.^2,2*b*x,2*b*y,b^2*ones(N,1)];

X=mean(x);
Y=mean(y);
XY=mean(Chi(:,2))/2;
X2=mean(Chi(:,1));
Y2=mean(Chi(:,3));

W=[6*X2  6*XY      X2+Y2 6*b*X 6*b*Y b^2;
   6*XY  4*(X2+Y2) 6*XY  4*b*Y 4*b*X 0;
   X2+Y2 6*XY      6*Y2  2*b*X 6*b*Y b^2;
   6*b*X 4*b*Y     2*b*X 4*b^2 0     0;
   2*b*Y 4*b*X     6*b*Y 0     4*b^2 0;
   b^2   0         b^2   0     0     0;];
opt.issym=true;
X=Chi'*Chi/N;
[alpha,~] = eigs(W,X,1,'lm',opt);
end

function phi=ellipseFitPhaseRetrievalInSection3B(V,a,b)
x=V(:,1)-b*(a(3)*a(4)-a(2)*a(5))/(a(2)^2-a(1)*a(3));
y=V(:,2)-b*(a(1)*a(5)-a(2)*a(4))/(a(2)^2-a(1)*a(3));
k=sqrt(4*a(2)^2+(a(1)-a(3))^2);
t=(a(1)-a(3)+k)/(2*a(2));

re=sqrt(abs(a(1)+a(3)+k))*(y+t*x);
im=sqrt(abs(a(1)+a(3)-k))*(x-t*y);
phi=angle(complex(re,im));
end