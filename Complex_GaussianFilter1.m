function [I]=Complex_GaussianFilter1(I,freq,freq1)

[NR, NC]=size(I);
[u,v]=meshgrid(1:NC, 1:NR);
u0=floor(NC/2)+1; 
v0=floor(NR/2)+1;           
u=u-u0; 
v=v-v0;       
H=1-exp(-(u.^2+v.^2)/(2*(freq)^2));
H1=exp(-(u.^2+v.^2)/(2*(freq1)^2));
H=H.*H1;
C=fft2(I);    
CH=C.*ifftshift(H);    
I=ifft2(CH);       