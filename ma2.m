function [mat2] =ma2(k)
N=1000;
n1=(1:N);
n2=(-N:0);
if k==0;
   mat2=round(10^5*(sum(1./(n1-0.5).^3) +sum(1./abs(n2-0.5).^3)))/10^5;
else
kna1=k*(n1-0.5);
kna2=k*(n2-0.5);
ck1=exp(-i*kna1);
ck2=exp(-i*kna2);
mat2=round(10^5*(sum(ck1./((n1-0.5).^3))+sum(ck2./(abs(n2-0.5).^3))))/10^5;
   
end




end

