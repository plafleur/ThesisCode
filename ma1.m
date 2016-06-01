function [mat1]=ma1(k);
N=1000;
n1=(1:N);
n2=(-1:-N);
if k==0;
    mat1=2*(sum(1./(n1.^3))+sum(1./(abs(n2).^3)));
    
else
    kna1=n1*k;
    kna2=n2*k;
    
    ck1=2*cos(kna1);
     %ck1=exp(-i*kna1);
    
    
    
    
    mat1=(sum(ck1./abs(n1).^3));
    
end

end

