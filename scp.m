clear;
n=100;
q=linspace(0,pi,n);
cosq=cos(q);
d=1;

%100 q-point parameters
C=3;
A=21;
B=3.1e23;
Q=8*1.602e-19;

V=(3.905e-10)^3;
e0=8.854e-12;

M=4e-26;
dx2old=(abs(A-12*C)+0.5)/B;
dx2new=1e10;
diff=abs(dx2old-dx2new);
wsum=0;
E=0;
hbar=1.054571e-34;
kb=1.308e-23;
con=Q^2/(e0*V);
hbar2kb=hbar/(2*kb);
con2=hbar/(2*M*n^3);
T=zeros(1,d);
omega=zeros(1,d);
chi=zeros(1,d);
con3=Q*E;	
z=0;
for k=1:d;
    T(k)=300;
    clc
    disp(k)
	while diff>1e-27;
  		for x=1:n;
		  	for y=1:n;
           			  F=cosq(x)+cosq(y)+cosq(:);     
			  		  w=sqrt(1/M)*sqrt(A-4*F*C+B*dx2old);
                      wsum=wsum+sum((1./w).*coth((hbar2kb.*w)./T(k)));		
	 	 end
  	end
    wq0=(A-12*C+B*dx2old);
  	dx2new=con2*wsum+(1/3)*(con3/wq0)^2;
    diff=abs(dx2old-dx2new);
  	dx2old=0.7*dx2old+0.3*dx2new;  
    %B*dx2old
	wsum=0;
    end
dx2(k)=dx2new;
  omega0=(A-12*C+B*dx2old);
  %omega=(A-4.*C.*(cosq+cos(0)+cos(0))+B*dx2old);
  chi(k)=con./(omega0+(2*B/3)*(Q*E/omega0)^2);
  diff=1;
end

