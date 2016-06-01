clear
n=101;
q=1.602e-19;
e0=8.854e-12;
l=1.05e-15;
p=3.21e-11;
m1=3.81754035e-26;
m2=5.81e-26;
M=m1*m2/(m1+m2);
mu1=sqrt(M/m1);
mu2=sqrt(M/m2);
c1=zeros(1,n);
c2=zeros(1,n);
k=zeros(1,n);
Q1=1;
Q2=-1;
D0=72.2867019009821;
R0=3.192e-10;
for c=1:n;
	k(c)=(-pi+0.001)+(2*pi)*(c-1)/(n-1);
	Phi1=D0-(q^2/(2*pi*e0*(2*R0)^3))*ma1(k(c));
	Phi2=[(q^2/(2*pi*e0*(2*R0)^3))*ma2(k(c))-2*(l/p^2)*exp(-(R0/p))*cos(k(c)/2)];
	D11=mu1^2*Phi1;
	D12=mu1*mu2*Phi2;
	D22=mu2^2*Phi1;
	A=[D11 D12; D12 D22];
	[b,lambda]=eig(A);
	c1(c)=Q1*mu1*b(1)+Q2*mu2*b(3);
	c2(c)=Q1*mu1*b(3)+Q2*mu2*b(4);

end
