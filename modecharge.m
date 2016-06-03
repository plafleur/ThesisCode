%%Calculate the effective charge for the acoustic and optical mode in the one-dimensional chain

clear
n=101; %Number of q-points
q=1.602e-19;
e0=8.854e-12;
l=1.05e-15;
p=3.21e-11;
m1=3.81754035e-26;
m2=5.81e-26;
M=m1*m2/(m1+m2);%Reduced mass
mu1=sqrt(M/m1);%Mass coefficient in dynamical matrix
mu2=sqrt(M/m2);%%Mass coefficient in dynamical matrix
c1=zeros(1,n);
c2=zeros(1,n);
k=zeros(1,n);
Q1=1;%Using charge of +1 and -1 instead of +e and -e for simplicity
Q2=-1;
D0=72.2867019009821;
R0=3.192e-10;%Equilibrium separation in the infinite one-dimensional chain
for c=1:n;
	k(c)=(-pi+0.001)+(2*pi)*(c-1)/(n-1);%Over the range of -pi to pi
	Phi1=D0-(q^2/(2*pi*e0*(2*R0)^3))*ma1(k(c));%Calculate spring constants
	Phi2=[(q^2/(2*pi*e0*(2*R0)^3))*ma2(k(c))-2*(l/p^2)*exp(-(R0/p))*cos(k(c)/2)];
	D11=mu1^2*Phi1;%Build the dynamical matrix
	D12=mu1*mu2*Phi2;
	D22=mu2^2*Phi1;
	A=[D11 D12; D12 D22];
	[b,lambda]=eig(A);%Find eigenvalues and eigenvectors of the dynamical matris
	c1(c)=Q1*mu1*b(1)+Q2*mu2*b(3);%Calculate the effective charge for each mode
	c2(c)=Q1*mu1*b(3)+Q2*mu2*b(4);

end
