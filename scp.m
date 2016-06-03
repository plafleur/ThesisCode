%%This code calculates the self-consistent phonon dispersion for the soft mode of strontium titanate
%%as a function of wavevector,electric field, and temperature. This requires a three-dimensional sum 
%%over wavevectors to determine the temperature fluctuation term. The q=0 value gives the differential
%%susceptibility

clear;
n=100; %Using 100 q-points in each of the three dimensions
q=linspace(0,pi,n);%Define q-points between 0 and pi (using the symmetry of the systems)
cosq=cos(q);%Define an array that contains the value of the cosines for each q-point
d=1;%Number of temperatures to be used

%100 q-point parameters
C=3;
A=21;
B=3.1e23;
Q=8*1.602e-19;

V=(3.905e-10)^3;%Volume of the unit cell
e0=8.854e-12;

M=4e-26;%Effective mass of the soft mode
dx2old=(abs(A-12*C)+0.5)/B;%Starting point for the self-consistent calculation (Initial guess)
dx2new=1e10;
diff=abs(dx2old-dx2new);
wsum=0;
E=0;
%Define constans rather than referencing individual values
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
    T(k)=300; %For one temperature, can set here. Or can set T(k)=k for a range of temperatures
    clc
    disp(k)%Displays which temperature is being calculated
	while diff>1e-27;%Criteria for self-consistency, comparing old fluctuation value to new fluctuation value
  		for x=1:n; %Three dimensional sum, with the qz value represented as an array to save computation time
		  	for y=1:n;
           			  F=cosq(x)+cosq(y)+cosq(:);     
			  		  w=sqrt(1/M)*sqrt(A-4*F*C+B*dx2old);
                      wsum=wsum+sum((1./w).*coth((hbar2kb.*w)./T(k)));%Summing over frequencies
	 	 end
  	end
    wq0=(A-12*C+B*dx2old);
  	dx2new=con2*wsum+(1/3)*(con3/wq0)^2;
    diff=abs(dx2old-dx2new);
  	dx2old=0.7*dx2old+0.3*dx2new;  %Calculate new fluctuation value, using more complicated mixing procedure
   	wsum=0;
    end
  omega0=(A-12*C+B*dx2old);
  %omega=(A-4.*C.*(cosq+cos(0)+cos(0))+B*dx2old);
  chi(k)=con./(omega0+(2*B/3)*(Q*E/omega0)^2);%Calculate the differential susceptibility
  diff=1;
end

