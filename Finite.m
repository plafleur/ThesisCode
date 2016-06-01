%#!/usr/local/bin/octave -qf
clear;
fid=fopen('C:\Users\Patrick\Dropbox\Thesis\Code\100atom.dat','w')
n=2e5;
at=100;
x=zeros(n,at);
v=zeros(n,at);
a=zeros(n,at);
t=zeros(n,1);
m=zeros(1,at);


for ma=2:2:at;
     m(ma)=5.81e-26;
end
for ma=1:2:at-1;
    m(ma)=3.818e-26;
end

for ab=1:at;
      x(1,ab)=(ab-1)*3.0437e-10;
      x(2,ab)=(ab-1)*3.0437e-10;
end

dt=1e-16;
qe=1.603e-19;
q=zeros(1,at);



for la=1:at;
    q(la)=((-1)^la)*qe;
end
    
con=1/(4*(3.14159)*8.854e-12);
l=1.05e-15;
p=3.21e-11;
lp=l/p;
gam=1e12;
atot=0;



   for k=2:n-1;
      
    
        for i=2:at;
    
           
             for j=[1:(i-1),(i+1):at];
             b=abs(i-j);
            
            if b==1
                 Ftot=con*((q(i)*q(j))/(abs(x(k,i)-x(k,j)))^2)+lp*exp(-(abs(x(k,i)-x(k,j)))/p);
                
            else
                 Ftot=con*((q(i)*q(j))/(abs(x(k,i)-x(k,j)))^2);                                
            end
            atot=atot+(1/m(i))*(i-j)/(abs(i-j))*Ftot;
           
            
                 
            end
            a(k,i)=atot-v(k,i)*gam;     
            Ftot=0;
            atot=0;
            
            
        end
    for i=1:at;      
    x(k+1,i)=x(k,i)+v(k,i)*(dt)+(1/6)*(4*a(k,i)-a(k-1,i))*(dt)^2;
    v(k+1,i)=v(k,i)+(1/6)*(5*a(k,i)-a(k-1,i))*dt;
    end
    t(k)=dt*k;
   end

       
       


    
   


pairs=zeros(1,at-1);
near=zeros(1,at-1);
pairs=[1:at-1];
for ac=1:at-1;
     near(ac)=x(k,ac+1)-x(k,ac);
end            
plot(pairs,near)
for b=1:at;
    fprintf(fid,'%e \n',x(k,b));
end
fclose(fid);