%This code integrates the PDE model using explicit finite difference method with forward differencing in time 
%and central differencing in space with neuman boundary conditions

function oscillation_status = PDE_finite_diff_flow(mu,m,c,d,g)

%mu=0.01;  
%m=5;
%c=2.8;
%d=0.1;
%g=0.6; 

tmax=100;
T=140000;  %number intervals in time
N=100; %number intervals in space
L=1;

dx = L/N;
x = 0:dx:L-dx; % define the mesh in space
y = 0:dx:L-dx; 
dt = tmax/T;
t = 0:dt:tmax; % define the mesh in time


D = mu*dt/dx^2; %diffusion rate should be less than 1/2 
B=zeros(length(x),1);
W=zeros(length(t),1);
B(:,1)=0.1; %initial conditions for butterfly
W(1,1)=0.9; %initial conditions for wasp

%Parameters
p=1+m.*x';

for j=1:T %time stepping
	  %left boundary point
  B(1,j+1) = (1-2*D).*B(1,j)+2*D.*B(2,j)+dt.*(p(1).*B(1,j).*(1-B(1,j))-W(j)*p(1).*B(1,j)./(c+p(1).*B(1,j)));
    
  for i=2:N-1
				%interior
    B(i,j+1) = (1-2*D).*B(i,j)+D.*B(i+1,j)+D*B(i-1,j)+dt.*(p(i).*B(i,j).*(1-B(i,j))-W(j).*p(i).*B(i,j)./(c+p(i).*B(i,j)));
    
  end
				%right boundary point
  B(N,j+1) = (1-2*D).*B(N,j)+2*D.*B(N-1,j)+dt.*(p(N).*B(N,j).*(1-B(N,j))-W(j).*p(N).*B(N,j)./(c+p(N)*B(N,j)));
        
  
  butterfly=g*W(j).*p.*B(:,j)./(c+p.*B(:,j));
  Int_w = simps(x,butterfly);
				%Int_w = trapz(x,butterfly);
  W(j+1)=W(j)+dt*(-d*W(j)+Int_w); %using Euler's method

end

osc1 = max(B(1,100000:140000))/min(B(1,100000:140000));
oscN = max(B(N,100000:140000))/min(B(N,100000:140000));
oscillation_status=1;
if ( osc1 < 1.2 & oscN < 1.2 )
    oscillation_status=-1;
end 

figure;
subplot(2,1,1); plot(t,B(1,:)); xlabel('time'); ylabel('B(1)')
title(['\mu=' num2str(mu) ', m=' num2str(m) ', c=' num2str(c) ', d=' num2str(d) ', g=' num2str(g)])
subplot(2,1,2); plot(t,B(100,:)); xlabel('time'); ylabel('B(100)')


