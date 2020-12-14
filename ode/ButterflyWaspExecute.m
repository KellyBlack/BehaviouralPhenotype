clear all
t=linspace(0,400,4001)';

%Define Parameters
m=12;
c=2.8;
d=0.1;
g=0.6; 
N=40;  %number of compartments

params  = [m,c,d,g,N]; %array of parameters to be passed into the ODE

%Create an array of initial conditions
y0=zeros(N+1,1);
for i=1:N
y0(i,:)=0.1;
end
y0(N+1,:)=0.9;   %initial condition for wasps


[t,y]=ode45(@(t,y0) ButterflyWaspModel(t,y0,params),t,y0) %run ODE

%Sum up butterflies across compartments
for i=1:length(t)
    total_but(i)=sum(y(i,[1:N]))';
end


%Plot wasp population over time
figure(1);
plot(t, y(:,N+1),'r-')
xlabel('time');
ylabel('wasp population');
title('Wasp Population')


%Plot total butterfly population over time
figure(2);
plot(t, total_but(1,:)/N,'b-')
xlabel('time');
ylabel('butterfly population');
title('Total Butterfly Population')



%Plot the last five butterfly compartments
C = {'k','g','y','c'}
figure(3);
plot(t, y(:,N),'b-')
xlabel('time');
ylabel('butterfly population');
title('Butterfly Population by Compartment')
hold on
for i=1:4
plot(t, y(:,N-i),C{i})
hold on
end
legend({'compartment N';'compartment N-1';'compartment N-2';'compartment N-3';'compartment N-4'})
