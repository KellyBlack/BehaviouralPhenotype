function tot=ButterflyWaspDiffusionModel(t,x,params)
 m=params(1);
 c=params(2);
 d=params(3);
 g=params(4);
 N=params(5);
 mu=params(6);

bdot(1)=(1+m/N)*x(1)*(1-x(1))-x(N+1)*((1+m/N)*x(1))/(c+x(1)*(1+m/N))+mu*(x(2)-x(1));
for i=2:N-1    
bdot(i)=(1+i*m/N)*x(i)*(1-x(i))-x(N+1)*((1+i*m/N)*x(i))/(c+x(i)*(1+i*m/N))+mu*(x(i+1)-2*x(i)+x(i-1));
end
bdot(N)=(1+N*m/N)*x(N)*(1-x(N))-x(N+1)*((1+N*m/N)*x(N))/(c+x(N)*(1+N*m/N))+mu*(-x(N)+x(N-1));
butterfly=0;
for i=1:N
 butterfly=butterfly+((1+i*m/N).*x(i))/(c+x(i)*(1+i*m/N));
end
wdot=-d*x(N+1)+(g/N)*x(N+1)*butterfly;
tot=[bdot';wdot];