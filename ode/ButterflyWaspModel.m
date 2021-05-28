function tot=ButterflyWaspModel(t,x,params)
 m=params(1);
 c=params(2);
 d=params(3);
 g=params(4);
 N=params(5);
 
butterfly=0;
for i=1:N    
bdot(i)=(1+i*m/N)*x(i)*(1-x(i))-x(N+1)*((1+i*m/N)*x(i))/(c+x(i)*(1+i*m/N));
butterfly=butterfly+((1+i*m/N).*x(i))/(c+x(i)*(1+i*m/N));
end
wdot=-d*x(N+1)+(g/N)*x(N+1)*butterfly;
tot=[bdot';wdot];