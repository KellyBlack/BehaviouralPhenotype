ode <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/rk45_c2.1.csv')
#ode <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/rk45.csv')
odeOrder <- sort(ode$m,index.return=TRUE)
cCoef = min(ode$c)

plot.new()
plot.window(xlim=c(0,1.05*max(ode$m)),
            ylim=c(0,1.16*max(ode$maxButterfly)))
title(main=
        paste("Max and Min Of The Butterfly Density\nAfter A Long Time Span, c=",cCoef),
      xlab='m',
      ylab='Butterfly Density')
axis(1)
axis(2)

##pdf("/tmp/behaviourPhenotypeSteadyState.pdf")
closeValues <- abs(ode$maxButterfly[odeOrder$ix]-ode$minButterfly[odeOrder$ix])< 3.71e-4
points(ode$m[odeOrder$ix],ode$minButterfly[odeOrder$ix],type='l',lty=1,lwd=3,col=4)
points(ode$m[odeOrder$ix],ode$maxButterfly[odeOrder$ix],type='l',lty=1,lwd=3,col=2)
points(ode$m[closeValues],ode$minButterfly[closeValues],type='l',lty=1,lwd=3,col=1)


#dev.off()
#dev.copy(pdf,'maxMinByM-mu-01-04.pdf')
#dev.off()
