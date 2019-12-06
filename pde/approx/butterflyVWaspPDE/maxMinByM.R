a <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingM-mu-0.2.csv')
b <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingM-mu-0.02.csv')
c <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingM-mu-0.1.csv')

highest <- max(c(a$maxButterfly,b$maxButterfly,c$maxButterfly))

#pdf("/tmp/behaviourPhenotypeSteadyState.pdf")
plot(b$m,b$minButterfly,type='l',
     ylim=c(0,1.05*highest),
     main="Max and Min Of The Butterfly Density After A Long Time Span",
     xlab='m',
     ylab='Butterfly Density',
     col=2,lwd=2)
points(b$m,b$maxButterfly,type='l',col=2,lwd=2)
points(a$m,a$maxButterfly,type='l',col=4,lwd=2)
points(a$m,a$minButterfly,type='l',col=4,lwd=2)
points(c$m,c$maxButterfly,type='l',col=6,lwd=2)
points(c$m,c$minButterfly,type='l',col=6,lwd=2)

legend(0.2,1.1,
       c(expression(paste(mu,'=0.2')),
         expression(paste(mu,'=0.1')),
         expression(paste(mu,'=0.02'))),
       lty=c(1,1,1),
       col=c(4,6,2),
       lwd=c(2,2,2)
       )
#dev.off()
