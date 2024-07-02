plotResults <- function(a)
{
    currentColour <- 0
    muLevels  <- sort(unique(a$mu),decreasing=FALSE)
    colours   <- c()
    labels    <- c()
    plotTypes <- c()
    pchTypes  <- c()
    nextPCH   <- 0
    for(mu in muLevels)
    {
        nextPCH       <- nextPCH + 1
        currentColour <- currentColour + 1
        labels    <- c(labels,as.expression(bquote(mu == .(mu))))
        ##labels  <- c(labels,parse(text=paste0('mu',"=",mu)))
        colours   <- c(colours,currentColour)
        plotTypes <- c(plotTypes,0)
        pchTypes  <- c(pchTypes,nextPCH)

        m <- a$m[a$mu==mu]
        maxButterfly <- a$maxButterfly[a$mu==mu]
        minButterfly <- a$minButterfly[a$mu==mu]
        m <- sort(m,index.return=TRUE)
        points(m$x,maxButterfly[m$ix],type='p',pch=nextPCH,col=currentColour,cex=0.75) # lwd=2,
        points(m$x,minButterfly[m$ix],type='p',pch=nextPCH,col=currentColour,cex=0.75) # lwd=2
    }

    return(list(labels=labels,
                plotTypes=plotTypes,
                pchTypes=pchTypes,
                colours=colours))
    
}

filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingM-multipleMu.csv'
filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMResults_1.csv'
#filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMResults_c=1.1.csv'
#filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMHysteresisReverse.csv'

filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMResults_c=2.5.csv'
filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMResults_c=2.6.csv'
filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMResults_c=2.7500.csv'
filename <- '../build-butterflyVWaspPDE-Desktop-Debug/changingMResults_c=2.7500_smallMu.csv'
filename <- '/tmp/changingMResults_c=2.0000.csv'
a <- read.csv(filename)

odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45_c1.1.csv'
odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45.csv'
odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45_c2.1.csv'

odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45_c2.5.csv'
odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45_c2.6.csv'
odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45_c-2.7500.csv'
odeFilename <- '../build-butterflyVWaspPDE-Desktop-Debug/rk45_c-2.7500.csv'
odeFilename <- '/tmp/rk45_c-2.0000.csv'
ode <- read.csv(odeFilename)

plot.new()
plot.window(xlim=c(0,1.05*max(a$m,ode$m)),
            ylim=c(0,1.05*max(ode$maxButterfly,a$maxButterfly)))
title(main="Butterfly Density Max/Min After A Long Time Span",
      xlab='m',
      ylab='Butterfly Density')
axis(1)
axis(2)
results <- plotResults(a)

## m <- a$m[a$mu==mu]
## tim <- a$time[a$mu==mu]
## points(m,tim,pch=3,col=3)


#plot(a$m[a$mu==0.02],a$time[a$mu==0.02])
#points(a$m[a$mu==0.04],a$time[a$mu==0.04],col=2)



#pdf("/tmp/behaviourPhenotypeSteadyState.pdf")
#plot(b$m,b$minButterfly,type='l',
#     ylim=c(0,1.05*highest),
#     main="Max and Min Of The Butterfly Density After A Long Time Span",
#     xlab='m',
#     ylab='Butterfly Density',
#     col=2,lwd=2)


odeOrder <- sort(ode$m,index.return=TRUE)
points(ode$m[odeOrder$ix],ode$minButterfly[odeOrder$ix],type='l',lty=1,lwd=3)
points(ode$m[odeOrder$ix],ode$maxButterfly[odeOrder$ix],type='l',lty=1,lwd=3)
results$labels <- c(results$labels,as.expression('ODE'))
results$plotTypes <- c(results$plotTypes,1)
results$pchTypes <- c(results$pchTypes,-1)
results$colours <- c(results$colours,1)

#legend(max(ode$m)*0.85,0.5,results$labels,lty=results$plotTypes,
#       pch=results$pchTypes,col=results$colours,lwd=2)
legend(30,.9,as.character(results$labels),lty=results$plotTypes,
       pch=results$pchTypes,col=results$colours,lwd=2)
#legend(0.2,1.1,
#       c(expression(paste(mu,'=0.2')),
#         expression(paste(mu,'=0.1')),
#         expression(paste(mu,'=0.02'))),
#       lty=c(1,1,1),
#       col=c(4,6,2),
#       lwd=c(2,2,2)
#       )
#dev.off()
#dev.copy(pdf,'maxMinByM-c-1.1-mu-01-04.pdf')
#dev.off()
