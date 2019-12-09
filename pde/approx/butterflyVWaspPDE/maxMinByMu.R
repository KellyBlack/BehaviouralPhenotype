#a <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingM-multipleMu.csv')
a <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingMResults.csv')

plot.new()
plot.window(xlim=c(0,1.05*max(a$mu)),
            ylim=c(0,1.06*max(a$maxButterfly)))
title(main="Max and Min Of The Butterfly Density After A Long Time Span",
      xlab='mu',
      ylab='Butterfly Density')
axis(1)
axis(2)
currentColour <- 0
mLevels  <- sort(unique(a$m),decreasing=TRUE)
colours   <- c()
labels    <- c()
plotTypes <- c()
pchTypes  <- c()
nextPCH   <- 0
for(m in mLevels)
{
    nextPCH       <- nextPCH + 1
    currentColour <- currentColour + 1
    labels    <- c(labels,as.expression(bquote(m == .(m))))
    ##labels  <- c(labels,parse(text=paste0('mu',"=",mu)))
    colours   <- c(colours,currentColour)
    plotTypes <- c(plotTypes,0)
    pchTypes  <- c(pchTypes,nextPCH)

    mu <- a$mu[a$m==m]
    maxButterfly <- a$maxButterfly[a$m==m]
    minButterfly <- a$minButterfly[a$m==m]
    mu <- sort(mu,index.return=TRUE)
    points(mu$x,maxButterfly[mu$ix],type='p',pch=nextPCH,col=currentColour,cex=0.75) # lwd=2,
    points(mu$x,minButterfly[mu$ix],type='p',pch=nextPCH,col=currentColour,cex=0.75) # lwd=2
}

## mu <- a$mu[a$m==m]
## tim <- a$time[a$m==m]
## points(mu,tim,pch=3,col=3)


#plot(a$mu[a$m==0.02],a$time[a$m==0.02])
#points(a$mu[a$m==0.04],a$time[a$m==0.04],col=2)



#pdf("/tmp/behaviourPhenotypeSteadyState.pdf")
legend(0.0,1.35,labels,lty=plotTypes,pch=pchTypes,col=colours,lwd=2)
#dev.off()
