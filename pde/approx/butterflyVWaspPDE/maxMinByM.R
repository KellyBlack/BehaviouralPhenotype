#a <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingM-multipleMu.csv')
a <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingMResults.csv')

plot.new()
plot.window(xlim=c(0,1.05*max(a$m)),
            ylim=c(0,1.06*max(a$maxButterfly)))
title(main="Max and Min Of The Butterfly Density After A Long Time Span",
      xlab='m',
      ylab='Butterfly Density')
axis(1)
axis(2)
currentColour <- 0
muLevels  <- sort(unique(a$mu),decreasing=TRUE)
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

legend(0.1,0.05,labels,lty=plotTypes,pch=pchTypes,col=colours,lwd=2)
#legend(0.2,1.1,
#       c(expression(paste(mu,'=0.2')),
#         expression(paste(mu,'=0.1')),
#         expression(paste(mu,'=0.02'))),
#       lty=c(1,1,1),
#       col=c(4,6,2),
#       lwd=c(2,2,2)
#       )
#dev.off()
