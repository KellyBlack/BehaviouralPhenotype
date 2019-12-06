a <- read.csv('../build-butterflyVWaspPDE-Desktop-Debug/changingM-multipleMu.csv')

plot.new()
plot.window(xlim=c(0,1.05*max(a$m)),
            ylim=c(0,1.05*max(a$maxButterfly)))
title(main="Max and Min Of The Butterfly Density After A Long Time Span",
      xlab='m',
      ylab='Butterfly Density')
axis(1)
axis(2)
currentColour <- 0
muLevels <- sort(unique(a$mu))
colours <- c()
labels <- c()
for(mu in muLevels)
{
    currentColour <- currentColour + 1
    labels <- c(labels,as.expression(bquote(mu == .(mu))))
    ##labels <- c(labels,parse(text=paste0('mu',"=",mu)))
    colours <- c(colours,currentColour)
    points(a$m[a$mu==mu],a$maxButterfly[a$mu==mu],
           type='l',col=currentColour,lwd=2)
    points(a$m[a$mu==mu],a$minButterfly[a$mu==mu],
           type='l',col=currentColour,lwd=2)
}



#pdf("/tmp/behaviourPhenotypeSteadyState.pdf")
#plot(b$m,b$minButterfly,type='l',
#     ylim=c(0,1.05*highest),
#     main="Max and Min Of The Butterfly Density After A Long Time Span",
#     xlab='m',
#     ylab='Butterfly Density',
#     col=2,lwd=2)

legend(0.1,1.15,labels,lty=1,col=colours,lwd=2)
#legend(0.2,1.1,
#       c(expression(paste(mu,'=0.2')),
#         expression(paste(mu,'=0.1')),
#         expression(paste(mu,'=0.02'))),
#       lty=c(1,1,1),
#       col=c(4,6,2),
#       lwd=c(2,2,2)
#       )
#dev.off()
