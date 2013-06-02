post.median <- qbeta(0.5, 20, 10) 
post.95ci <- qbeta(c(0.025, 0.975), 20, 10) 

plot(seq(0,1, by=0.0001), dbeta(seq(0,1, by=0.0001), 20,10), type='l', axes=F,
xlab='Detection probability', ylab='Probability density',cex.lab=1.6, lwd=2, mgp=c(2.75,1,0))
axis(side=1, lwd=2, lwd.ticks=2, cex.axis=1.6, xpd=TRUE)
axis(side=2, lwd=2, lwd.ticks=2, cex.axis=1.6, xpd=TRUE)
box(lwd=2)
axis(side=1, pos=c(0,0), labels =F, lwd=2, lwd.ticks=0, col="lightgrey")
arrows(x0=post.median, y0=0, y1= dbeta(post.median, 20,10), lwd=2, length = 0, lty=2)
arrows(x0=post.median, y0=0, y1= dbeta(post.median, 20,10), lwd=2, length = 0, lty=2)
arrows(x0=post.95ci[1], y0=0, y1= dbeta(post.95ci[1], 20,10), lwd=2, length = 0, lty=3)
arrows(x0=post.95ci[2], y0=0, y1= dbeta(post.95ci[2], 20,10), lwd=2, length = 0, lty=3)