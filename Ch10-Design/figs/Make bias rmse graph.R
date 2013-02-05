###make bias/precision graph for design chapter
plottab<-read.csv('Out_Sims.csv')
plot(plottab[,1:2], type='l', ylim=c(-0.1, 0.5),ylab='Relative RMSE/Bias', xlab='trap spacing [sigma]', lwd=2, cex.axis=1.3, cex.lab=1.3)
points(plottab[,c(1,3)], type='l', col="red", lwd=2)
legend(x=0.5, y=0.5, legend=c("RMSE", "Bias"), col=c("black", "red"), lwd=2, cex=1.3)