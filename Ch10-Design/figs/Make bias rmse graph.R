###make bias/precision graph for design chapter
plottab<-read.csv('c:/work/scrbook/Ch10-Design/figs/Out_Sims.csv')
plot(plottab[,1:2], type='l', ylim=c(-0.1, 0.5),ylab='Relative RMSE/Bias', xlab='Trap spacing [sigma]', lwd=2, cex.axis=1.5, cex.lab=1.5)
points(plottab[,c(1,3)], type='l', col="grey", lwd=2)
legend(x=0.5, y=0.5, legend=c("RMSE", "Bias"), col=c("black", "grey"), lwd=2, cex=1.5)