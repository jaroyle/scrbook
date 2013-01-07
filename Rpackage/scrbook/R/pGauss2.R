pGauss2 <-
function(parms,Dmat){
  a0<-parms[1]
  sigma<-parms[2]
  lp<-  parms[1] -(1/(2*parms[2]*parms[2]))*Dmat*Dmat
  p<- 1-exp(-exp(lp))
  p
}
