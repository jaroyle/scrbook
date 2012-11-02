pGauss2 <-
function(parms,D){
  a0<-parms[1]
  sigma<-parms[2]
  lp<-  parms[1] -(1/(2*parms[2]*parms[2]))*D*D
  p<- 1-exp(-exp(lp))
  p
}
