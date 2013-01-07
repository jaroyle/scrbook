pGauss1 <-
function(parms,Dmat){
  a0<-parms[1]
  sigma<-parms[2]
  p<-  plogis(parms[1])*exp( -(1/(2*parms[2]*parms[2]))*Dmat*Dmat )
  p
}
