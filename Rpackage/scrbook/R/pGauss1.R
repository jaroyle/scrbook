pGauss1 <-
function(parms,D){
  a0<-parms[1]
  sigma<-parms[2]
  p<-  plogis(parms[1])*exp( -(1/(2*parms[2]*parms[2]))*D*D )
  p
}
