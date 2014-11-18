\name{intlik3}
\alias{intlik3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Computes marginal likelihood for model SCR0 using either the Binomial
  or Poisson observation models.  
  %%  ~~function to do ... ~~
}
\description{
Use this to obtain MLEs for the SCR model based on binomial encounter
frequencies (out of K samples) or based on the Poisson encounter
frequency model. The input data here are the 2-d 
array of individual x trap encounter frequencies, therefore the data
must be summarized over the "occasions". 
}
\usage{
intlik3(start = NULL, y = y, K = NULL, delta = 0.3, X = traplocs, ssbuffer = 2,predict=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{
    Starting values of model parameters, a vector of length 3. The
    paramerization is (alpha0, log(alpha1), log(n0) ) where probcap =
    expit(alpha0)* exp(-alpha1*D*D).  alpha1 = (1/(2*sigma*sigma)) in
    the notation of the Gaussian kernel.
    Population size, N = n + n0. 
}
  \item{y}{
    "nindividuals x ntraps" encounter frequency matrix.
    y[i,j] = number of times individual i was encountered in trap j.
}
  \item{K}{
Number of sampling occasions. 
}
  \item{delta}{
Grid spacing for the numerical integration grid, in units of the trap
coordinate system (see next argument).
}
  \item{X}{
Trap coordinates, should be "ntraps x 2" matrix. 
}
  \item{ssbuffer}{
Width of the buffer around the minimum area rectangle enclosing the
traps. Should be in units of the trap coordinate system. E.g., if trap
locations are in UTM coordinates, meters, and you specify ssbuffer=1,
this is a 1 meter buffer. Normally this would not be satisfactory. 
%%     ~~Describe \code{ssbuffer} here~~
}
\item{model}{
  "B" = binomial encounter model,
  "P" = poisson encounter model
}
\item{predict}{
If "TRUE" then the posterior probability surface of each activity center
(for the captured individuals) is returned for the given parameter
values. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
SCR book XXXX%% ~put references to the literature/web site here ~
}
\author{
Andy Royle, aroyle@usgs.gov 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library("scrbook")
data("wolverine")

traps<-wolverine$wtraps
traplocs<-traps[,2:3]/10000
K.wolv<-apply(traps[,4:ncol(traps)],1,sum)
y3d<-SCR23darray(wolverine$wcaps,wolverine$wtraps)
y2d<-apply(y3d,c(1,2),sum)


starts<-c(-1.5,0,3)
frog<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.2,ssbuffer=2)

#
#test out different integration grid densities. This takes a long time
#
#unix.time(a1<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.3,ssbuffer=2)$estimate)
#unix.time(a2<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.2,ssbuffer=2)$estimate)
#unix.time(a3<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.1,ssbuffer=2)$estimate)
#unix.time(a4<-nlm(intlik3,starts,hessian=TRUE,y=y2d,K=K.wolv,X=traplocs,delta=.05,ssbuffer=2)$estimate)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
