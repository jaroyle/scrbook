\name{secr.bear}
\alias{secr.bear}
\title{
analysis of the Ft. Drum black bear data using secr
}
\description{
analysis of the Ft. Drum black bear data using a range of models in secr
}
\usage{
secr.bear()
}
\details{
secr.bear() runs the null model (bear.0), the null model with negative exponential detection function (bear.0exp),
time as factor on detection (bear.t), a global and a local trap response model (bear.B and bear.b, respectively), 
a time (factor) plus global behavioral response model (bear.Bt), a model with sex specific detection and movement (bear.sex) 
and a heterogeneity model with 2-class mixtures on detection and movement (bear.h2). The function further compares all models by their AIC.
}
\value{
List with AIC table and output for all models:
\item{AIC.tab}{AIC table}
\item{bear.0}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~1, g0~1, sigma~1),buffer = 20000)}
\item{bear.0exp}{secr output of model bear.0exp=secr.fit (bear.cap, model=list(D~1, g0~1, sigma~1),buffer = 20000,detectfn=2)}
\item{bear.t}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~1, g0~t, sigma~1),buffer = 20000)}
\item{bear.B}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~1, g0~b, sigma~1),buffer = 20000)}
\item{bear.b}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~1, g0~bk, sigma~1),buffer = 20000)} 
\item{bear.Bt}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~1, g0~b+t, sigma~1),buffer = 20000)}
\item{bear.sex}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~session, g0~session, sigma~session),buffer = 20000)}
\item{bear.h}{secr output of model bear.0=secr.fit (bear.cap, model=list(D~1, g0~h2, sigma~h2),buffer = 20000)}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rahel Sollmann, Beth Gardner
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

bear.ou<-secr.bear()
secr.bear$AIC.tab

## returns:
                                       model    detectfn npar    logLik
bear.b                     D~1 g0~bk sigma~1  halfnormal    4 -641.7215
bear.h2           D~1 g0~h2 sigma~h2 pmix~h2  halfnormal    6 -653.8382
bear.0exp                   D~1 g0~1 sigma~1 exponential    3 -663.9152
bear.B                      D~1 g0~b sigma~1  halfnormal    4 -677.6175
bear.Bt                 D~1 g0~b + t sigma~1  halfnormal   11 -668.3044
bear.sex  D~session g0~session sigma~session  halfnormal    6 -677.7151
bear.t                      D~1 g0~t sigma~1  halfnormal   10 -674.4134
bear.0                      D~1 g0~1 sigma~1  halfnormal    3 -686.2455
               AIC     AICc  dAICc AICwt
bear.b    1291.443 1292.395  0.000     1
bear.h2   1319.676 1321.776 29.381     0
bear.0exp 1333.830 1334.389 41.994     0
bear.B    1363.235 1364.187 71.792     0
bear.Bt   1358.609 1366.152 73.757     0
bear.sex  1367.430 1369.530 77.135     0
bear.t    1368.827 1374.938 82.543     0
bear.0    1378.491 1379.049 86.654     0

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ secr }
\keyword{ AIC }% __ONLY ONE__ keyword per line
