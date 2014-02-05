\name{shad.cjs}
\alias{shad.cjs}
\title{
analysis of the shad data provided in the scrbook package using JAGS
}
\description{
Function to analyze the American shad data under the non-spatial Cormack-Jolly-Seber model, the multi-state model, and the spatial Cormack-Jolly-Seber model.
}
\usage{
ovenbirds.js(model=c('SMS', 'NSJS', 'SJS'), n.chains, n.adapt, n.iter)
}
\details{
shad.js allows you to run a multi-state (2 state) model ('MS-CJS'), a non-spatial Cormack-Jolly-Seber model ('NS-CJS'), and fully spatial capture-recapture
Cormack-Jolly-Seber model ('S-CJS') using a subset of the shad data described by Raabe et al. 2014.  
Model results are displayed and discussed in Ch16 of Royle et al. 2014. Spatial Capture-Recapture.
}
\value{
JAGS model output (for details, see the rjags manual, under
coda.samples())
}


\references{
Raabe, J. K., Gardner, B., & Hightower, J. E. (2014). A spatial capture-recapture model to estimate fish 
survival and location from linear continuous monitoring arrays. Canadian Journal of Fisheries and Aquatic Sciences, 71(1), 120-130.
}

\author{
Beth Gardner
}

\seealso{
%% ~~objects to See Also as \code{\link{Ch16shaddata}}, ~~~
}
\examples{

##runs S-CJS model with 3 chains, 400 adaptive iterations an 2000 iterations
##NOTE that S-CJS models can take a long time to run

out<-shad.cjs(model="S-CJS", 3, 400, 2000)

##look at summary output

summary(out)

Iterations = 401:2400
Thinning interval = 1 
Number of chains = 3 
Sample size per chain = 2000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

            Mean      SD  Naive SE Time-series SE
lam0[1]   5.6024 0.22165 0.0028615      0.0055856
lam0[2]   4.4599 0.15606 0.0020148      0.0042420
lam0[3]   1.9025 0.06777 0.0008750      0.0013330
lam0[4]   1.1306 0.05479 0.0007073      0.0009147
lam0[5]   0.9513 0.05765 0.0007443      0.0009662
lam0[6]   0.3591 0.04075 0.0005261      0.0006822
lam0[7]   0.1881 0.03104 0.0004008      0.0004843
lam0[8]   0.3100 0.04405 0.0005687      0.0007690
lam0[9]   0.3628 0.05255 0.0006785      0.0007999
lam0[10]  0.6306 0.07241 0.0009349      0.0011029
lam0[11]  1.6266 0.10882 0.0014048      0.0017272
lam0[12]  0.9436 0.14016 0.0018094      0.0030365
phi       0.7844 0.01200 0.0001549      0.0002189
sigma    13.8424 0.19784 0.0025541      0.0066899

2. Quantiles for each variable:

            2.5%     25%     50%     75%   97.5%
lam0[1]   5.1769  5.4541  5.5950  5.7483  6.0489
lam0[2]   4.1543  4.3545  4.4569  4.5641  4.7700
lam0[3]   1.7765  1.8558  1.9009  1.9471  2.0400
lam0[4]   1.0237  1.0943  1.1298  1.1669  1.2424
lam0[5]   0.8411  0.9118  0.9504  0.9906  1.0662
lam0[6]   0.2852  0.3312  0.3575  0.3844  0.4457
lam0[7]   0.1321  0.1661  0.1862  0.2083  0.2529
lam0[8]   0.2297  0.2805  0.3075  0.3371  0.4026
lam0[9]   0.2686  0.3257  0.3595  0.3979  0.4734
lam0[10]  0.4983  0.5807  0.6282  0.6773  0.7800
lam0[11]  1.4229  1.5548  1.6225  1.6963  1.8510
lam0[12]  0.6921  0.8463  0.9364  1.0301  1.2503
phi       0.7600  0.7763  0.7847  0.7927  0.8070
sigma    13.4624 13.7073 13.8384 13.9736 14.2431




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ JAGS }
\keyword{ open models }% __ONLY ONE__ keyword per line
\keyword{ Cormack-Jolly-Seber }% __ONLY ONE__ keyword per line
\keyword{ shad data }% __ONLY ONE__ keyword per line