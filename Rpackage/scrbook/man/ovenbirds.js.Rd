\name{ovenbirds.js}
\alias{ovenbirds.js}
\title{
analysis of the ovenbird data provided in the secr package using JAGS
}
\description{
Function to analyze  the ovenbird data under the spatial multi-session model, the non-spatial Jolly-Seber model, and the spatial Jolly-Seber model.
}
\usage{
ovenbirds.js(model=c('SMS', 'NSJS', 'SJS'), n.chains, n.adapt, n.iter)
}
\details{
ovenbirds.js allows you to run a spatial multi-session model ('SMS'), a non-spatial Jolly-Seber model ('NSJS'), and fully spatial capture-recapture
Jolly-Seber model ('SJS') using the overbird data provided in the R package secr.  
Model results are displayed and discussed in Ch16 of Royle et al. 2014. Spatial Capture-Recapture.
}
\value{
JAGS model output (for details, see the rjags manual, under
coda.samples())
}

\author{
Beth Gardner
}

\seealso{
%% ~~objects to See Also as \code{\link{ovenbird}}, ~~~
}
\examples{

##runs SJS model with 3 chains, 300 adaptive iterations an 1000 iterations
##NOTE that the SMS and SJS models can take a long time to run

out=ovenbirds.js(model="SJS", 3, 300, 1000)
##look at summary output
##Density given in animals per meter squared
##shown in book as animals per hectare

summary(out)
Iterations = 301:1300
Thinning interval = 1 
Number of chains = 3 
Sample size per chain = 1000 

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

               Mean        SD  Naive SE Time-series SE
D[1]      9.648e-05 1.461e-05 2.668e-07      6.313e-07
D[2]      1.023e-04 1.239e-05 2.262e-07      5.375e-07
D[3]      1.141e-04 1.372e-05 2.505e-07      6.540e-07
D[4]      1.071e-04 1.428e-05 2.607e-07      6.239e-07
D[5]      7.999e-05 1.272e-05 2.322e-07      5.906e-07
N[1]      3.349e+01 5.073e+00 9.261e-02      2.191e-01
N[2]      3.552e+01 4.300e+00 7.851e-02      1.866e-01
N[3]      3.961e+01 4.762e+00 8.695e-02      2.270e-01
N[4]      3.719e+01 4.956e+00 9.049e-02      2.166e-01
N[5]      2.777e+01 4.415e+00 8.060e-02      2.050e-01
R2        1.463e+01 4.852e+00 8.859e-02      2.348e-01
R3        1.895e+01 4.831e+00 8.821e-02      2.084e-01
R4        1.954e+01 4.849e+00 8.853e-02      2.075e-01
R5        8.428e+00 3.731e+00 6.812e-02      1.718e-01
alpha0   -2.901e+00 1.283e-01 2.342e-03      6.439e-03
alpha1    1.245e-04 1.433e-05 2.616e-07      1.005e-06
gamma[1]  5.001e-01 2.892e-01 5.279e-03      5.280e-03
gamma[2]  9.238e-02 3.563e-02 6.505e-04      1.679e-03
gamma[3]  1.302e-01 4.051e-02 7.397e-04      1.733e-03
gamma[4]  1.519e-01 4.654e-02 8.497e-04      1.990e-03
gamma[5]  8.154e-02 3.955e-02 7.222e-04      1.712e-03
phi       5.384e-01 5.914e-02 1.080e-03      2.542e-03
psi       1.706e-01 3.626e-02 6.619e-04      1.303e-03
sigma     6.368e+01 3.652e+00 6.667e-02      2.497e-01

2. Quantiles for each variable:

               2.5%        25%        50%        75%      97.5%
D[1]      7.202e-05  8.642e-05  9.507e-05  1.066e-04  0.0001268
D[2]      8.354e-05  9.218e-05  1.008e-04  1.095e-04  0.0001296
D[3]      9.218e-05  1.037e-04  1.123e-04  1.239e-04  0.0001440
D[4]      8.354e-05  9.795e-05  1.066e-04  1.152e-04  0.0001383
D[5]      5.762e-05  7.202e-05  7.778e-05  8.642e-05  0.0001095
N[1]      2.500e+01  3.000e+01  3.300e+01  3.700e+01 44.0250000
N[2]      2.900e+01  3.200e+01  3.500e+01  3.800e+01 45.0000000
N[3]      3.200e+01  3.600e+01  3.900e+01  4.300e+01 50.0000000
N[4]      2.900e+01  3.400e+01  3.700e+01  4.000e+01 48.0000000
N[5]      2.000e+01  2.500e+01  2.700e+01  3.000e+01 38.0000000
R2        5.000e+00  1.100e+01  1.500e+01  1.800e+01 25.0000000
R3        1.000e+01  1.600e+01  1.900e+01  2.200e+01 30.0000000
R4        1.000e+01  1.600e+01  2.000e+01  2.300e+01 29.0000000
R5        1.000e+00  6.000e+00  8.000e+00  1.100e+01 16.0000000
alpha0   -3.160e+00 -2.991e+00 -2.900e+00 -2.812e+00 -2.6486385
alpha1    9.796e-05  1.146e-04  1.235e-04  1.335e-04  0.0001551
gamma[1]  2.801e-02  2.462e-01  5.054e-01  7.471e-01  0.9732959
gamma[2]  3.280e-02  6.659e-02  8.943e-02  1.161e-01  0.1689490
gamma[3]  6.283e-02  1.012e-01  1.264e-01  1.557e-01  0.2179725
gamma[4]  6.713e-02  1.201e-01  1.494e-01  1.818e-01  0.2485741
gamma[5]  1.536e-02  5.240e-02  7.813e-02  1.054e-01  0.1673875
phi       4.230e-01  4.982e-01  5.372e-01  5.796e-01  0.6504095
psi       1.075e-01  1.444e-01  1.680e-01  1.943e-01  0.2481473
sigma     5.677e+01  6.119e+01  6.364e+01  6.604e+01 71.4432416



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ JAGS }
\keyword{ open models }% __ONLY ONE__ keyword per line
\keyword{ Jolly-Seber }% __ONLY ONE__ keyword per line
\keyword{ ovenbird data }% __ONLY ONE__ keyword per line