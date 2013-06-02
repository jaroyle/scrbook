
\begin{center}
\begin{verbatim}
  pixel ID                 Cost
 4 8 12 16            100   1   1  1
 3 7 11 15            100 100   1  1
 2 6 10 14            100 100 100  1
 1 5  9 13            100 100   1  1
\end{verbatim}
\end{center}

# THESE CALCS NEED UPDATED

Then we assigned low cost of 1 to ``good habitat'' pixels (or pixels
we think of as ``highly connected'' by virtue of being in good
habitat) and, conversely, we assign high cost (100) to ``bad
habitat''. So the shortest cost-weighted distance between pixels 5 and
9 in this example is just 1 unit, the shortest cost-distance between
pixels 5 and 10 is $\sqrt{2}(1+1)/2 = 1.414214$ units, the shortest
distance between pixels 4 and 8 is 100 units, while the shortest
cost-distance between 4 and 12 is 150.5. A tough one is: what is the
shortest distance between 7 and 16? An individual at pixel 7 can move
diagonal and pay $sqrt(2)*(100+1)/2 + 1 =72.41778$.  This simple cost
raster is shown in Fig. \ref{ecoldist.fig.raster}.


library(raster)
r<-raster(nrows=4,ncols=4)
projection(r)<- "+proj=utm +zone=12 +datum=WGS84" #sets the projection
#We use UTM here because distances and directions are correct
#i.e., they are adjusted for the earth's curvature
extent(r)<-c(.5,4.5,.5,4.5) #sets the extent of the raster
costs1<- c(100,100,100,100,1,100,100,100,1,1,100,1,1,1,1,1)
values(r)<-matrix(costs1,4,4,byrow=FALSE) #assign the costs to the raster



Once the cost raster is created, the
least-cost path distances are computed with just a couple {\bf R}
commands, and those can be inserted directly into the likelihood
construction for an ordinary spatial capture-recapture model (Appendix
1). The {\bf R} package \mbox{\tt gdistance} uses the implementation
of Dijkstra's algorithm \citep{dijkstra:1959} found in the \mbox{\tt
  igraph} package \citep{csardi:2010}.  Using \mbox{\tt gdistance},
we define the incremental cost of moving from one pixel
to another as the distance-weighted {\it average} of the 2 pixel
costs.

To compute the least-cost path, or the minimum cost-weighted
distances between every pixel and every other pixel, we make use of 
the helper functions \mbox{\tt transition}, which
calculates the cost of moving between neighboring pixels, and
\mbox{\tt geoCorrection} which modifies the costs of moving diagonally
by the additional distance, produce output which feeds into the
function \mbox{\tt costDistance} to compute the pair-wise distance
matrix. For that, we define the center points of each raster.  The
commands altogether are as follows:

# The transition function specifies that we'd like to use the mean of
# 2 neighboring pixels as the weight and then takes the inverse to convert
# the costs into conductances, which algorithms within gdistance and
# igraph require. We use the argument directions= 8 to specify that
# all 8 touching pixels are neighbors.)

library(gdistance)
tr1<-transition(r,transitionFunction=function(x) 1/mean(x),directions=8)

#The geoCorrection function corrects the conductances for the diagonal neighbors
#and, if the raster is in latitude longitude, also corrects for 
#curvature of the earth's surface

tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)

#here we specify the locations that we'd like to calculate distances among
# 
pts<-cbind( sort(rep(1:4,4)),rep(1:4,4))

#the costDistance function calculates the least cost distane among pts using
#the conductances specified in tri1CorrC
costs1<-costDistance(tr1CorrC,pts)
#here we convert costs1 into a matrix
outD<-as.matrix(costs1)




Now we can look at the result and see if it makes sense to us. Here we
print the first 4 columns of this distance matrix and illustrate a
couple of examples of calculating the minimum cost-weighted distance
between points:
\begin{verbatim}
> outD[1:5,1:5]
         1        2        3        4        5
1   0.0000 100.0000 200.0000 205.2426 100.0000
2 100.0000   0.0000 100.0000 200.0000 141.4214
3 200.0000 100.0000   0.0000 100.0000 126.1604
4 205.2426 200.0000 100.0000   0.0000 105.2426
5 100.0000 141.4214 126.1604 105.2426   0.0000
\end{verbatim}


### Not in book chapter
So, to calculate the distance from the first cell in the lower left 
corner to the cell above it, the least cost distance is just the 
mean of the two cells because that is both the shortest line and 
all cells surrounding the first the cost in the first cell is 100 and
the cost in the cell immediately above it is also 100.  So the distance
is (100+100)/2 = 100, which is the same as the mean of the 2 pixels.
%% NOTE: This uses Tabitha's number of pixels
\begin{verbatim}
 plot(r)
 points(pts[1,1],pts[1,2],col="red")
 points(pts[2,1],pts[2,2],col="blue")
 lines(c(1,1),c(1,2))
\end{verbatim}
To move from the pixel in the lower left corner to the upper left corner
points 1 \& 4, the shortest distance is to go directly up.
\begin{verbatim}
 points(pts[4,1],pts[4,2],col="blue")
 lines(c(1,1),c(1,4))
#but the least cost distance is to go around the outside of the high cost area.
 lines(c(1,2),c(1,1)) #with cost = (100 + 100)/2 = 100
 lines(c(2,3),c(1,1)) #with cost = (100 + 1)/2 = 50.5
#To account for the increased distance along the diagonal, we multiply by sqrt(2)
 lines(c(3,4),c(1,2)) #with cost = sqrt(2) * (1+1)/2 = 1.414214
 lines(c(4,3),c(2,3)) #with cost = sqrt(2) * (1+1)/2 = 1.414214
 lines(c(3,2),c(3,4)) #with cost = sqrt(2) * (1+1)/2 = 1.414214
 lines(c(2,1),c(4,4)) #with cost = (100 + 1)/2 = 50.5
 100+50.5* 2+(1.414214)* 3 # = 205.2426 
\end{verbatim}
This matches the distance in our distance matrix between points 1 and 4. 


