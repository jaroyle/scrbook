spatial.plot<-
function (x, y, add = FALSE, cx = 1,col="gray")
{
    nc <- as.numeric(cut(y, 10))
    if (!add)
        plot(x, pch = " ",asp=1)
        if(col=="gray"){
        cc<-seq(3,17,,10)/20
        cc<-gray(cc)
}
else
 cc<-terrain.colors(10)
  points(x, pch = 20, col = cc[nc], cex = cx)
  image.scale(y, col = cc)
}
