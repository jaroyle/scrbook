area <-
function(x, y)
{
        if(missing(y)) {
                if(is.matrix(x) && ncol(x) == 2.) {
                        y <- x[, 2.]
                        x <- x[, 1.]
                }
                else if(!is.null(x$x) && !is.null(x$y)) {
                        y <- x$y
                        x <- x$x
                }
        }
        x <- c(x, x[1.])
        y <- c(y, y[1.])
        i <- 2.:length(x)
        return(0.5 * sum(x[i] * y[i - 1.] - x[i - 1.] * y[i]))
}
