rot <-
function (m) 
{
# takes a matrix m and rotates it such that image(rot(m)) is in
# real space.
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}
