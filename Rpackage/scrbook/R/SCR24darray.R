SCR24darray <-
function (edf, tdf)
{


# Check for dups
    uq<- paste(edf[,2],edf[,3],edf[,4])
    uq<- unique(uq)
    if(any(uq>1)) cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used",fill=TRUE)


    nind <- max(edf[, 2])
    ntraps <- nrow(tdf)
    nperiods <- ncol(tdf) - 3

    # check no individual captured in occasion > nperiods
    if(max(edf[,3])> nperiods) {
        cat("edf[,3] indicates more sample periods than trap operation columns",fill=TRUE)
        return(NULL)
        }

#individual ID should be sequential integer
ind.id<- edf[,2]
if(length(unique(ind.id)) != max(ind.id)){
 cat("individual ID needs to be numbered sequentially, converting...,",fill=TRUE)
 individual.names<- levels(as.factor(edf[,2]))
 edf[,2]<-  as.numeric(as.factor(edf[,2]))
}
else{
    individual.names<- 1:nind
    }

    session<- edf[,1]
    u.session<- unique(session)
    u.session<- u.session[order(u.session)]
    n.session<-length(u.session)
    session.id<- match(session,u.session)


# Trap ID should be sequential integer
trap.id<- as.numeric(tdf[,1])
if(any(is.na(trap.id))){
    cat("some trap ID are not numeric, converting all to sequential integers",fill=TRUE)
    real.trap.id<- match(edf[,4],tdf[,1])
    edf[,4] <- real.trap.id
    tdf[,1] <- 1:nrow(tdf)
 }


    per.id <- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
    if(any(is.na(per.id))){
        cat("some period ID's not numeric. Coercing to sequential integers",fill=TRUE)
        per.id<- 1:nperiods
        }

    ind.id <- edf[, 2]
    trap.id <- edf[, 4]
    if (length(per.id) != length(min(per.id):max(per.id))) {
        x <- 1:nperiods
        names(x) <- as.character(per.id)
        per.id <- x[as.character(edf[, 3])]
    }
    else {
        per.id <- edf[, 3]
    }
    y <- array(0, c(n.session,nind, ntraps, nperiods))
    tmp <- cbind(session.id, ind.id, trap.id, per.id)
    y[tmp] <- 1
    attr(y, "session.names")<- u.session
    attr(y, "individual.names")<- individual.names
    y
}
