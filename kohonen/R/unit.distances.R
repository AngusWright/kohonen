### Calculate distances in a Kohonen map. Crude and
### slow implementation, but hey.

unit.distances <- function(grid, toroidal)
{
  if (missing(toroidal)) toroidal <- grid$toroidal

  if (!toroidal) {
    if (grid$topo == "hexagonal") {
      return(as.matrix(stats::dist(grid$pts)))
    } else {
      return(as.matrix(stats::dist(grid$pts, method="maximum")))
    }
  }

  ## only for toroidal maps:
  np <- nrow(grid$pts)
  maxdiffx <- grid$xdim/2
  maxdiffy <- max(grid$pts[,2])/2
  
  result <- matrix(0, np, np)
  for (i in 1:(np-1)) {
    for (j in (i+1):np) {
      diffs <- abs(grid$pts[j,] - grid$pts[i,])
      if (diffs[1] > maxdiffx)
        diffs[1] <- 2*maxdiffx - diffs[1]
      if (diffs[2] > maxdiffy)
        diffs[2] <- 2*maxdiffy - diffs[2]
      
        if (grid$topo == "hexagonal") {
          result[i,j] <- sum(diffs^2)
        } else {
          result[i,j] <- max(diffs)
        }
    }
  }

  if (grid$topo == "hexagonal") {
    sqrt(result + t(result))
  } else {
    result + t(result)
  }
}

unit.distances.fast<-function (grid, toroidal) {
    if (missing(toroidal)) 
        toroidal <- grid$toroidal
    if (!toroidal) {
        if (grid$topo == "hexagonal") {
            return(as.matrix(stats::dist(grid$pts)))
        } else {
            return(as.matrix(stats::dist(grid$pts, method = "maximum")))
        }
    }
    dx<-stats::dist(grid$pts[,1],diag=TRUE,upper=TRUE)
    dy<-stats::dist(grid$pts[,2],diag=TRUE,upper=TRUE)
    max.dx<-max(dx)
    max.dy<-max(dy)
    dx.diff<-max.dx-dx+1
    dy.diff<-max.dy-dy+1
    dx<-pmin(dx,dx.diff)
    dy<-pmin(dy,dy.diff)
    if (grid$topo == "hexagonal") {
        result <- sqrt(dx^2+dy^2)
    } else {
        result <- pmax(dx,dy)
    }
    return(as.matrix(result))
}
