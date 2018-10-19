#' @title Conditional and joint multivariate quantiles
#'
#' @description
#' Quantiles of observations in multivariate space.
#'
#' @param x vector, matrix, or data.frame of observations.
#'
#' @param pltype character, plot type, one of
#'      \code{c('none','cdf','pairs','rgl','persp')}.
#'
#' @param ngrid vector, number of grid points.
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' Plots to device, or else object of class \code{'zzz'}.
#'
#' @details
#' Where:\cr
#'      X = n x m matrix,\cr
#'      EPDF = empirical probability distribution function (density),
#'      and\cr
#'      ECDF = empirical cumulative distribution function;\cr\cr
#'      Then three possible multivariate quantiles are:\cr
#'      [1] Marginal quantile: from ECDF of raw data for EACH axis
#'      independently (so yields m separate vectors each of length
#'      n).\cr
#'      [2] Joint quantile: from ECDF of raw data across ALL m axes
#'      simultaneously (so yields 1 vector of length n); always
#'      monotonically increasing toward higher axis values.\cr
#'      [3] Conditional quantile: from the ECDF of the EPDF of ALL m
#'      axes simultaneously (so yields 1 vector of length n);
#'      monotonically increasing toward lower \emph{density} values,
#'      but may vary with respect to \emph{axis} values.\cr\cr
#' Conditional quantiles are analogous to data depth methods; current
#'      implementation allows 1-6 dimensions. Joint quantiles are
#'      analogous to Pareto frontiers; current  implementation allows
#'      1-3 dimensions.
#'
#' @examples
#' # iris data
#' x <- iris[,1:3]
#' cond_cdf(x, 'pairs')
#' joint_cdf(x, 'pairs')
#' cond_cdf(x, 'rgl')
#' joint_cdf(x, 'rgl')
#'
#' # dustbunny data
#' set.seed(23)
#' x <- data.frame(q = rnorm(99,0,5)^2,
#'                 r = rnorm(99,0,5)^2,
#'                 s = rnorm(99,10,10)^2)
#' cond_cdf(x, 'pairs')
#' joint_cdf(x, 'pairs')
#' cond_cdf(x, 'rgl')
#' joint_cdf(x, 'rgl')
#'
#' @export
#' @rdname mquantile
###   conditional CDF: percentiles of an EPDF across 1-6 dimensions
`cond_cdf` <- function(x, pltype, ngrid, ...){
     if(!is.data.frame(x)) x <- as.data.frame(x)
     na <- anyNA(x)
     if(na){
          nr <- nrow(x)
          x  <- na.omit(x)
          ix <- attr(x, 'na.action')
          cat(length(ix), 'NA rows were removed!')
     }
     nc <- dim(x)[[2]]
     pp <- c('none','cdf','pairs','rgl','persp')
     if(missing(pltype)) {
          pltype <- 'none'
     } else {
          pltype <- pp[pmatch(pltype, pp)]
     }
     if(pltype=='persp' & nc==2){
          if(missing(ngrid))  ngrid <- 44
          f <- kde(x=x, gridsize=ngrid) # fits across uniform grid
          persp(f$estimate,
                xlab=f$names[1],
                ylab=f$names[2],
                zlab='Prob. density',
                col=ecole::surfcol(f$estimate, ngrid=ngrid, ...), ...)
          return(NULL)
     }
     if(pltype=='persp' & nc>2){
          message('more than 2 columns, using `pairs` not `persp`')
          pltype <- 'pairs'
     }
     f <- kde(x=x, eval.points=x)  # EPDF of orig pts
     z <- f$estimate               # density values of raw data
     if(na){                       # put values in appropriate rows
          v  <- vector('numeric', nr)
          v[ 1:nr %in% ix] <- NA
          v[!1:nr %in% ix] <- z
          z <- v
          m <- matrix(NA, nrow=nr, ncol=nc)
          m[!1:nr %in% ix, ] <- as.matrix(f$eval.points)
          f$eval.points <- as.data.frame(m)
     }
     e <- ecdf(z)   # ECDF  of the *DENSITY* values! (not raw data)
     p <- e(z)      # %iles of the *DENSITY* values! (not raw data)
     p <- (1-p)     # take complement (higher values more 'extreme')
     o <- cbind(f$eval.points, den=z, p=p)
     if(pltype=='none'){
          return(o)
     }
     ccc <- ecole::colvec(p, ...)
     if(pltype=='cdf'){
          plot(z, p, pch=16, col=ccc, bty='l',
               las=1, xlab='Density values',
               ylab='Joint percentile of density values ')
     }
     if(pltype=='rgl' & nc==3){
          rgl::plot3d(o[,1], o[,2], o[,3], pch=16, size=10, col=ccc,
                      xlab=f$names[1], ylab=f$names[2],
                      zlab=f$names[3], box=F)
     } else {
          if(pltype=='rgl' & nc!=3){
               cat('`rgl` needs 3 columns, plotting `pairs` instead')
               pltype <- 'pairs'
          }
     }
     if(pltype=='pairs' & nc>2){
          pairs(f$eval.points, upper.panel=NULL, pch=16, col=ccc, ...)
     }else{
          if(pltype=='pairs' & nc==2){
               plot(f$eval.points, pch=16, col=ccc, ...)
          }
     }
}
#' @export
#' @rdname mquantile
###   joint CDF: over ALL dimensions simultaneously, for 1-3 dimensions
`joint_cdf` <- function(x, pltype, ngrid, ...){
     if(!is.data.frame(x)) x <- as.data.frame(x)
     nc <- dim(x)[[2]]
     if(nc > 3) stop('only works for 1-3 dimensions')
     na <- anyNA(x)
     if(na){
          nr <- nrow(x)
          x  <- na.omit(x)
          ix <- attr(x, 'na.action')
          cat(length(ix), 'NA rows were removed!')
     }
     pp <- c('none','cdf','pairs','rgl','persp')
     if(missing(pltype)) {
          pltype <- 'none'
     } else {
          pltype <- pp[pmatch(pltype, pp)]
     }
     if(pltype=='persp' & nc==2){
          if(missing(ngrid))  ngrid <- 44
          f <- kcde(x=x, gridsize=ngrid)# fits across uniform grid
          persp(f$estimate,
                col=ecole::surfcol(f$estimate, ngrid=ngrid, ...), ...)
          return(NULL)
     }
     if(pltype=='persp' & nc>2){
          message('more than 2 columns, using `pairs` not `persp`')
          pltype <- 'pairs'
     }
     f <- kcde(x=x, eval.points=x) # ECDF of orig pts
     p <- f$estimate  # %iles of raw data, higher values more extreme
     if(na){                       # put values in appropriate rows
          v  <- vector('numeric', nr)
          v[ 1:nr %in% ix] <- NA
          v[!1:nr %in% ix] <- p
          p <- v
          m <- matrix(NA, nrow=nr, ncol=nc)
          m[!1:nr %in% ix, ] <- as.matrix(f$eval.points)
          f$eval.points <- as.data.frame(m)
     }
     o <- cbind(f$eval.points, p=p)
     if(pltype=='none'){
          return(o)
     }
     ccc <- ecole::colvec(p, ...)
     if(pltype=='cdf'){
          auto_rowcol <- function(n = nc) {
               if (n <= 3)
                    c(1, n)
               else if (n <= 6)
                    c(2, (n + 1)%/%2)
               else if (n <= 12)
                    c(3, (n + 2)%/%3)
               else c(ceiling(n/(nr <- ceiling(sqrt(n)))), nr)
          }
          op <- par(mfrow=auto_rowcol())
          for(i in 1:nc){
               plot(x[,i], p, pch=16, col=ccc, bty='l', las=1,
                    xlab=f$names[i], ylab='Marginal percentile')
          }
          par(op)
     }
     if(pltype=='rgl' & nc==3){
          rgl::plot3d(o[,1],o[,2],o[,3], pch=16, size=10, col=ccc,
                      xlab=f$names[1], ylab=f$names[2],
                      zlab=f$names[3], box=F)
     } else {
          if(pltype=='rgl' & nc!=3){
               cat('`rgl` needs 3 columns, plotting `pairs` instead')
               pltype <- 'pairs'
          }
     }
     if(pltype=='pairs' & nc>2){
          pairs(f$eval.points, upper.panel=NULL, pch=16, col=ccc, ...)
     }else{
          if(pltype=='pairs' & nc==2){
               plot(f$eval.points, pch=16, col=ccc, ...)
          }
     }
}