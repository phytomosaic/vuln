#' @title Noise
#'
#' @description
#' Add random uncertainty for sensitivity analysis.
#'
#' @param z vector, matrix, or data.frame of observations.
#'
#' @param perc numeric, percent of data range.
#'
#' @param unif logical, use default uniform distribution? Alternative
#'      is to use normal distribution.
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' Matrix of noisified data of same dimensions as \code{'z'}.
#'
#' @details
#' Adds noise, ignoring NA values.  If \code{unif=TRUE}, then
#'      uniformly-distributed noise is bounded within \code{perc}
#'      percent of the original data range.  If \code{unif=FALSE},
#'      then normally-distributed noise has standard deviation equal
#'      to \code{perc} percent of the original data range.  Useful for
#'      perturbing values for sensitivity analysis.
#'
#' @examples
#' # iris data
#' x  <- iris[,1:3]
#' xn <- noise(x, perc=1)
#' pairs(x, lower.panel=NULL)
#' pairs(xn, lower.panel=NULL)
#'
#' @export
#' @rdname noise
### add random uncertainty, bounded by some % of the range
`noise` <- function(z, perc=0.0001, unif=TRUE, ...){
     if(is.vector(z)) wasvec <- TRUE else wasvec <- FALSE
     if(!is.matrix(z)) z <- as.matrix(z)
     nr  <- dim(z)[1]
     rng <- diff(range(z, na.rm=T))*perc/100
     n   <- prod(dim(z))
     if(unif){
          zr <- matrix(runif(n,rng*(-1),rng),nrow=nr)
     }else{
          zr <- matrix(rnorm(n,0,rng),nrow=nr)
     }
     out <- zr + z
     if(wasvec) as.vector(out) else out
}