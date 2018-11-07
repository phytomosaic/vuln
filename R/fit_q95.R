#' @title Fit 95th percentile as nonlinear 2D surface
#'
#' @description
#' Hard-coded to fit nonlinear regression surface.
#'
#' @param data data.frame of observations, with columns named
#'     \code{x1}, \code{x2}, and \code{y}.
#'
#' @param m  fitted model from \code{crs}.
#'
#' @param zlo,zhi upper and lower z-limits for plotting.
#'
#' @param zlab z-axis label.
#'
#' @param txt character, text to add as plot label.
#'
#' @param ... further arguments to \code{persp}.
#'
#' @return
#' \code{fit_q95} returns a fitted \code{crs} object, \code{persp_q95}
#'     plots its perspective plot, and \code{qlab} is a convenience
#'     labeling function.
#'
#' @details
#' Hard-coded, so only for one custom plotting purpose.
#'
#' @import crs
#' @export
#' @rdname fit_q95
`fit_q95` <- function(data, ...){
     m <- crs::crs(
          y~x1+x2, data=data, degree = c(3,3), segments = c(6,6),
          knots = 'quantiles', basis = 'glp', cv='none',
          kernel=FALSE, tau=0.95)
     cat('R2 =', m$r.squared, '\n')
     m
}
`persp_q95` <- function(m, data, zlo, zhi, zlab='\nV1', ...){
     op <- par(cex.axis=0.7, cex.lab=0.9, mar=c(1,0,0,0))
     on.exit(par(op))
     ngrd  <- 23
     xlab  <- '\nElevation (m)'
     ylab  <- '\nLatitude (°)'
     x1seq <- seq(min(data$x1),max(data$x1),length=ngrd)
     x2seq <- seq(min(data$x2),max(data$x2),length=ngrd)
     xgrid <- expand.grid(x1=x1seq, x2=x2seq)
     newd  <- data.frame(x1=xgrid[,1], x2=xgrid[,2])
     z     <- matrix(predict(m, newdata=newd), ngrd, ngrd)
     if(missing(zlo)) zlo <- -30
     if(missing(zhi)) zhi <- max(z)
     persp(x=x1seq, y=x2seq, z=z, xlab=xlab, ylab=ylab, zlab=zlab,
           zlim=c(zlo, zhi), ticktype='detailed', border='#4D4D4D',
           theta=134, phi=35, expand=0.7, shade=0.7, ...)
}
`qlab` <- function(txt='A)', ...){
     mtext(txt,3,-0.5,adj=0.1)
}