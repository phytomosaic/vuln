#' @title Utility functions for plotting
#'
#' @description
#' Various utility functions.
#'
#' @param x vector of climate values at which a species occurs.
#'
#' @param d density object.
#'
#' @param e ECDF object.
#'
#' @param tau specified quantile level.
#'
#' @param ht height of text label in relative coordinates.
#'
#' @param e1,e2,e3 three ECDF objects.
#'
#' @param xval specified x value to return quantile values.
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' See below.
#'
#' @details
#' Lightweight functions to plot.
#'
#' @export
#' @rdname plotutils
`make_d` <- function(x, ...){
     mn   <- min(x)*0.9
     mx   <- max(x)*1.1
     e    <- ecdf(x)         # ECDF
     xx   <- unique(sort(c(seq(mn,mx,.01), knots(e))))
     col  <- ecole::colvec(e(xx)) # color by percentiles
     d    <- density(x, adjust=1.3, n=length(xx)) # EPDF
     d95  <- data.frame(x=d$x, y=d$y) # for polygon > 95th percentile
     d95  <- d95[d95$x > quantile(e,0.95), ]
     d95  <- rbind(d95, c(min(d95$x), 0))
     list(e=e, d=d, d95=d95, col=col)
}
#' @export
#' @rdname plotutils
`plot_d` <- function(d, plotd=FALSE, ...){
     if(plotd) polygon(d$d, col='#00000020')      # full EPDF
     polygon(d$d95, col='#00000050', border=NA)   # EPDF > 95th
     points(d$d, pch=16, col=d$col, cex=0.3)      # EPDF
}
#' @export
#' @rdname plotutils
`ftxt_tau` <- function(e, tau=0.95, ht, ...){
     xval <- quantile(e,tau)  # value at tau
     arrows(x0=xval, y0=0, x1 = xval, y1 = ht, code=1,
            length = 0.12, angle = 20, lwd=2)
     txt <- bquote(
          paste(italic(x[i]) ==.(sprintf('%.2f', round(xval,2))),
                ', ', tau ==.(sprintf('%.2f', tau)),'       '))
     sw   <- strwidth(txt) *0.4
     sh   <- strheight(txt)*0.6
     rect(xval - sw/2 - 0.0, ht - sh/2 - 0.01,
          xval + sw/2 + 0.0, ht + sh/2 + 0.01,
          col = 'white',  border = NA)
     text(xval, ht, txt, cex=0.7)
}
#' @export
#' @rdname plotutils
`ftxt_val` <- function(e1, e2, e3, xval=18.5, ht, ...){
     arrows(x0=xval, y0=0, x1 = xval, y1 = ht, code=1,
            length = 0.12, angle = 20, lwd=2)
     t0 <- bquote(italic(x[i]) ==.(sprintf('%.2f', round(xval,2))))
     t1 <- bquote(tau[1] == .(sprintf('%.2f', round(e1(xval),2))))
     t2 <- bquote(tau[2] == .(sprintf('%.2f', round(e2(xval),2))))
     t3 <- bquote(tau[3] == .(sprintf('%.2f', round(e3(xval),2))))
     sw   <- strwidth(t0) *0.5
     sh   <- strheight(t0)*0.6
     rect(xval - sw/2 - 0.0, ht - sh*5 - 0.01,
          xval + sw/2 + 0.0, ht + sh/2 + 0.01, col='white', border=NA)
     text(xval, ht,           t0, cex=0.7)
     text(xval, ht-(0.07*ht), t1, cex=0.7)
     text(xval, ht-(0.14*ht), t2, cex=0.7)
     text(xval, ht-(0.21*ht), t3, cex=0.7)
}
