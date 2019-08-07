#' @title Utility functions for plotting conceptual figures
#'
#' @description
#' Utility plotting functions.
#'
#' @param x vector of climate values at which a species occurs.
#'
#' @param d density object.
#'
#' @param plotd logical, plot the density as polygon?
#'
#' @param e ECDF object.
#'
#' @param tau specified quantile level.
#'
#' @param ht height of text label in relative coordinates.
#'
#' @param i character, text in subscript.
#'
#' @param e1,e2,e3 three ECDF objects.
#'
#' @param xval specified x value to return quantile values.
#'
#' @param panels number of figure panels in `set_par`
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' See below.
#'
#' @details
#' Lightweight functions to plot conceptual figures.
#'
#' @export
#' @rdname utils_plot
`make_d` <- function(x, ...){
     mn   <- min(x, ...) * 0.9
     mx   <- max(x, ...) * 1.1
     e    <- ecdf(x)
     xx   <- unique(sort(c(seq(mn,mx,.01), knots(e))))
     col  <- colvec(e(xx))            # color by percentiles
     d    <- density(x, adjust=1.3, n=length(xx)) # EPDF
     d95  <- data.frame(x=d$x, y=d$y) # for polygon > 95th percentile
     d95  <- d95[d95$x > quantile(e, 0.95), ]
     d95  <- rbind(d95, c(min(d95$x), 0))
     list(e=e, d=d, d95=d95, col=col)
}
#' @export
#' @rdname utils_plot
`plot_d` <- function(d, plotd=FALSE, ...){
     if(plotd) polygon(d$d, col='#00000020')      # full EPDF
     polygon(d$d95, col='#00000050', border=NA)   # EPDF > 95th
     points(d$d, pch=16, col=d$col, cex=0.3)      # EPDF
}
#' @export
#' @rdname utils_plot
`txt_tau` <- function (e, tau = 0.95, ht, i, ...) {
     xval <- quantile(e, tau)
     arrows(x0 = xval, y0 = 0, x1 = xval, y1 = ht, code = 1,
            length = 0.12, angle = 20, lwd = 2)
     txt <- bquote(italic(Q[.(i)])(0.95) ==
                        .(sprintf('%.1f', round(xval, 2))))
     sw <- strwidth(txt) * 0.6
     sh <- strheight(txt) * 0.6
     rect(xval - sw/2 - 0, ht - sh/2 - 0.01, xval + sw/2 + 0,
          ht + sh/2 + 0.01, col = 'white', border = NA)
     text(xval, ht, txt, cex = 0.7)
}
#' @export
#' @rdname utils_plot
`txt_val` <- function(e1, e2, e3, xval=18.5, ht, ...){
     arrows(x0=xval, y0=0, x1 = xval, y1 = ht, code=1,
            length = 0.12, angle = 20, lwd=2)
     t0 <- bquote(italic(x[i]) ==.(sprintf('%.1f', round(xval,2))))
     t1 <- bquote(tau[a] == .(sprintf('%.2f', round(e1(xval),2))))
     t2 <- bquote(tau[b] == .(sprintf('%.2f', round(e2(xval),2))))
     t3 <- bquote(tau[c] == .(sprintf('%.2f', round(e3(xval),2))))
     sw <- strwidth(t0) * 0.5
     sh <- strheight(t0) * 0.6
     rect(xval - sw/2 - 0.0, ht - sh*5 - 0.01,
          xval + sw/2 + 0.0, ht + sh/2 + 0.01, col='white', border=NA)
     text(xval, ht,           t0, cex=0.7)
     text(xval, ht - (0.07 * ht), t1, cex=0.7)
     text(xval, ht - (0.14 * ht), t2, cex=0.7)
     text(xval, ht - (0.21 * ht), t3, cex=0.7)
}
#' @export
#' @rdname utils_plot
`set_par` <- function (panels = 1, ...) {
     pty <- 's'
     mgp <- c(2.5, 0.7, 0)
     mar <- c(4, 4, 1.2, 1)
     oma <- c(0, 0, 0, 0)
     auto_rowcol <- function(n = panels) {
          if (n <= 3)
               c(1, n)
          else if (n <= 6)
               c(2, (n + 1)%/%2)
          else if (n <= 12)
               c(3, (n + 2)%/%3)
          else c(ceiling(n/(nr <- ceiling(sqrt(n)))), nr)
     }
     mfrow <- auto_rowcol()
     if (panels > 4)
          panels <- 4
     switch(as.character(panels),
            `4` = par(mfrow = mfrow, mgp = mgp,
                      mar = mar, pty = pty, oma = oma, bty = 'L',
                      las = 1, cex.lab = 1.2, ...),
            `3` = par(mfrow = mfrow, mgp = mgp,
                      mar = mar, pty = pty, oma = oma, bty = 'L',
                      las = 1, cex.lab = 1.4, cex.axis = 1.2, ...),
            `2` = par(mfrow = mfrow,
                      mgp = mgp, mar = mar, pty = pty, oma = oma,
                      bty = 'L', las = 1, cex.axis = 0.85, ...),
            `1` = par(mfrow = mfrow,
                      mgp = mgp, mar = mar, pty = pty, oma = oma,
                      bty = 'L', las = 1, cex.axis = 0.85, ...))
}
#' @export
#' @rdname utils_plot
`add_text` <- function (x, y, labels, bold = FALSE, ...) {
     if (bold) font <- 2 else font <- 1
     text(x = graphics::grconvertX(x, from = 'npc', to = 'user'),
          y = graphics::grconvertY(y, from = 'npc', to = 'user'),
          labels = labels, adj = 0, font = font, ...)
}
#' @export
#' @rdname utils_plot
`colvec` <- function (x, n = 99, alpha = 0.6, begin = 0.2, end = 0.9,
                      dir = 1, pal, ...) {
     if (is.factor(x)) {
          n <- nlevels(x)
     }
     if (missing(pal)) {
          pal <- viridis::inferno(n = n, alpha = alpha, begin = begin,
                                  end = end, direction = dir)
     }
     pal[cut(as.numeric(x), breaks=length(pal), include.lowest=T)]
}
#' @export
#' @rdname utils_plot
`surfcol` <- function (x, ngrid, alpha = 0.6, begin = 0.2, end = 0.9,
                       dir = 1, pal, ...) {
     if (missing(pal)) {
          pal <- viridis::inferno(n = ngrid, alpha = alpha,
                                  begin = begin, end = end,
                                  direction = dir)
     }
     xavg <- (x[-1, -1] +
                   x[-1, -(ncol(x) - 1)] +
                   x[-(nrow(x) - 1), -1] +
                   x[-(nrow(x) - 1), -(ncol(x) - 1)]) / 4
     pal[cut(xavg, breaks = length(pal), include.lowest = T)]
}
#' @export
#' @rdname utils_plot
`plot_loess` <- function (x, y, col = '#00000040', lcol = '#FF0000BF',
                          cex=0.8, pch=16, las=1, bty='L', ...) {
     ow   <- getOption('warn')
     options(warn = -1)
     plot(x = x, y = y, col = col, pch = pch, las = las, bty = bty,
          cex = cex, ...)
     f    <- stats::loess(y ~ x, ...)
     fhat <- predict(f)
     o    <- order(x)
     lines(x[o], fhat[o], col = lcol, lwd = 2)
     options(warn = ow)
}
