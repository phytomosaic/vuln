#' @title Utility functions for vulnerability analysis
#'
#' @description
#' Various utility functions.
#'
#' @param spe matrix or data.frame, species occurrences or abundances.
#'
#' @param y vector, environmental values matching \code{ybin}.
#'
#' @param ybin data.frame, of class \code{'bingrid'}.
#'
#' @param na.rm logical, should NAs be removed?
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' See below.
#'
#' @details
#' Lightweight functions to calculate vulnerability indices
#'      (\code{litetvi}), specific quantiles (\code{quantileonly}), or
#'      all percentiles (\code{percentileonly}).
#'
#' @examples
#' ### TODO
#'
#' @export
#' @rdname utils
`litetvi` <- function(spe, y, ybin, ...){
     spe <- as.matrix(spe)
     spe[spe>0] <- 1
     spe[spe==0] <- NA
     nr <- nrow(spe)
     nc <- ncol(spe)
     v  <- p <- data.frame(spe * y)
     p[] <- NA
     vest <- data.frame(as.matrix(ybin))
     e <- apply(vest,2,function(x)ecdf(x))
     for(j in 1:nc) p[,j] <- e[[j]](v[,j])        # current %iles
     p95 <- apply(vest, 2, quantile, probs=0.95, ...) # upper limit
     vo <- p                     # percentiles of vulnerab spp
     vo[vo<0.95] <- NA             # omit non-vulnerable
     vo[vo>0]    <- 1              # vulnerable occurrences mx
     t1 <- as.vector(rowSums(vo, ...)/rowSums(spe, ...))*100
     t2 <- rowMeans(p,...)*100   # site avg of spp percentiles
     vc <- sweep(spe, 2, p95, `*`)  # weighted UPPER tols (STIs)
     t3 <- y - rowMeans(vc, ...)# devn of CWM from local MWMT
     cbind(t1=t1, t2=t2, t3=t3)
}
#' @export
#' @rdname utils
`quantileonly` <- function(spe, ybin, na.rm=TRUE, ...){
     tmp <- as.matrix(spe)
     tmp[tmp>0]  <- 1
     tmp[tmp==0] <- NA
     nr  <- nrow(tmp)
     nc  <- ncol(tmp)
     vest<- data.frame(as.matrix(ybin)) # cli vals for niche
     ncb <- ncol(vest)
     if(nc!=ncb)stop('species mismatch between pt, bin data')
     p50 <- apply(vest, 2, median, na.rm=na.rm)
     iqr <- apply(vest, 2, IQR, na.rm=na.rm)
     p95 <- apply(vest, 2, quantile, probs=0.95, na.rm=na.rm)
     p05 <- apply(vest, 2, quantile, probs=0.05, na.rm=na.rm)
     srf <- colSums(tmp,na.rm=T)/nr     # spp rel freq
     out <- cbind(p05=p05, p50=p50, p95=p95, iqr=iqr, srf=srf)
     out <- apply(out, 2, round, digits=2)
     out
}
#' @export
#' @rdname utils
`percentileonly` <- function(spe, y, ybin, ...){
     spe   <- as.matrix(spe)
     spe[spe>0]  <- 1
     spe[spe==0] <- NA
     nr    <- nrow(spe)
     nc    <- ncol(spe)
     v     <- p00 <- data.frame(spe * y)
     p00[] <- NA
     vest  <- data.frame(as.matrix(ybin))
     e     <- apply(vest,2,function(x)ecdf(x))
     for(j in 1:nc) p00[,j] <- e[[j]](v[,j])     # current %iles
     gc()
     p00
}

###################################################################
### unexported utility functions:
###################################################################

`describe` <- function (x, na.rm = TRUE, digits = 2, type = 1, ...) {
     if (!is.numeric(x)) {
          return(NULL)
     }
     m <- mean(x, na.rm = na.rm)
     s <- stats::sd(x, na.rm = na.rm)
     v <- stats::var(x, na.rm = na.rm)
     na <- sum(is.na(x))
     n <- length(x) - na
     se <- s/sqrt(n - 1)
     cv <- ifelse(m != 0, s/m * 100, 0)
     skw <- e1071::skewness(x, na.rm = na.rm, type = type)
     krt <- e1071::kurtosis(x, na.rm = na.rm, type = type)
     out <- data.frame(mean = m, sd = s, var = v, sem = se, cv = cv,
                       n = n, NAs = na, skw = skw, krt = krt)
     out <- round(out, digits = digits)
     mode(out) <- 'numeric'
     return(out)
}