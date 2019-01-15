#' @title Sensitivity
#'
#' @description
#' Sensitivity of vulnerability outputs to uncertainty in climate
#'    variable inputs.
#'
#' @param xvec vector of climate values.
#'
#' @param ybin data.frame, of class \code{'bingrid'}..
#'
#' @param perc numeric, amount of uniformly distributed noise to add,
#'    expressed as a percentage of the range of \code{xvec}
#'
#' @param nrep integer, number of replicates.
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' data.frame of class \code{'sens'}.
#'
#' @details
#' Evaluates effects of uncertainty in climate inputs on the
#'     vulnerability values.
#'
#' @seealso
#' \link[vuln]{utils} for internal functions.
#'
#' @export
#' @rdname sens
`sens` <- function(xvec, ybin, perc, nrep, ...){
     x0 <- litetvi(spe, xvec, ybin, na.rm=T) # baseline iteration
     `f` <- function(xvec, perc, ybin, ...){
          xx <- litetvi(spe,noise(z=xvec,perc=perc),ybin=ybin,na.rm=T)
          c(mae(x0[,1], xx[,1], stdz=T),
            mae(x0[,2], xx[,2], stdz=T),
            mae(x0[,3], xx[,3], stdz=T))
     }
     Q <- matrix(nrow=nrep, ncol=3) # initiate sensitivity matrix
     for(i in 1:nrep){
          cat('round',i,'of',nrep,'\n')
          Q[i,] <- f(xvec, perc, ybin)  # add 5% uncertainty, get Q
     }
     Q
}
