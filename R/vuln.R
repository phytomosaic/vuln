#' @title Multidimensional vulnerability index
#'
#' @description
#' Multidimensional vulnerability index, using truncated or extended
#'     niches.
#'
#' @param spe species abundance or occurrence matrix (response
#'     variables).
#'
#' @param y vector, matrix, or list of explanatory variables (climate,
#'     fire, etc), where rows correspond to rows in \code{spe}.
#'
#' @param ybin optional, data.frame or list of data.frames from
#'     \code{bingrid}, with binned explanatory variables from
#'     rasterizing.  Function gives extended niches if supplied,
#'     possibly truncated niches if not supplied.
#'
#' @param yinv vector of logicals, in order of \code{y}, defining
#'     which explanatory variables y should be `inverted` (e.g., low
#'     moisture = high vulnerability). Defaults to
#'     \code{rep(FALSE, NCOL(y))}.
#'
#' @param na.rm logical, indicating whether to strip NA values before
#'     computation. Default is \code{na.rm = TRUE}
#'
#' @param fullinfo logical, default \code{fullinfo=TRUE} propagates
#'     any \code{NA} cases that occur in any \code{y} variable;
#'     when \code{fullinfo=TRUE}, percentiles will be based only on
#'     non-NA \code{y} cases.
#'
#' @param obj object of class \code{'vuln'}.
#'
#' @param pick numeric, indicating which vulnerability metric (1,2,3).
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' List object of class \code{'vuln'}, containing:
#' \describe{
#'   \item{\code{v}}{\code{y} values at all sites in \code{spe}.}
#'   \item{\code{ve}}{\code{y} values at all sites in \code{ybin}.}
#'   \item{\code{vo}}{\strong{vulnerable} occurrences matrix.}
#'   \item{\code{va}}{\strong{vulnerable} abundances matrix.}
#'   \item{\code{t1}}{values for 1st vulnerability index, describing
#'       the proportion of observed species which are vulnerable.}
#'   \item{\code{t2}}{values for 2nd vulnerability index, describing
#'       the community-mean percentile at the observed site.}
#'   \item{\code{t3}}{values for 3rd vulnerability index, describing
#'       the safety margin as the deviation of community-mean upper
#'       \code{y} limits from the observed \code{y} value.}
#' }
#'
#' @examples
#' ## Load data
#' \donttest{
#' data(braun, package='ecole')
#' spe <- braun$spe
#' env <- braun$env
#' ## Vulnerability estimates
#' v <- vuln(spe, y=env$bio1) # w/o bingrid, may be 'truncated' niche
#' ## Mapping
#' d <- data.frame(env, t1 = gv(v,1), t2 = gv(v,2), t3 = gv(v,3))
#' f <- function(d, j='bio1', ...){
#'      cc <- c("#420A6899","#65156E99","#89226A99",
#'              "#AD305D99","#CC424899","#E75D2E99",
#'              "#F7821299","#FCAA0F99","#F6D64599")
#'      plot(d$lat, d$lon, col=cc[cut(d[,j],breaks=9)],pch=16,main=j)
#' }
#' par(mfrow=c(2,2), bty='L', las=1)
#' f(d, 'bio1')
#' f(d, 't1')
#' f(d, 't2')
#' f(d, 't3')
#' }
#'
#' @export
#' @rdname vuln
###  general vulnerability function
`vuln` <- function(spe, y, ybin, yinv, na.rm=TRUE, fullinfo=TRUE,...){
     y   <- f.as.matrix(y)              # force as matrix
     tmp_abu <- tmp <- f.as.matrix(spe) # force as matrix
     tmp <- (tmp>0)*1                   # binary occurrences matrix
     tmp[tmp==0]         <- NA          # replace 0 with NA
     tmp_abu[tmp_abu==0] <- NA          # replace 0 with NA
     nr  <- nrow(tmp)                   # n sites
     nc  <- ncol(tmp)                   # n species
     ncy <- ncol(y)                     # n explanatory variables

     ### compatibility checks
     if (nr!=nrow(y)) stop('incompatible dims: `spe` and `y`')
     if (!missing(ybin)){
          if (!inherits(ybin, 'bingrid'))
               stop('`ybin` must be of class `bingrid`')
          if (is.vector(ybin[[1]])) ybin <- list(ybin)
          if (ncy!=length(ybin))
               stop('incompatible dims: `ybin`, `y`')
          if (ncy!=1){
               if (length(unique(lapply(ybin,dimnames)))!=1)
                    stop('incompatible dimnames in `ybin` list items')
               if (!identical(dimnames(y)[[2]],names(ybin)))
                    stop('verify colnames identical for `ybin`, `y`')
          }
     }

     ### initialize arrays
     v <- pp <- array(NA, dim=c(nr, nc, ncy))
     dimnames(v) <- dimnames(pp) <- dimnames(tmp)

     ### function to calculate percentiles
     `ptile` <- function(j,a,b){
          a <- a[,j]
          b <- b[,j]
          if (all(is.na(a)) | all(is.na(b))){
               message(paste0('no info for variable ', jj,
                              ', species ',j))
               rep(NA, times=length(b))
          } else {
               ecdf(a)(b)
          }
     }

     ### for each explanatory variable, calculate ECDF, percentile
     for(jj in 1:ncy){
          v[,,jj] <- tmp * c(y[,jj])    # y vals at occupied sites
          if (missing(ybin)){
               vest <- v[,,jj]          # y vals for niche estimate
               message(paste0('ECDF from pointwise data: ',
                              nr, ' sites, ', nc ,
                              ' species. Niche may be truncated.'))
          } else {
               vest <- f.as.matrix(ybin[[jj]]) # y vals for niche est
               nrb <- nrow(vest)
               ncb <- ncol(vest)
               if (nc!=ncb) stop('incompatible dims: `spe` and `ybin`')
               message(paste0('ECDF from rasterized data: ',
                              nrb, ' sites, ', ncb,
                              ' species. Niche is extended.'))
          }
          message(paste0(' variable ', jj, ' of ', ncy))
          pp[,,jj] <- sapply(1:nc, FUN=ptile, a=vest, b=v[,,jj])
     }

     ### invert by taking complement (if higher values more 'extreme')
     if (!missing(yinv)){
          if (!identical(length(yinv), ncy)){
               stop('length(yinv) != NCOL(y)')
          }
          r <- which(yinv)
          pp[,,r] <- (1-pp[,,r]) # complement
     }

     ### calculate site average percentile across ALL y variables
     if (ncy!=1){
          p <- rowMeans(pp, na.rm = !fullinfo, dims = 2)
     } else {
          p <- drop(pp)
          v <- drop(v)
     }

     ### occurrences and abundances of species that are locally >95%
     vo <- p                       # initialize
     vo[vo<0.95] <- NA             # omit non-vulnerable
     vo[vo>0]    <- 1              # vulnerable occurrences mx
     va <- vo * tmp_abu            # vulnerable abundances mx

     ### calculate vulnerability indices
     t1  <- as.vector(
          rowSums(vo,na.rm=na.rm)/rowSums(tmp,na.rm=na.rm))*100 # TVIa
     t2  <- rowMeans(p,na.rm=na.rm)*100 # TVIb, site avg of spp %tiles
     # calc vulnerability index C only if y is univariate
     if (ncy==1){
          p95 <- apply(vest, 2, quantile, probs=0.95, na.rm=na.rm)
          vc  <- sweep(tmp, 2, p95, `*`) # weighted UPPER tols (STIs)
          cwm <- rowMeans(vc, na.rm=na.rm) # CWM(CTI)=site avg STI
          t3  <- y - cwm           # TVIc, devn of CWM from local MWMT
     }

     ### collect output
     out <- list(v=v,    # y values at all occupied sites
                 ve=NULL, # extended y values from `ybin`
                 vo=vo,    # vulnerable occurrences mx
                 va=va,     # vulnerable abundances mx
                 t1=t1,      # vulnerability indices
                 t2=t2,       # vulnerability indices
                 t3=NULL       # vulnerability indices
     )
     if (exists('t3'))       out$t3 <- t3
     if (!missing(ybin))  out$ve <- vest
     class(out) <- 'vuln'
     out
}
#' @export
#' @rdname vuln
### extractor function to get vulnerability values
`gv` <- function(obj, pick=1, ...){
     stopifnot(class(obj) == 'vuln')
     return(unlist(obj[c('t1','t2','t3')[pick]], use.names=FALSE))
}

#####   UNEXPORTED FUNCTIONS   ##################################

### function always returns a matrix, even if passed a vector
`f.as.matrix` <- function(x, ...){
     if (is.vector(x) & is.numeric(x)){
          return(as.matrix(x, rownames.force=T, dimnames=dimnames(x)))
     }
     if (!all(apply(x, 1, is.numeric)))
          stop('data has non-numeric values')
     dn <- dimnames(x)
     dm <- dim(x)
     x <- as.matrix(x, rownames.force=TRUE, dimnames=dimnames(x))
     if (!identical(dim(x),dm))
          stop('unequal column lengths')
     if (!identical(dimnames(x),dn))
          stop('dimnames not retained')
     x
}