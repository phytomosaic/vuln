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
#' Object of class \code{'vuln'}.
#'
#' @examples
#' \dontrun{
#' # VI from 'truncated' niches based on point FIA data
#' v1 <- vuln(spe, y=id$mwmt)
#' v2 <- vuln(spe, y=id$cmd)
#' v3 <- vuln(spe, y=id$fire)
#' v4 <- vuln(spe, y = cbind(mwmt=id$mwmt,cmd=id$cmd,fire=id$fire))
#' # VI from 'extended' niches based on rasterized CNALH herb. data
#' v1 <- vuln(spe, y=id$mwmt, ybin=m_mwmt)
#' v2 <- vuln(spe, y=id$cmd,  ybin=m_cmd)
#' v3 <- vuln(spe, y=id$fire, ybin=m_fire)
#' v4 <- vuln(spe,
#'            y = cbind(mwmt=id$mwmt,cmd=id$cmd,fire=id$fire),
#'            ybin = list(mwmt=m_mwmt, cmd=m_cmd, fire=m_fire))
#' plot(gv(v1,1), gv(v2,1))
#' }
#' @export
#' @rdname vuln
###  general vulnerability function
`vuln` <- function(spe, y, ybin, yinv, na.rm=TRUE, fullinfo=TRUE,
                   ...){
     y  <- f.as.matrix(y)               # force as matrix
     tmp_abu <- tmp <- f.as.matrix(spe) # force as matrix
     tmp <- (tmp>0)*1                   # binary occurrences matrix
     tmp[tmp==0]         <- NA          # replace 0 with NA
     tmp_abu[tmp_abu==0] <- NA          # replace 0 with NA
     nr  <- nrow(tmp)                   # n sites
     nc  <- ncol(tmp)                   # n species
     ncy <- ncol(y)                     # n explanatory variables

     ### compatibility checks
     if(nr!=nrow(y)) stop('incompatible dims: `spe` and `y`')
     if(!missing(ybin)){
          if (!inherits(ybin, 'bingrid'))
               stop('`ybin` must be of class `bingrid`')
          if(is.vector(ybin[[1]])) ybin <- list(ybin)
          if(ncy!=length(ybin))
               stop('incompatible dims: `ybin`, `y`')
          if(ncy!=1){
               if(length(unique(lapply(ybin,dimnames)))!=1)
                    stop('incompatible dimnames in `ybin` list items')
               if(!identical(dimnames(y)[[2]],names(ybin)))
                    stop('verify colnames identical for `ybin`, `y`')
          }
     }

     ### initialize arrays
     v <- pp <- array(NA, dim=c(nr, nc, ncy))
     dimnames(v) <- dimnames(pp) <- dimnames(tmp)

     ### function to calculate percentiles
     `percentile` <- function(j,a,b){
          a <- a[,j]
          b <- b[,j]
          if( all(is.na(a)) |  all(is.na(b))  ){
               cat('no info for variable',jj,'species',j,'\n')
               rep(NA, times=length(b))
          } else {
               ecdf(a)(b)
          }
     }

     ### for each explanatory variable, calculate ECDF, percentile
     for(jj in 1:ncy){
          v[,,jj] <- tmp * c(y[,jj])    # y vals at occupied sites
          if(missing(ybin)){
               vest <- v[,,jj]          # y vals for niche estimate
               cat(paste0('ECDF is from pointwise data: ',
                          nr, ' sites, ', nc ,' species\n'))
          } else {
               vest <- f.as.matrix(ybin[[jj]]) # y vals for niche est
               nrb <- nrow(vest)
               ncb <- ncol(vest)
               if(nc!=ncb)stop('incompatible dims: `spe` and `ybin`')
               cat(paste0('ECDF is from rasterized data: ',
                          nrb, ' sites, ',ncb,' species\n'))
          }
          cat(' variable', jj, 'of', ncy, '\n')
          pp[,,jj] <- sapply(1:nc, FUN=percentile, a=vest, b=v[,,jj])
     }

     ### invert by taking complement (if higher values more 'extreme')
     if(!missing(yinv)){
          if(!identical(length(yinv), ncy)){
               stop('length(yinv) != NCOL(y)')
          }
          r <- which(yinv)
          pp[,,r] <- (1-pp[,,r]) # complement
     }

     ### calculate site average percentile across ALL y variables
     if(ncy!=1){
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
     if(ncy==1){
          p95 <- apply(vest, 2, quantile, probs=0.95, na.rm=na.rm)
          vc  <- sweep(tmp, 2, p95, `*`) # weighted UPPER tols (STIs)
          cwm <- rowMeans(vc, na.rm=na.rm) # CWM(CTI)=site avg STI
          t3  <- y - cwm           # TVIc, devn of CWM from local MWMT
     }

     ### collect output
     out <- list(v=v,     # y values at all occupied sites
                 vo=vo,    # vulnerable occurrences mx
                 va=va,     # vulnerable abundances mx
                 t1=t1,       # vulnerability indices
                 t2=t2,        # vulnerability indices
                 t3=NULL        # vulnerability indices
     )
     if(exists('t3')) out$t3 <- t3
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
     if(is.vector(x) & is.numeric(x)){
          return(as.matrix(x, rownames.force=T, dimnames=dimnames(x)))
     }
     if(!all(apply(x, 1, is.numeric)))
          stop('data has non-numeric values')
     dn <- dimnames(x)
     dm <- dim(x)
     x <- as.matrix(x, rownames.force=TRUE, dimnames=dimnames(x))
     if(!identical(dim(x),dm))
          stop('unequal column lengths')
     if(!identical(dimnames(x),dn))
          stop('dimnames not retained')
     x
}