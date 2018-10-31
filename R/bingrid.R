#' @title Bin to grid
#'
#' @description
#' Bin observations into grid cells by calculating mean of
#'     environmental values per species per grid cell of a raster
#'     (i.e., \link[raster]{rasterize}). Useful for aggregating many
#'     herbarium records into larger areas, given that herbarium
#'     sampling is often geographically unrepresentative.
#'
#' @param x long-format data.frame where each row is one observation,
#'      perhaps from herbarium records or plot sampling; minimally
#'      contains columns named \code{'spp'}, \code{'lat'},
#'      \code{'lon'}; typically also has further columns for site
#'      descriptors such as climate or fire history.
#'
#' @param field character, name of column in \code{'spp'} for which to
#'      aggregate values.
#'
#' @param nr,nc integer > 0, number of rows or columns.
#'
#' @param xmn minimum x coordinate (left border).
#'
#' @param xmx maximum x coordinate (right border).
#'
#' @param ymn minimum y coordinate (bottom border).
#'
#' @param ymx maximum y coordinate (top border).
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' data.frame of class \code{'bingrid'}.
#'
#' @details
#' For each species in column 'spp', this function rasterizes values
#'      by taking the aggregate mean of values per grid cell.
#'
#' @examples
#' # iris data
#' set.seed(1959)
#' x <- iris
#' n <- NROW(x)
#' x <- data.frame(x,
#'                 lon = seq(-80,-120,len=n) + rnorm(n,0,2),
#'                 lat = seq(30,45,len=n) + rnorm(n,0,2))
#' colnames(x)[colnames(x) == 'Species'] <- 'spp'
#' # (mmap(x, 'Sepal.Length', alpha=1, name='Sepal length'))
#' bin_sepal <- bingrid(x, nr=30, nc=50, field='Sepal.Length',
#'                      xmn=-125,xmx=-75,ymn=20,ymx=55)
#' r <- raster(nrows=30, ncols=50,
#'                     xmn=-125, xmx=-75, ymn=20, ymx=55)
#' r[] <- 0
#' s <- stack(r,r,r)
#' for (i in 1:3){ s[[i]][] <- bin_sepal[[i]] }
#' plot(s, col=viridis::inferno(99))
#'
#' @seealso
#' \link[raster]{rasterize} for core function.
#'
#' @export
#' @rdname bingrid
`bingrid` <- function(x, field='', nr, nc, xmn, xmx, ymn, ymx, ...){
     nm <- names(x)
     hascoord <- ('lat' %in% nm) & ('lon' %in% nm) & ('spp' %in% nm)
     r <- raster::raster(nrows=nr, ncols=nc,
                         xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx)
     if (!hascoord) stop('need names `lat`,`lon`,`spp` in data')
     out <- plyr::dlply(
          x,
          'spp',
          .fun=function(x){
               x <- data.frame(x[,!names(x)%in%'spp',drop=T])
               sp::coordinates(x) <- c('lon','lat')
               raster::rasterize(x, r, field=field, fun=mean, na.rm=T)
          },
          .progress='text',
          .drop=F)
     out <- data.frame(t(plyr::laply(out, .fun=raster::getValues)))
     class(out) <- c(class(out), 'bingrid')
     out
}
