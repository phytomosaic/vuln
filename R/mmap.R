#' @title Regional mapping
#'
#' @description
#' Map geographic locations as points colored by values.
#'
#' @param x data.frame of observations.
#'
#' @param field character, name of column in \code{x} by which to
#'      color points.
#'
#' @param name character, name for color scale.
#'
#' @param title character, main title.
#'
#' @param alpha numeric vector, values used for relative alpha
#'      transparency.
#'
#' @param sizept,sizetxt numeric, relative size for points and text.
#'
#' @param xlim numeric vector, x limits (longitude in decimal
#'     degrees).
#'
#' @param ylim numeric vector, y limits (latitude in decimal degrees).
#'
#' @param bg character, background fill for land area (defaults to
#'      transparent).
#'
#' @param colorscale character, one of \code{c('bw','bw2','inferno',
#'     'discretebw')} defining the color scale.
#'
#' @param ... further arguments currently ignored.
#'
#' @return
#' Returns or plots an object inheriting from class \code{'ggplot'}.
#'
#' @details
#' Plot continuous or discrete values mapped across regions of North
#'     America.
#'
#' @examples
#' # iris data
#' set.seed(1959)
#' x <- iris[,1:3]
#' n <- NROW(x)
#' # add geographic structure
#' x <- data.frame(
#'      x,
#'      lon = seq(-80,-120,len=n) + rnorm(n,0,2),
#'      lat = seq(30,45,len=n) + rnorm(n,0,2),
#'      Sepal.Big = as.factor(ifelse(x$Sepal.Length > 6,1,0)))
#' (mmap(x, 'Sepal.Length', alpha=1, name='Sepal length'))
#' (mmap(x, 'Sepal.Length', alpha=1, name='Sepal length',
#'       colorscale='inferno'))
#' (mmap(x, 'Sepal.Big',  name='Sepal > 6 mm'))
#'
#' @import ggplot2
#' @export
#' @rdname mmap
### plot regional maps (continuous values)
`mmap` <- function(x, field, name=field, title='', alpha=NA, sizept=1,
                   sizetxt=1, xlim=c(-122,-75), ylim=c(25,49.5),
                   bg=NA, colorscale, ...){
     hascoord <- ('lat' %in% names(x)) & ('lon' %in% names(x))
     if (!hascoord) stop('need `lat` and `lon` in `x`')
     eb <- element_blank()
     isnum <- is.numeric(x[ ,field])
     v  <- c('bw','bw2','inferno','discretebw')
     # select appropriate colorscale
     if (!isnum){
          message('values not numeric, using discrete scale')
          colorscale <- 'discretebw'
     } else if (isnum & missing(colorscale)){
          colorscale <- 'bw'
     } else {
          colorscale <- v[pmatch(colorscale, v)]
     }
     # select appropriate alpha level
     isalf <- is.numeric(alpha)
     if (!isnum & !isalf){
          alf <- NA
          # alf <- c(0.1,1)[as.numeric(x[,field])+1]
     } else if (isnum & !isalf){
          alf <- ecole::standardize(x[,field])
     } else if (isnum & isalf){
          alf <- alpha
     }
     p <- ggplot(x, aes(x=lon, y=lat)) + labs(x='',y='',title=title) +
          guides(alpha=F) +
          borders('state', colour='black', fill=bg,
                  size=rel(sizetxt)*0.4) +
          coord_map('albers', 37, 49.5, xlim=xlim, ylim=ylim) +
          geom_point(aes_string(colour=field), alpha=alf,
                     shape=16, size=rel(sizept)) +
          theme_classic() +
          theme(plot.margin = unit(c(-5,1,1,1),'mm'),
                plot.title=element_text(hjust=0.5, vjust=-2,
                                        size=rel(sizetxt)*1),
                plot.background=eb,
                panel.background=eb,
                legend.background=eb,
                legend.position=c(0.2,0.15),
                legend.direction='horizontal',
                legend.key.height=unit(0.03,'npc'),
                legend.key.width=unit(0.02,'npc'),
                legend.title=element_text(size=rel(sizetxt)*0.8),
                legend.text=element_text(size=rel(sizetxt)*0.7),
                axis.line=eb,
                axis.ticks=eb,
                axis.text=eb,
                axis.title=eb)
     if (colorscale == 'discretebw'){
          p <- p +  scale_colour_manual(
               name=name, na.value='transparent',
               values=c('#00000003','#000000E6'),
               guide=guide_legend(
                    title.position='top',
                    keywidth=1, keyheight=.7,
                    override.aes=list(
                         size=rel(sizept)*2.5,
                         colour=c('#0000001A','#000000E6'))))
     } else if (colorscale == 'bw'){
          p <- p + scale_colour_gradient(
               name=name, low='white', high='black',
               na.value='transparent',
               guide=guide_colorbar(title.position='top',
                                    barwidth=rel(sizetxt)*4,
                                    barheight=rel(sizetxt)*0.5))
     } else if (colorscale == 'bw2'){
          message('colorscale = `bw2` is for data in [0,100]')
          p <- p + scale_colour_gradient2(
               name=name, low='white', mid='grey', high='black',
               midpoint=55, na.value='transparent',
               labels=c('0','25','50','75','100'),
               breaks = c(0, 25, 50, 75, 100),
               guide=guide_colorbar(title.position='top',
                                    barwidth=rel(sizetxt)*4,
                                    barheight=rel(sizetxt)*0.5))
     } else if (colorscale == 'inferno'){
          p <- p + scale_color_viridis(
               name=name, na.value='transparent',
               guide=guide_colorbar(title.position='top',
                                    barwidth=rel(sizetxt)*4,
                                    barheight=rel(sizetxt)*0.5),
               begin=0.1, end=0.9, option='B')
     }
}
