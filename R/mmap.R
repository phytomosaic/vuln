#' @title Custom mapping functions
#'
#' @description
#' Map geographic locations as points colored by values.
#'
#' @param x data.frame of observations.
#'
#' @param ptcol character, name of column in \code{x} by which to
#'      color points.
#'
#' @param name character, name for color scale.
#'
#' @param title character, main title.
#'
#' @param alpha numeric vector, values used for relative alpha
#'      transparency.
#'
#' @param size numeric, point size.
#'
#' @param xlim numeric vector, x limits (longitude).
#'
#' @param ylim numeric vector, x limits (latitude).
#'
#' @param tt logical, two-tone color scale?
#'
#' @param bg character, background fill for land area (defaults to
#'      transparent).
#'
#' @param low character, lower end of color scale.
#'
#' @param high character, upper end of color scale.
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' Plots to device, or else object of class \code{'zzz'}.
#'
#' @details
#' Plot continuous (\code{'mmap'})or discrete values (\code{'dmap'}).
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
#' (dmap(x, 'Sepal.Big',  name='Sepal > 6 mm'))
#'
#' @import ggplot2
#' @export
#' @rdname mmap
### plot regional maps (continuous values)
`mmap` <- function(x, ptcol, name=ptcol, title='', alpha=NULL, size=1,
                   xlim=c(-122,-75), ylim=c(25,49.5), tt=F,
                   bg=NA, low='white', high='black', ...){
     eb <- element_blank()
     hascoord <- ('lat' %in% names(x)) & ('lon' %in% names(x))
     if (!hascoord) stop('need `lat` and `lon` in `x`')
     if (!is.numeric(alpha)){
          alf <- ecole::standardize(x[,ptcol])
     } else {
          alf <- alpha
     }
     p <- ggplot(x, aes(x=lon, y=lat)) + labs(x='',y='',title=title)+
          guides(alpha=F) + borders('state', colour='black', fill=bg)+
          coord_map('albers', 37, 49.5, xlim=xlim, ylim=ylim) +
          geom_point(aes_string(colour=ptcol), alpha=alf,
                     shape=16, size=rel(size)) +
          theme_classic() +
          theme(plot.title=element_text(hjust=0.5),
                plot.background=eb,
                panel.background=eb,
                legend.background=eb,
                legend.position=c(0.2,0.15),
                legend.direction='horizontal',
                legend.key.height=unit(0.03,'npc'),
                legend.key.width=unit(0.02,'npc'),
                legend.title=element_text(size=9),
                legend.text=element_text(size=7),
                axis.line=eb,
                axis.ticks=eb,
                axis.text=eb,
                axis.title=eb)
     if(tt){
          p <- p + scale_colour_gradient2(
               name=name, low='white', mid='grey', high='black',
               midpoint=55, na.value='transparent',
               labels=c('0','25','50','75','100'),
               guide=guide_colorbar(title.position='top',
                                    barwidth=5, barheight=.7))
     }else{
          p <- p + scale_colour_gradient(
               name=name, low=low, high=high, na.value='transparent',
               guide=guide_colorbar(title.position='top',
                                    barwidth=5, barheight=.7))
     }
}
#' @export
#' @rdname mmap
### plot regional maps (discrete)
`dmap` <- function(x, ptcol, name=ptcol, title='', alpha=NULL, size=1,
                   xlim=c(-122,-75), ylim=c(25,49.5), bg=NA, ...){
     eb <- element_blank()
     hascoord <- ('lat' %in% names(x)) & ('lon' %in% names(x))
     if (!hascoord) stop('need `lat` and `lon` in `x`')
     if (!is.numeric(alpha)){
          alf <- c(1,0.1)[as.numeric(x[,ptcol])+1]
     } else {
          alf <- alpha
     }
     p <- ggplot(x, aes(x=lon, y=lat)) +
          borders('state', colour='black', fill=bg) +
          coord_map('albers', 37, 49.5, xlim=xlim, ylim=ylim) +
          # geom_point(aes_string(colour=ptcol),
          #            shape=16, size=rel(size)) +
          geom_point(aes_string(colour=ptcol), alpha=alf,
                     shape=16, size=rel(size)) +
          labs(x='',y='', title=title) +
          scale_colour_manual(name=name, values=c('grey','black'),
                              na.value='transparent',
                              guide=guide_legend(
                                   title.position='top',
                                   keywidth=1, keyheight=.7,
                                   override.aes=list(size=2.5))) +
          # scale_alpha_discrete(range=c(.1,1), guide='none') +
          theme_classic() +
          theme(plot.title=element_text(hjust=0.5),
                plot.background=eb,
                panel.background=eb,
                legend.background=eb,
                legend.position=c(0.2,0.15),
                legend.direction='horizontal',
                legend.key.height=unit(0.03,'npc'),
                legend.key.width=unit(0.02,'npc'),
                legend.title=element_text(size=9),
                legend.text=element_text(size=7),
                axis.line=eb,
                axis.ticks=eb,
                axis.text=eb,
                axis.title=eb)
}