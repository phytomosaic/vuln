###################################################################
# Community Vulnerability from Species' Niche Constraints
#      Reviewer-only materials for J Biogeogr manuscript
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 06 Aug 2019
##  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)


### NOTE: this script makes use of packages and local files not
###   available in the online Github package


###   begin preamble   ############################################
rm(list=ls())
if(!require(devtools)) install.packages('devtools')
devtools::install_github('phytomosaic/vuln', upgrade=F)
require(vuln)
require(ecole)
require(raster)
require(vegan)
require(MASS)
require(vuln)
require(ggplot2)
require(grid)
require(gridExtra)

### load data
data(lichen_spe, lichen_id, lichen_mex)
id  <- lichen_id        # FIA sites
spe <- lichen_spe       # FIA species
a <- x <- lichen_mex    # CNALH occurrences
rm(lichen_id, lichen_mex)
wd <- 'C:/Users/Rob/Desktop/'
setwd(wd)


###################################################################
### climate DENSITIES in N Am, CNALH, FIA approximately similar ###
### make regular grid of points to query ClimateNA
range(x$lat) ; range(x$lon)
g <- expand.grid(lat=seq(14.5, 83.5, len=500),
                 long=seq(-178, -52.5, len=500))
ecole::prep_climatena(ID1=1:NROW(g), lat=g$lat, long=g$long,
                      file='./to_climatena.csv')
### ... go to ClimateNA to query values...
gg <- read.csv('./from_climatena.csv')
names(gg) <- tolower(names(gg))
gg <- gg[gg$mwmt != (-9999), names(gg) %in% c('mwmt','cmd'), ]
gg <- rbind(gg, cbind(mwmt=x$mwmt,cmd=x$cmd),
            cbind(mwmt=id$mwmt,cmd=id$cmd))
# tiff('./climate_density.tif', 8.5, 3.5,
#      unit='in', compr='lzw+p',res=700)
yl <- 'CMD (mm)'
xl <- 'MWMT (\u00B0C)'
ym <- c(0,1500)
xm <- c(0,max(gg$mwmt)+1)
lv <- c(1,2,4,8)*0.00001
cc <- colvec(1:5, alpha=1)
ecole::set_par(3)
par(mgp=c(2.7, 0.5, 0))
plot(0,type='n',ylim=ym,xlim=xm,ylab=yl,xlab=xl,main='All N. America')
contour(kde2d(gg$mwmt,gg$cmd),levels=lv,drawlabels=F,add=T,col=cc,lwd=2)
plot(0,type='n',ylim=ym,xlim=xm,ylab=yl,xlab=xl,main='Herbarium data')
contour(kde2d(x$mwmt,x$cmd),levels=lv,drawlabels=F,add=T,col=cc,lwd=2)
plot(0,type='n',ylim=ym,xlim=xm,ylab=yl,xlab=xl,main='FIA plot data')
contour(kde2d(id$mwmt,id$cmd, n=20),levels=lv,drawlabels=F,add=T,
        col=cc,lwd=2)
# dev.off()
rm(gg,yl,xl,ym,xm,lv,cc)
###################################################################


###################################################################
###   FIA COMPLETENESS   ##########################################
# how does COMPLETENESS of species capture vary by grain size?

### FIA COMPLETENESS among 6,474 FIA plots
### estimate sampling completeness at FIA sites
id$o <- vegan::specnumber(spe) # observed species richness
id$e <- vegan::estimateR(spe)['S.ACE',] # Chao expected richness
id$e[is.na(id$e)] <- id$o[is.na(id$e)]
id$cscore <- id$o/id$e # completeness = observed / expected
# ### plot
# tiff('./FIA_Cscores.tif', 4, 4.5, unit='in',
#     compr='lzw+p', res=500)
set_par(1)
hist(id$cscore, breaks=44, xlab='C-statistic', col='#00000050',
     main='Completeness of sampling\namong 6,474 FIA plots', las=1)
box(bty='L')
# dev.off()
###################################################################



###################################################################
###   CNALH COMPLETENESS   ########################################
# how does COMPLETENESS of species capture vary by grain size?

### from CNALH keep only coordinates, remove duplicates
x <- x[,c('spp','lat','lon')]
x <- x[!duplicated(x),,]

### create reference raster in mollweide equal-area projection
trg_prj <- CRS('+proj=moll +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84
               +units=m +no_defs')
trg_rez <- c(25000,25000) # 25 x 25 km resolution
trg_e   <- c(-7045000,5045000,0,9025000) # target extent (m)

### expand grain sizes
r25  <- raster(x=extent(trg_e), crs=trg_prj, res=trg_rez)
r25[1:ncell(r25)] <- NA            # 25 km
r50  <- aggregate(r25, fact=2)     # 50 km
r100 <- aggregate(r25, fact=4)     # 100 km
r200 <- aggregate(r25, fact=8)     # 200 km
r400 <- aggregate(r25, fact=16)    # 400 km

### reproject points to mollweide equal-area projection
sp::coordinates(x) <- c('lon', 'lat')
projection(x) <- '+init=epsg:4326'
x <- spTransform(x, projection(r25))

### function estimates C-score of completeness from species richness
`f` <- function(r, x, ...){
        # empty rasters
        if(!all(is.na(r[]))) r[1:ncell(r)] <- NA
        ro <- re <- r
        x$cell <- raster::cellFromXY(re, coordinates(x))  # cell index
        x$abu  <- 1                                       # presence
        # abundances (not pres/abs), as counts per species per cell
        xa <- aggregate(x$abu,by=list(spp=x$spp,cell=x$cell), FUN=sum)
        names(xa) <- c('spp','cell','abu')
        # convert to binary cell x species wide matrix
        xa <- labdsv::matrify(
                data.frame(cell=xa$cell, spp=xa$spp, abu=xa$abu))
        # observed species richness
        o <- vegan::specnumber(xa)
        ro[as.numeric(names(o))] <- o
        # expected species richness
        e <- vegan::estimateR(xa)['S.ACE',] # Chao estimator
        re[as.numeric(names(e))] <- e
        # estimate bin completeness = observed / expected
        rc <- ro/re
        # stack as output
        s <- stack(ro,re,rc)
        names(s) <- c('ro','re','rc')
        s
}
### C-score for each grain size {25,50,100,200,400}
f25  <- f(r25, x)
f50  <- f(r50, x)
f100 <- f(r100, x)
f200 <- f(r200, x)
f400 <- f(r400, x)
cs   <- list(f25,f50,f100,f200,f400)
names(cs) <- c('25 km','50 km','100 km','200 km','400 km')
### plot C-scores
cc <- viridisLite::inferno(999)
# tiff('C_score.tif', 9, 5.25, unit='in', compr='lzw+p', res=500)
ecole::set_par(15)
par(mfrow=c(3,5))
lapply(seq_along(cs), function(x)
        plot(cs[[x]][[3]], col=cc, legend=F, main=names(cs)[[x]]))
lapply(cs, function(x){
        rr <- x[[3]]
        hist(rr, breaks=44, col='#00000060', main='', xlab='C-score',
             xlim=c(0,1), bty='L')
        box(bty='L')
        abline(v=median(rr[], na.rm=TRUE), col=2, lwd=2)
})
lapply(cs, function(x){
        rr <- x[[3]]
        yy <- yFromCell(rr, 1:ncell(rr))
        plot(yy, rr[], col='#00000060', pch=16, cex=0.4,
             xlab='Latitude (m)', ylab='C-score')
        m <- lm(rr[]~yy)
        if(summary(m)$coefficients[1,4] >= 0.5) { abline(m, col=2) }
})
# dev.off()
###################################################################



###################################################################
###  F-values: loss-of-information based on variance   ############
`pvar` <- function(x, na.rm = TRUE, ...){
        if(all(is.na(x))){ return(NA) }
        n <- length(x)
        sum((x - mean(x, na.rm=na.rm))^2, na.rm=na.rm) / n
}
`bing` <- function (x, field = 'mwmt', r, ...){
        x <- by(data    = x[, !colnames(x) %in% 'spp'],
                INDICES = x[, 'spp'],
                FUN = function(j) {
                        sp::coordinates(j) <- c('lon', 'lat')
                        projection(j) <- '+init=epsg:4326'
                        j <- spTransform(j, projection(r))
                        raster::rasterize(j, r, field = field,
                                          fun = pvar, na.rm = TRUE)
                }, simplify = FALSE)
        x <- sapply(x, FUN = raster::getValues)
        class(x) <- c(class(x), 'bingrid')
        x
}
### TIME WARN ! ! ! calculate within-cell variance per species
tot  <- pvar(a$mwmt, na.rm=TRUE) # total variance across all points
v25  <- bing(a, 'mwmt', r25)
v50  <- bing(a, 'mwmt', r50)
v100 <- bing(a, 'mwmt', r100)
v200 <- bing(a, 'mwmt', r200)
v400 <- bing(a, 'mwmt', r400)
### calc F-stat per species
`f` <- function(v,tot) (tot-colMeans(v, na.rm=T))/colMeans(v,na.rm=T)
fval <- rbind(c25  = f(v25, tot),
              c50  = f(v50, tot),
              c100 = f(v100, tot),
              c200 = f(v200, tot),
              c400 = f(v400, tot))
fval[fval == Inf] <- min(fval[fval != Inf]) # correct any Inf values
plot_heatmap(fval, xord='mean', logbase=10)
ff <- rowMeans(fval, na.rm=TRUE) # average F-value per grain-size
### plot
# tiff('./F-values.tif', 6.5, 3, unit='in',
#      compr='lzw+p', res=700)
ecole::set_par(3)
par(mgp=c(2.7, 0.5, 0))
cc <- colvec(1:5, alpha=1)
plot(c(25,50,100,200,400), ff, pch=16,
     xlab='Grain (km)', ylab='F-statistic', col=cc)
plot(c(25,50,100,200,400),
     c(0.4535307,0.5190303,0.6042733,0.6780822,0.7440347), pch=16,
     xlab='Grain (km)', ylab='C-statistic', col=cc)
abline(h=0.50, lty=2)
plot(ff, c(0.4535307,0.5190303,0.6042733,0.6780822,0.7440347),
     pch=16, xlab='F-statistic', ylab='C-statistic', col=cc)
abline(h=0.50, lty=2)
# dev.off()
###################################################################



###################################################################
###  Compare niche truncation when using all WORLD data   #########
rm(list=ls())
p <- '~/_prj/4_vuln/ms/ms_jbiogeogr_2019/03_minor/'
x <- read.csv(paste0(p, 'occurrences_1.csv'), stringsAsFactors=F)
y <- read.csv(paste0(p, 'occurrences_2.csv'), stringsAsFactors=F)
if (identical(names(x), names(y))) {
        x <- rbind(x,y) ; rm(y)
}
x <- x[,c('lat','lon')]
x <- x[!(is.na(x$lat) | is.na(x$lon)),,]
x <- x[!duplicated(x[,c('lat','lon')]),,]
x <- x[!(x$lon > 180 | x$lon < (-180)),,]
x <- x[!(x$lat > 90 | x$lat < (-90)),,]
x <- x[order(x$lat, x$lon),]
dimnames(x)[[1]] <- 1:dim(x)[[1]]
w <- getData('worldclim',var='bio',res=10) # bioclim world data
w <- w[[c(5,12)]]
names(w) <- c('mwmt','map')
pts <- SpatialPoints(cbind(x=x$lon,y=x$lat), proj4string=w@crs)
v <- as.data.frame(extract(w, pts))
v$mwmt <- v$mwmt/10
v <- cbind(x, v) ; rm(x)
v <- v[!is.na(v$mwmt),,]
head(v)
dim(v)

### N. America lichen data
data(lichen_mex)
head(lichen_mex)
x <- lichen_mex ; rm(lichen_mex)
x <- x[,c('lat','lon')]
x <- x[!(is.na(x$lat) | is.na(x$lon)),,]
x <- x[!duplicated(x[,c('lat','lon')]),,]
x <- x[!(x$lon > 180 | x$lon < (-180)),,]
x <- x[!(x$lat > 90 | x$lat < (-90)),,]
pts <- SpatialPoints(cbind(x=x$lon,y=x$lat), proj4string=w@crs)
y <- as.data.frame(extract(w,pts))
y$mwmt <- y$mwmt/10
y <- cbind(x, y) ; rm(x)
y <- y[!is.na(y$mwmt),,]
head(y)
dim(y)

### USA-only lichen data
data(lichen_id)
x <- lichen_id ; rm(lichen_id)
x <- x[,c('lat','lon')]
x <- x[!(is.na(x$lat) | is.na(x$lon)),,]
x <- x[!duplicated(x[,c('lat','lon')]),,]
x <- x[!(x$lon > 180 | x$lon < (-180)),,]
x <- x[!(x$lat > 90 | x$lat < (-90)),,]
pts <- SpatialPoints(cbind(x=x$lon,y=x$lat), proj4string=w@crs)
z <- as.data.frame(extract(w,pts))
z$mwmt <- z$mwmt/10
z <- cbind(x, z) ; rm(x)
z <- z[!is.na(z$mwmt),,]
head(z)
dim(z)

### remove BIOCLIM raster files
unlink('./wc10', recursive=TRUE)

### plot in environmental space
`f` <- function(x, ptcol=NULL, ...){
        if(is.null(ptcol)) ptcol <- '#00000005'
        plot(x[,'lon'], x[,'lat'], pch=16, cex=0.5, col=ptcol,
             ylab='Latitude (\u00B0)', xlab='Longitude (\u00B0)', ...)
}
`g` <- function(x, ptcol=NULL, ...){
        if(is.null(ptcol)) ptcol <- '#00000005'
        plot(x[,'mwmt'], x[,'map'], pch=16, cex=0.5, col=ptcol,
             ylim=ym, xlim=xm, ylab=yl, xlab=xl, ...)
}
`h` <- function(x, lv = c(1,2,4,8)*0.000007, ...){
        cc <- colvec(1:5, alpha=1)
        plot(0, type='n', ylim=ym, xlim=xm, ylab=yl, xlab=xl, ...)
        contour(kde2d(x[,'mwmt'], x[,'map'], n=44), levels=lv,
                drawlabels=F, add=T, col=cc, lwd=2)
}
# tiff('./env_overlap.tif', 6.5, 6.5,
#      unit='in', compr='lzw+p', res=700)
yl <- 'MAP (mm)'
xl <- 'MWMT (\u00B0C)'
ym <- c(0,4000)
xm <- c(0,max(v$mwmt)+1)
ecole::set_par(9)
par(mgp=c(2.7, 0.5, 0))
f(v, main='World') ; f(y, main='N. America')
f(z, ylim=c(25,65), main='USA')
g(v) ; g(y) ; g(z, ptcol='#00000025')
h(v) ; h(y) ; h(z, lv=c(1,2,4,8)*0.000009)
# dev.off()

### niche overlap at varying geographic extents
round(overlap(v$mwmt, y$mwmt), 3)
round(overlap(y$mwmt, z$mwmt), 3)
round(overlap(z$mwmt, v$mwmt), 3)
round(overlap(v$map, y$map), 3)
round(overlap(y$map, z$map), 3)
round(overlap(z$map, v$map), 3)
###################################################################



###################################################################
###   Twitter image   #############################################
data(lichen_spe, lichen_id, lichen_mex)
spe <- lichen_spe   # FIA species
id  <- lichen_id    # FIA sites
mex <- lichen_mex   # CNALH occurrences
rm(lichen_spe, lichen_id, lichen_mex)
m_mwmt <- bingrid(mex, field='mwmt', nr=155, nc=200,
                  xmn=-177.9, xmx=-52.63, ymn=15.37, ymx=82.4)
m_cmd  <- bingrid(mex, field='cmd', nr=155, nc=200,
                  xmn=-177.9, xmx=-52.63, ymn=15.37, ymx=82.4)
rm(mex) ; gc()
v1  <- vuln(spe, y=id$mwmt, ybin=m_mwmt)
v2  <- vuln(spe, y=id$cmd,  ybin=m_cmd)
z   <- cbind(id, Vulnerability=v1$t1)
z$Vulnerability <- ecole::standardize(log1p(z$Vulnerability))
eb  <- element_blank()
theme_set(theme_classic() + theme(
        plot.title=element_text(hjust=0.5,size=rel(0.5),face='bold'),
        axis.text=element_text(
                colour='black',size=rel(0.5),face='bold'),
        axis.title=eb, axis.ticks=eb))
p1  <- mmap(z, 'Vulnerability',
            alpha=ifelse(z$Vulnerability == 0, 0.4, 1),
            title = 'Community vulnerability', colorscale='inferno') +
        theme(legend.position=c(0.7,0.05),
              legend.text=element_text(size=rel(0.5)),
              legend.title=eb)
b    <- as.numeric(v2$v[,'Bryfus'])
dens <- density(b, na.rm=T)
df   <- data.frame(x=dens$x, y=dens$y)
quantiles <- quantile(b, prob=c(0.5, 0.75, 0.90, 0.95, 0.99), na.rm=T)
df$quant <- factor(findInterval(df$x,quantiles))
p2   <- ggplot(df, aes(x,y)) +
        geom_ribbon(aes(ymin=0, ymax=y, fill=quant, colour=quant)) +
        scale_x_continuous(breaks=quantiles, expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_viridis_d(option='B',guide='none',begin=.1,end=.9)+
        scale_colour_viridis_d(option='B',guide='none',begin=.1,end=.9)+
        geom_line() + theme_classic() +
        theme(axis.ticks=eb, axis.text=eb,
              axis.title=element_text(size=rel(0.5), colour='black'),
              plot.title=element_text(size=rel(0.5), face='bold'),
              plot.background=eb, panel.background=eb) +
        ylab('Density') +
        xlab(expression(Temperature~symbol('\256'))) +
        annotate('text', x = 750, y = 0.0015,
                 label='Single-\nspecies\nniche', size=rel(2.5))
# tiff('./twitter_image.tif', wid=4, hei=3,
#      unit='in', compr='lzw+p', res=700)
print(p1)
print(p2, vp=viewport(width=0.35,height=0.35,x=0,y=0,just=c(0,0)))
# dev.off()
####################################################################
####   END   ####
####################################################################