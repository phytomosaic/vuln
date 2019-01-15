#####################################################################
# Community Vulnerability from Species' Niche Constraints
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 15 Jan 2019
##  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)


###   begin preamble   ##############################################
if(!require(devtools)) install.packages('devtools')
devtools::install_github('phytomosaic/vuln', upgrade=F)
require(vuln)
require(ggplot2)
require(grid)
require(gridExtra)
###   end preamble   ################################################


###   begin load data   #############################################
data(lichen_spe, lichen_id, lichen_mex)
spe <- lichen_spe   # FIA species
id  <- lichen_id    # FIA sites
mex <- lichen_mex   # CNALH occurrences
rm(lichen_spe, lichen_id, lichen_mex) # clean up
all(identical(row.names(spe), row.names(id)))
prod(dim(spe))-sum(spe==0)              # 75740 unique occurrences
dim(spe)                                # of 443 spp at 6474 sites
dim(mex)[1]                             # 172127 unique occurrences
length(unique(mex$spp))                 # 443 species
dim(unique(cbind(mex$lat,mex$lon)))[1]  # 46343 unique sites
lab1 <- expression(italic(v[1]))
lab2 <- expression(italic(v[2]))
lab3 <- expression(italic(v[3]))
###   end load data  ################################################


###   begin rasterize    ############################################
# rasterize CNALH per 3 most biologically important vars MWMT,CMD,HI
#  (155 rows x 200 columns = approx 50 x 50 km)
m_mwmt <- bingrid(mex, field='mwmt', nr=155, nc=200,
                  xmn=-177.9, xmx=-52.63, ymn=15.37, ymx=82.4)
m_cmd  <- bingrid(mex, field='cmd', nr=155, nc=200,
                  xmn=-177.9, xmx=-52.63, ymn=15.37, ymx=82.4)
###   end rasterize    ##############################################


###   begin calculate vulnerability   ###############################
# VI from 'truncated' niches based on point FIA data
vt1 <- vuln(spe, y=id$mwmt)
vt2 <- vuln(spe, y=id$cmd)
# VI from 'extended' niches based on rasterized CNALH herb. data
v1 <- vuln(spe, y=id$mwmt, ybin=m_mwmt)
v2 <- vuln(spe, y=id$cmd,  ybin=m_cmd)
###   end calculate vulnerability   #################################


###   begin Appendix materials   ####################################
# names of all climate variables under consideration:
nm <- c('mat','mwmt','mcmt','td','map','msp','ahm','shm','ffp','pas',
        'emt','ext','eref','cmd','rh','rh_sm','hi')
lab <- c( 'Latitude (°)',
          'Longitude (°)',
          'Mean annual temperature (°C)',
          'Mean warmest-month temperature (°C)',
          'Mean coldest-month temperature (°C)',
          'Continentality (°C)',
          'Mean annual precipitation (mm)',
          'Mean summer precipitation (mm)',
          'Annual heat moisture index (unitless)',
          'Summer heat moisture index (unitless)',
          'Frost-free period (d)',
          'Precipitation as snow (mm)',
          '30-y extreme minimum temperature (°C)',
          '30-y extreme maximum temperature (°C)',
          'Reference evaporation (mm)',
          'Climatic moisture deficit (mm)',
          'Mean annual relative humidity (%)',
          'Mean warm-season relative humidity (%)',
          'Summer heat index (°C)')
lab <- enc2utf8(lab)

### Appendix S1: summary statistics of raw climate data
# Table S1. Summaries of seventeen climate variables used to construct
# 'extended' niches for North American macrolichen species. Also
# included are latitude and longitude. Climate data (from ClimateNA
# database) was extracted for each location in the herbarium data
# (from the Consortium of North American Lichen Herbaria). Locations
# included the entirety of Canada, United States and Mexico.
cli <- mex[,c('lat','lon',nm)] # CNALH data
tab <- data.frame(
     lab,
     t(sapply(cli, function(x)
          round(c(min=min(x,na.rm=T), max=max(x,na.rm=T)),1))),
     t(sapply(cli, describe, digits = 1, na.rm=T))
)
tab <- tab[,!colnames(tab) %in% c('n','NAs','sem','var')]
tab <- cbind(abbrev=dimnames(tab)[[1]], tab)
names(tab) <- c('Abbreviation','Variable','Minimum','Maximum','Mean',
                'Standard_deviation','Coefficient_of_variation',
                'Skewness','Kurtosis')
tab
write.csv(tab, './fig/S1_climatesummary.csv', row.names = FALSE)
rm(cli, lab, tab)

### Appendix S2: table of climate niche quantiles for ALL (!) 17 vars
###   Climate abbreviations follow table S1.
tab <- lapply(
     nm,
     function(x){ # TIME WARN ! about 10 s per iter
          cat('now doing `', x, '` from', length(nm), 'variables\n')
          quantileonly(
               spe,
               ybin=bingrid(mex, field=x, nr=155, nc=200,
                            xmn=-177.9,xmx=-52.63,ymn=15.37,ymx=82.4))
     })
tab <- Reduce(cbind, tab)
tab <- tab[,!1:NCOL(tab)%in%(grep('srf',colnames(tab))[-1])]
tab <- tab[,c(5,(1:NCOL(tab))[-5])] # move relfreq to 1st column
colnames(tab)[1]  <- 'Relative_frequency'
colnames(tab)[-1] <- paste(colnames(tab)[-1],rep(nm, each=4),sep='_')
tab <- as.data.frame(tab)
data(crosswalk)
cw <- crosswalk
tab$Abbreviation  <- dimnames(tab)[[1]]
tab$Species <- cw[match(tab$Abbreviation,cw$abbrev),'species']
tab <- tab[,c(NCOL(tab)-1, NCOL(tab), 1:(NCOL(tab)-2))]
head(tab, 21)
write.csv(tab, './fig/S2_nichesummary.csv', row.names = FALSE)
rm(nm, cw, tab, mex)

### Appendix S3: niche overlap btwn truncated and extended
tab <- data.frame(
     MWMT_overlap = sapply(names(spe),
                           function(j) overlap(v1$v[,j], v1$ve[,j])),
     CMD_overlap  = sapply(names(spe),
                           function(j) overlap(v2$v[,j], v2$ve[,j])))
data(crosswalk)
cw <- crosswalk
tab$abbrev  <- dimnames(tab)[[1]]
tab$species <- cw[match(tab$abbrev,cw$abbrev),'species']
tab         <- tab[,c(3,4,1,2)] # reorder
tab$relfreq <- round(colSums((spe>0)*1, na.rm = T)/NROW(spe),3)
sum(tab$MWMT_overlap < 0.75)/443 # 260 of 443 extended vs trunc niches
sum(tab$CMD_overlap < 0.75)/443 # 262 of 443 extended vs trunc niches
names(tab) <- c('Abbreviation','Species','MWMT_overlap','CMD_overlap',
                'Relative_frequency')
write.csv(tab, './fig/S3_nicheoverlap.csv', row.names = FALSE)
rm(cw, tab, crosswalk)

### Appendix S4: plot IQR vs tau=0.95
d <- data.frame(quantileonly(spe, ybin=m_mwmt, na.rm=T),
                quantileonly(spe, ybin=m_cmd, na.rm=T))
# tiff('./fig/S4_nichebreadth.tif',wid=6.5,hei=6,units='in',
#      compr='lzw+p', res=700)
pdf('./fig/S4_nichebreadth.pdf', wid=6.5, hei=6)
set_par(6)
par(mfrow=c(2,3), oma=c(0,1,3,0), mar=c(4,3,0,0), mgp=c(2.0,0.6,0))
`f` <- function(...){ vuln:::plot_loess(..., ylab='Niche breadth') }
f(d$p05,d$iqr,xlab=expression(italic(Q)*'(0.05) of MWMT (°C)'))
mtext(expression(bolditalic(tau)=='0.05' ), 3, 1)
mtext('MWMT', 2, 4, las=3, font=2)
f(d$p50, d$iqr, xlab=expression(italic(Q)*'(0.50) of MWMT (°C)'))
mtext(expression(bolditalic(tau)=='0.50' ), 3, 1)
f(d$p95, d$iqr, xlab=expression(italic(Q)*'(0.95) of MWMT (°C)'))
mtext(expression(bolditalic(tau)=='0.95' ), 3, 1)
f(d$p05.1,d$iqr.1,xlab=expression(italic(Q)*'(0.05) of CMD (mm)'))
mtext('CMD', 2, 4, las=3, font=2)
f(d$p50.1,d$iqr.1,xlab=expression(italic(Q)*'(0.50) of CMD (mm)'))
f(d$p95.1,d$iqr.1,xlab=expression(italic(Q)*'(0.95) of CMD (mm)'))
dev.off()
rm(f, d, y)
###   end Appendices material  ######################################


###   begin conceptual figure   #####################################
# make data for environmental values at which 3 species occur
set.seed(423)
x1 <- ((rbeta(400,3,5)*5)+15.11)[sample(1:400,size=100,rep=F)]
x2 <- ((rbeta(400,3,3)*6)+14.71)[sample(1:400,size=100,rep=F)]
x3 <- ((rbeta(400,1.5,1.1)*7)+14.01)[sample(1:400,size=100,rep=F)]
# calculate EPDF and ECDF for each species
d1 <- make_d(x1)
d2 <- make_d(x2)
d3 <- make_d(x3)
### manually calculate vulnerability:
# number of species at site `i`
nspecies <- 3
# upper limits, at each species 95th percentile
(q1 <- quantile(x1, 0.95))
(q2 <- quantile(x2, 0.95))
(q3 <- quantile(x3, 0.95))
# local climate value at site `i`
(xi <- 18.55)
# climate percentile of each species at `xi`
(tau1 <- ecdf(x1)(xi))
(tau2 <- ecdf(x2)(xi))
(tau3 <- ecdf(x3)(xi))
# three vulnerability indices
V1 <- sum(tau1 > 0.95, tau2 > 0.95, tau3 > 0.95) / nspecies
V2 <- mean(c(tau1,tau2,tau3))
V3 <- xi - mean(c(q1,q2,q3))
# summary of vulnerability at site `i`
round(data.frame(V1,V2,V3), 3)
### Figure 1 - conceptual figure
# tiff('./fig/fig_01_concept.tif', wid=3, hei=6, units='in',
#      compr='lzw+p', res=700)
pdf('./fig/fig_01_concept.pdf', wid=3, hei=6)
par(mfrow=c(2,1), oma=c(0,0,1,0), mar=c(2.8,2.8,0,0), mgp=c(1.7,0.5,0),
    bty='L', las=1, cex.axis = 0.85, pty='s')
# upper thermal tolerances differ per species
plot(d1$d,ylim=c(0,max(d1$d$y)*1.75),xlim=c(min(x1)*0.95,max(x1)*1.2),
     main='', col=NA, ylab='Probability density',
     xlab=bquote(italic(X)))
plot_d(d1) ; plot_d(d2) ; plot_d(d3)
txt_tau(d1$e, 0.95, ht=0.67, i='a')
txt_tau(d2$e, 0.95, ht=0.55, i='b')
txt_tau(d3$e, 0.95, ht=0.43, i='c')
mtext(text = 'a)', side = 3, adj = -0.25)
# at a given x-value, calc vulnerability
plot(d1$d,ylim=c(0,max(d1$d$y)*1.75),xlim=c(min(x1)*0.95,max(x1)*1.2),
     main='', col=NA, ylab='Probability density',
     xlab=bquote(italic(X)))
plot_d(d1) ; plot_d(d2) ; plot_d(d3)
txt_val(d1$e, d2$e, d3$e, xval=18.55, ht=0.5)
add_text(.05,.95,bquote(italic(v[1]) == .(sprintf('%.2f',V1))),cex=.9)
add_text(.05,.87,bquote(italic(v[2]) == .(sprintf('%.2f',V2))),cex=.9)
add_text(.05,.79,bquote(italic(v[3]) == -.(sprintf('%.2f',V3*-1))),
         cex=.9)
mtext(text = 'b)', side = 3, adj = -0.25)
dev.off()
rm(x1,x2,x3,d1,d2,d3,nspecies,xi,tau1,tau2,tau3,q1,q2,q3,V1,V2,V3)
###   end conceptual figure   #######################################


###   begin truncated vs extended niches   ##########################
### Figure 2 - truncated vs extended niches
`f` <- function(x, y, pick=1, ...){ # MAE for trunc v extd
     x <- gv(x, pick=pick)
     y <- gv(y, pick=pick)
     l <- range(c(x,y), na.rm=TRUE)
     plot(y, x, xlab='Truncated', ylab='Extended',
          xlim=l, ylim=l, pch=16, col='#00000010', cex=0.6)
     abline(0,1)
     add_text(0.8, 0.2,
              paste0('MAE =\n',round(mae(y,x,stdz=T),3)*100,'%'),
              cex=0.8)
}
# tiff('./fig/fig_02_truncation.tif', wid=6.5, hei=4.25, units='in',
#      compr='lzw+p', res=700)
pdf('./fig/fig_02_truncation.pdf', wid=6.5, hei=4.25)
par(mfcol=c(2,3), oma=c(0,1,3,0), mar=c(4,3,0,0), mgp=c(2.0,0.6,0),
    pty='s', bty='L', las=1, cex.lab=1.2)
f(v1,vt1,1) ; mtext(expression(bolditalic(V[1])), 3, 1)
mtext('MWMT', 2, 4, las=3, font=2)
f(v2,vt2,1) ; mtext('CMD', 2, 4, las=3, font=2)
f(v1,vt1,2) ; mtext(expression(bolditalic(V[2])), 3, 1)
f(v2,vt2,2)
f(v1,vt1,3) ; mtext(expression(bolditalic(V[3])), 3, 1)
f(v2,vt2,3)
dev.off()
rm(f)
# Conclusion: extended niches give generally lower V for all 3 indices
#     and all 3 climate variables; MAE is acceptable.
###   end truncated vs extended niches   ############################


###   begin plot nationwide vulnerability   #########################
### Figure 3 - mapping nationwide vulnerability
tmp <- cbind(id,   # V based on *extended* niches hereafter
             va = v1$va,
             ### MWMT
             t11 = gv(v1,1),
             t12 = gv(v1,2),
             t13 = gv(v1,3),
             t1c = as.factor(ifelse(gv(v1,3)>0,'pos','neg')),
             ### CMD
             t21 = gv(v2,1),
             t22 = gv(v2,2),
             t23 = gv(v2,3),
             t2c = as.factor(ifelse(gv(v2,3)>0,'pos','neg'))
)
names(tmp) <- gsub('va.', '', names(tmp))
# tiff('./fig/fig_03_mapvuln.tif',
#      wid=6.5, hei=3.25, units='in', compr='lzw+p', res=700)
pdf('./fig/fig_03_mapvuln.pdf', wid=6.5, hei=3.25)
gp <- grid::gpar(fontsize=10, fontface='bold')
grid.arrange(
     ggplot() + theme_void() ,
     grid::textGrob(expression(bolditalic(v[1])), gp=gp),
     grid::textGrob(expression(bolditalic(v[2])), gp=gp),
     grid::textGrob(expression(bolditalic(v[3])), gp=gp),
     grid::textGrob('MWMT', rot=90, gp=gp),
     mmap(tmp, 't11', name=lab1, sizept=.5, sizetxt=.7),
     mmap(tmp, 't12', name=lab2, sizept=.5, sizetxt=.7, col='bw2'),
     mmap(tmp, 't1c', name=lab3, sizept=.5, sizetxt=.7),
     grid::textGrob('CMD', rot=90, gp=gp),
     mmap(tmp, 't21', name=lab1, sizept=.5, sizetxt=.7),
     mmap(tmp, 't22', name=lab2, sizept=.5, sizetxt=.7, col='bw2'),
     mmap(tmp, 't2c', name=lab3, sizept=.5, sizetxt=.7),
     ncol=4, widths=c(.1,1,1,1), heights=c(.07,.5,.5))
dev.off()
# Conclusion: vulnerability generally greatest in Calif, Ariz,
#     Georgia, for all 3 indices and all climate variables.
###   end plot nationwide vulnerability   ##########################


###   begin vuln by elevation and latitude   #######################
### Figure 4 - vuln by elevation and latitude
# tiff('./fig/fig_04_elev-latitude.tif', wid=6.5, hei=3.75,
#      units='in', compr='lzw+p', res=700)
pdf('./fig/fig_04_elev-latitude.pdf', wid=6.5, hei=3.75)
set_par(6)
par(oma=c(0,1,1.5,0))
isna <- is.na(tmp$elev) # first remove NA cases
### MWMT
qtmp <- data.frame(x1=tmp$elev, x2=tmp$lat, y=tmp$t11)[!isna,]
persp_q95(fit_q95(qtmp), qtmp, zlab='\nV1') #; qlab('A)')
mtext(expression(bolditalic(V[1])), 3, 1)
mtext('MWMT', 2, 5, las=3, font=2)
qtmp <- data.frame(x1=tmp$elev, x2=tmp$lat, y=tmp$t12)[!isna,]
persp_q95(fit_q95(qtmp), qtmp, zlab='\nV2') #; qlab('B)')
mtext(expression(bolditalic(V[2])), 3, 1)
qtmp <- data.frame(x1=tmp$elev, x2=tmp$lat, y=tmp$t13)[!isna,]
persp_q95(fit_q95(qtmp), qtmp, zlab='\nV3') #; qlab('C)')
mtext(expression(bolditalic(V[3])), 3, 1)
### CMD
qtmp <- data.frame(x1=tmp$elev, x2=tmp$lat, y=tmp$t21)[!isna,]
persp_q95(fit_q95(qtmp), qtmp, zlab='\nV1') #; qlab('D)')
mtext('CMD', 2, 5, las=3, font=2)
qtmp <- data.frame(x1=tmp$elev, x2=tmp$lat, y=tmp$t22)[!isna,]
persp_q95(fit_q95(qtmp), qtmp, zlab='\nV2', zlo=c(0)) #; qlab('E)')
qtmp <- data.frame(x1=tmp$elev, x2=tmp$lat, y=tmp$t23)[!isna,]
persp_q95(fit_q95(qtmp), qtmp, zlab='\nV3',zlo=c(-750)) #; qlab('F)')
dev.off()
rm(isna, qtmp)
###   end vuln by elevation and latitude   #########################


###   begin sensitivity analysis   #################################
# TIME WARN ! 30 min for 600 reps, 5 min for 100 reps, ~3 s per rep
nreps <- 100
s1 <- sens(xvec=id$mwmt, ybin=m_mwmt, perc=5, nrep=nreps)
s2 <- sens(xvec=id$cmd,  ybin=m_cmd,  perc=5, nrep=nreps)
s3 <- sens(xvec=id$mwmt, ybin=m_mwmt, perc=10, nrep=nreps)
s4 <- sens(xvec=id$cmd,  ybin=m_cmd,  perc=10, nrep=nreps)
s  <- t(sapply(list(s1,s2,s3,s4), colMeans))
dimnames(s) <- list(
     c('MWMT_05','CMD_05','MWMT_10','CMD_10'), c('V1','V2','V3'))
round(s,3)  ### result
write.csv(s, './fig/zzz_sensitivity.csv')
rm(s,s1,s2,s3,s4,m_cmd)
###   end sensitivity analysis   ###################################


###   begin future scenarios analysis ( MWMT only )   #############

#  add incremental warming....
x0.5 <- litetvi(spe, (id$mwmt+0.5), m_mwmt, na.rm=T)
x1.0 <- litetvi(spe, (id$mwmt+1.0), m_mwmt, na.rm=T)
x2.6 <- litetvi(spe, id$mwmt_rcp45_2055, m_mwmt, na.rm=T)
x3.5 <- litetvi(spe, id$mwmt_rcp85_2055, m_mwmt, na.rm=T)
rm(m_mwmt)

# temporary dataframe for plotting
tmp <- data.frame(id, t1=v1$t1, t2=v1$t2, t3=v1$t3,
                  x0.5, x1.0, x2.6, x3.5,
                  t3cut =as.factor(ifelse(v1$t3>0,'pos','neg')),
                  t3cut1=as.factor(ifelse(x0.5[,3]>0,'pos','neg')),
                  t3cut2=as.factor(ifelse(x1.0[,3]>0,'pos','neg')),
                  t3cut3=as.factor(ifelse(x2.6[,3]>0,'pos','neg')),
                  t3cut4=as.factor(ifelse(x3.5[,3]>0,'pos','neg')))
tmp <- tmp[!(is.na(x2.6[,3]) | is.na(id$mwmt)),]
nr  <- length(v1$t1)

### Figure 5 - map future warming effects
p  <- c('Current,  ','+0.5°C,  ','+1.0°C,  ','+2.5°C,  ','+3.6°C,  ')
nm <- colnames(tmp)
mm <- theme(plot.margin=unit(c(1,1,1,1),'mm'),legend.position='none')
pk <- which(nm %in% c('t1','t1.1','t1.2','t1.3','t1.4'))
t1crit <- apply(tmp[,pk], 2, function(x)
     formatC(sum(x>=50,na.rm=T)/nr*100, digits=1, format='f'))
titl <- sapply(1:5, function(i) paste0(p[i], t1crit[i], '%'))
p10 <- mmap(tmp, nm[pk[1]], '', titl[1], tmp[,pk[1]]/100, .5, .7) + mm
p11 <- mmap(tmp, nm[pk[2]], '', titl[2], tmp[,pk[2]]/100, .5, .7) + mm
p12 <- mmap(tmp, nm[pk[3]], '', titl[3], tmp[,pk[3]]/100, .5, .7) + mm
p13 <- mmap(tmp, nm[pk[4]], '', titl[4], tmp[,pk[4]]/100, .5, .7) + mm
p14 <- mmap(tmp, nm[pk[5]], '', titl[5], tmp[,pk[5]]/100, .5, .7) + mm
pk <- which(nm %in% c('t2','t2.1','t2.2','t2.3','t2.4'))
t2crit <- apply(tmp[,pk], 2, function(x)
     formatC(sum(x>=95,na.rm=T)/nr*100, digits=1, format='f'))
titl <- sapply(1:5, function(i) paste0(p[i], t2crit[i], '%'))
p20 <- mmap(tmp, nm[pk[1]], '', titl[1], tmp[,pk[1]]/100, .5, .7) + mm
p21 <- mmap(tmp, nm[pk[2]], '', titl[2], tmp[,pk[2]]/100, .5, .7) + mm
p22 <- mmap(tmp, nm[pk[3]], '', titl[3], tmp[,pk[3]]/100, .5, .7) + mm
p23 <- mmap(tmp, nm[pk[4]], '', titl[4], tmp[,pk[4]]/100, .5, .7) + mm
p24 <- mmap(tmp, nm[pk[5]], '', titl[5], tmp[,pk[5]]/100, .5, .7) + mm
pk <- which(nm %in% c('t3cut','t3cut1','t3cut2','t3cut3','t3cut4'))
t3crit <- apply(tmp[,pk], 2, function(x)
     formatC(length(which(x=='pos'))/nr*100, digits=1, format='f'))
titl <- sapply(1:5, function(i) paste0(p[i], t3crit[i], '%'))
p30 <- mmap(tmp, nm[pk[1]], '', titl[1], NA, .5, .7) + mm
p31 <- mmap(tmp, nm[pk[2]], '', titl[2], NA, .5, .7) + mm
p32 <- mmap(tmp, nm[pk[3]], '', titl[3], NA, .5, .7) + mm
p33 <- mmap(tmp, nm[pk[4]], '', titl[4], NA, .5, .7) + mm
p34 <- mmap(tmp, nm[pk[5]], '', titl[5], NA, .5, .7) + mm

# tiff('./fig/fig_05_warming.tif', wid=6.5, hei=7.5, units='in',
#      compr='lzw+p', res=500)
pdf('./fig/fig_05_warming.pdf', wid=6.5, hei=7.5)
(grid.arrange(layout_matrix = matrix(1:18, ncol=3, nrow=6, byrow=F),
              textGrob(lab1,gp=gp),p10,p11,p12,p13,p14,
              textGrob(lab2,gp=gp),p20,p21,p22,p23,p24,
              textGrob(lab3,gp=gp),p30,p31,p32,p33,p34,
              ncol=3, heights=c(.1,1,1,1,1,1)))
dev.off()


### Figure 6 - percentage of *very* vulnerable US communities
tmp <- data.frame(warm = rep(c(0,0.5,1.0,2.5,3.6),3),
                  crit = sapply(c(t1crit,t2crit,t3crit),as.numeric),
                  v    = as.factor(rep(c('tvi1','tvi2','tvi3'),ea=5)))
# tiff('./fig/fig_06_veryvuln.tif', wid=3.5, hei=3.5, units='in',
#      compr='lzw+p', res=400)
pdf('./fig/fig_06_veryvuln.pdf', wid=3.5, hei=3.5)
ggplot(tmp, aes(warm,crit,group=v)) + lims(x=c(0,4),y=c(0,35)) +
     labs(x='Warming scenario (°C)',
          y='Percentage of very vulnerable\nU.S. communities') +
     geom_point() + geom_line(aes(linetype=v)) +
     theme_classic() +
     theme(axis.text = element_text(colour='black'),
           legend.background = element_blank(),
           legend.position=c(0.4,0.8),
           legend.title=element_text(size=8, colour='black'),
           legend.text=element_text(size=8, colour='black'),
           legend.key.height=unit(0.03,'npc')) +
     scale_linetype_manual(name='Critical values:', values=1:3,
                           labels=c(' >50% species vulnerable (V1)',
                                    ' >95th thermal percentile (V2)',
                                    ' >0°C safety margin (V3)'))
dev.off()
###   end future scenarios analysis   ##############################

###     END     ####################################################