# vuln

Niche-based vulnerability of species and communities.


## What

This package permits calculating single-species vulnerabilities as well as three metrics of multi-species community vulnerability.  Extensions from single to multiple environmental factors are provided. 


## Why

Species' occurrences indicate their realized niche, i.e., their tolerance to climate, disturbance, and other factors at those geographic locations.  Probability densities in niche space indicate habitat suitability.  If niche extremes (say, >95th percentile) indicate conditions beyond which a species will decline to local extinction, then occurrences near these extremes may be considered "vulnerable".  Mapped across landscapes, this will permit the targeting of locations where the greatest shifts in community composition are expected.


## Run the entire script to reproduce Smith, Jovan and McCune (2019)

Reproduce analyses from the publication.  Warning!  May take several minutes!
```r
source('./demo/script_for_publication.R')
```


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/vuln')
```


## Load data

Get some spatial data, with species and environment:
```r
require(vuln)
require(ecole)
data(braun)
spe <- braun$spe
env <- braun$env
```


## Vulnerability estimates

Calculate vulnerability values:
```r
v <- vuln(spe, y=env$bio1)  # w/o bingrid, may be 'truncated' niche
```


## Mapping

Map vulnerability values (and original environmental values) in geographic space:
```r
d <- data.frame(env, t1 = gv(v,1), t2 = gv(v,2), t3 = gv(v,3))
ecole::set_par(4)
plot(d$lat, d$lon, col=ecole::colvec(d$bio1), pch=16, main='bio1')
plot(d$lat, d$lon, col=ecole::colvec(d$t1),   pch=16, main='t1')
plot(d$lat, d$lon, col=ecole::colvec(d$t2),   pch=16, main='t2')
plot(d$lat, d$lon, col=ecole::colvec(d$t3),   pch=16, main='t3')
```


## Vulnerability estimates, permitting 'extended' niches

Alternative using potentially extended niches from, e.g., herbarium records or more extensive sampling throughout environmental space:
```r
# first, make species binary
pa  <- (spe>0)*1
dimnames(pa) <- dimnames(spe)

# next, get environmental values at each species occurrence, 
#     then collect in long format
`dematrify` <- function (taxa, thresh = -999) {
     # straight from Dave Roberts' labdsv::dematrify
     tmp <- which(taxa > thresh, arr.ind = TRUE)
     ii <- dimnames(tmp)[[1]]
     jj <- dimnames(taxa)[[2]][tmp[, 2]]
     abu <- taxa[tmp]
     ord <- order(tmp[, 1], tmp[, 2])
     res <- data.frame(ii[ord], jj[ord], abu[ord])
     names(res) <- c("site", "spp", "value")
     return(res)
}
x      <- dematrify(pa)
x$lon  <- dematrify(pa*env$lon)$value
x$lat  <- dematrify(pa*env$lat)$value
x$bio1 <- dematrify(pa*env$bio1)$value

# bin to grid
xrng <- range(x$lon)
yrng <- range(x$lat)
ybin <- bingrid(x, field = 'bio1', nr=10, nc=15, 
                xmn=xrng[1], xmx=xrng[2], ymn=yrng[1], ymx=yrng[2])

# calculate vulnerability values
v <- vuln(spe, y=env$bio1, ybin=ybin)

# mapping
d <- data.frame(env, t1 = gv(v,1), t2 = gv(v,2), t3 = gv(v,3))
ecole::set_par(4)
plot(d$lat, d$lon, col=ecole::colvec(d$bio1), pch=16, main='bio1')
plot(d$lat, d$lon, col=ecole::colvec(d$t1),   pch=16, main='t1')
plot(d$lat, d$lon, col=ecole::colvec(d$t2),   pch=16, main='t2')
plot(d$lat, d$lon, col=ecole::colvec(d$t3),   pch=16, main='t3')

```

Please contact `smithr2@oregonstate.edu` for updates.
