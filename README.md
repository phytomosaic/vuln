# vuln

Niche-based vulnerability of species and communities.


## Motivation

Species' occurrences indicate their realized niche, i.e., their tolerance to climate, disturbance, and other factors at those geographic locations.  Probability densities in niche space indicate habitat suitability. If niche extremes (say, >95th percentile) indicate conditions beyond which a species will decline to local extinction, then occurrences near these extremes may be considered "vulnerable". This package permits calculating single-species vulnerabilities as well as three metrics of multi-species community vulnerability.  Extensions from single to multiple environmental factors are provided. Mapped across landscapes, this will permit the targeting of locations where the greatest shifts in community composition are expected.


## Installation

Install the package from github as follows:
```r
install.packages('devtools')
devtools::install_github('phytomosaic/vuln')
```


## Load data

```r
require(vuln)
?veg
data(veg)

```

## Visualize

```r
# plot
plot(xy, pch=19, col='#00000050')

#species
plot_heatmap(spe, logbase=10)


```