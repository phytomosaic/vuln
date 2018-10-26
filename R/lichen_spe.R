#' @name lichen_spe
#' @title North American epiphytic macrolichen data
#' @aliases lichen_spe lichen_id lichen_mex
#' @docType data
#' @description
#' This data set gives information about spatial coordinates, species,
#'     and site descriptors for North American macrolichens.
#'
#' @usage
#' data(lichen_spe)
#' data(lichen_id)
#' data(lichen_mex)
#'
#' @format Multiple objects:\cr
#'     - \code{spe} dataframe, 443 lichen species (columns) at 6474
#'           sites (rows), with 75740 unique occurrences. Each cell is
#'           lichen abundance measured on an approximately logarithmic
#'           0-4 scale. Source: FIA.\cr
#'     - \code{id} dataframe, 57 site descriptors (columns) at 6474
#'           sites (rows).  Spatial coordinates are fuzzed for privacy
#'           of landowners. Source: FIA.\cr
#'     - \code{mex} dataframe, combined species and site descriptors
#'           (columns) for each of 172127 unique occurrences (rows)
#'           from herbarium records across 46343 unique locations.
#'           Source: CNALH.
#'
#' @details
#' Species data were sourced either from the US Forest Service's
#'      Forest Inventory and Analysis program (FIA 2011), or the
#'      Consortium for North American Lichen Herbaria (CNALH 2018),
#'      as collated by Smith et al. (2019).  Species should be
#'      identical among the FIA and CNALH datasets. Climate data
#'      originated from ClimateNA (Wang et al. 2016).
#'
#' @references
#' Consortium of North American Lichen Herbaria [CNALH]. 2018.
#'      Consortium of North American Lichen Herbaria.
#'      URL: http://lichenportal.org/portal/.
#'
#' Smith, R.J., S. Jovan, and B. McCune. 2019. Vulnerability of forest
#'      lichen communities to species loss under climatic warming.
#'      Unpublished manuscript.
#'
#' US Forest Service, Forest Inventory and Analysis [FIA]. 2011. Phase
#'      3 Field Guide – Lichen Communities, Version 5.1.
#'
#' Wang, T., A. Hamann, D. Spittlehouse, and C. Carroll. 2016. Locally
#'      downscaled and spatially customizable climate data for
#'      historical and future periods for North America. PLoS ONE
#'      11:e0156720.
#'
#' @keywords datasets
"lichen_spe"
