#' @title Heat index
#'
#' @description
#' Heat index based on temperature and relative humidity. Useful for
#'     describing simultaneously hot and moist conditions, such as
#'     those that may limit positive carbon balance for lichens.
#'
#' @param TC vector of temperature values (in Celsius)
#'
#' @param RH vector of relative humidity values (as a percentage
#'     [0,100])
#'
#' @param ... further arguments passed to other functions.
#'
#' @return
#' vector of heat index values, in Celsius.
#'
#' @details
#' The heat index approximates \sQuote{experienced} heat.
#'
#' @examples
#' set.seed(99)
#' TC <- sort(runif(99, 5, 35)+1)
#' RH <- sort(runif(99, 50, 90))
#' par(mar=c(4,4,0,4), pty='s')
#' plot(RH, TC, ylim=c(0,80))
#' points(RH, heatindex(TC,RH), col=2)
#' axis(4, at = seq(0,80,len=5), labels = seq(0,80,len=5), col.axis=2)
#' mtext('Heat index', 4, line=2, col=2)
#'
#' @seealso
#' https://www.wpc.ncep.noaa.gov/html/heatindex_equation.shtml for the
#'     algorithm, and \link[weathermetrics]{heat.index} for a more
#'     complicated adaptation.
#'
#' @export
#' @rdname heatindex
`heatindex` <- function (TC, RH, round = 3){
     if (length(RH) != length(TC)) {
          stop('TC and RH of different lengths')
     }
     if (length(RH[!is.na(RH) & (RH > 100 | RH < 0)]) > 0) {
          RH[!is.na(RH) & (RH > 100 | RH < 0)] <- NA
          warning('RH values beyond [0,100];
                  heat index for these was set to NA')
     }
     `f` <-  function (TF = NA, RH = NA){
          if (is.na(RH) | is.na(TF)) {
               hi <- NA
          } else if (TF <= 40) {
               hi <- TF
          } else {
               alf <- 61 + ((TF - 68) * 1.2) + (RH * 0.094)
               hi <- 0.5*(alf + TF)
               if (hi > 79) {
                    hi <- (-42.379 +
                                2.04901523 * TF +
                                10.14333127 * RH -
                                0.22475541 * TF * RH -
                                6.83783 * 10^-3 * TF^2 -
                                5.481717 * 10^-2 * RH^2 +
                                1.22874 * 10^-3 * TF^2 *
                                RH + 8.5282 * 10^-4 * TF * RH^2 -
                                1.99 * 10^-6 * TF^2 * RH^2)
                    if (RH <= 13 & TF >= 80 & TF <= 112) {
                         adj1 <- (13 - RH)/4
                         adj2 <- sqrt((17 - abs(TF - 95))/17)
                         totadj <- adj1 * adj2
                         hi <- hi - totadj
                    } else if (RH > 85 & TF >= 80 & TF <= 87) {
                         adj1 <- (RH - 85)/10
                         adj2 <- (87 - TF)/5
                         totadj <- adj1 * adj2
                         hi <- hi + totadj
                    }
               }
          }
          return(hi)
     }
     TF <- (9/5) * TC + 32    # from Celsius to Fahrenheit
     HI <- mapply(f, TF = TF, RH = RH) # calc Heat Index
     HI <- (5/9) * (HI - 32)  # from Fahrenheit to Celsius
     HI
}
