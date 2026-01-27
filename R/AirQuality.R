#' Air Quality Data Set
#'
#' Air quality measurements collected from a gas multisensor device deployed in an Italian city.
#'
#' @format A data frame with 9358 rows and 14 variables:
#' \describe{
#'   \item{`DateTime`}{Date and time of the measurement (\code{POSIXct}).}
#'   \item{`CO_GT`}{True hourly averaged CO concentration in \eqn{mg/m^3} (reference analyzer).}
#'   \item{`PT08_S1_CO`}{PT08.S1 (tin oxide) hourly averaged sensor response (CO targeted).}
#'   \item{`NMHC_GT`}{True hourly averaged non-methanic hydrocarbons concentration in \eqn{\mu g/m^3}.}
#'   \item{`C6H6_GT`}{True hourly averaged benzene concentration in \eqn{\mu g/m^3}.}
#'   \item{`PT08_S2_NMHC`}{PT08.S2 (titania) hourly averaged sensor response (NMHC targeted).}
#'   \item{`NOx_GT`}{True hourly averaged NOx concentration in ppb.}
#'   \item{`PT08_S3_NOx`}{PT08.S3 (tungsten oxide) hourly averaged sensor response (NOx targeted).}
#'   \item{`NO2_GT`}{True hourly averaged NO2 concentration in \eqn{\mu g/m^3}.}
#'   \item{`PT08_S4_NO2`}{PT08.S4 (tungsten oxide) hourly averaged sensor response (NO2 targeted).}
#'   \item{`PT08_S5_O3`}{PT08.S5 (indium oxide) hourly averaged sensor response (O3 targeted).}
#'   \item{`T`}{Temperature in degrees Celsius.}
#'   \item{`RH`}{Relative humidity (percent).}
#'   \item{`AH`}{Absolute humidity.}
#' }
#'
#' @source
#' UCI Machine Learning Repository:
#' \url{http://archive.ics.uci.edu/ml/datasets/Air+Quality}
#'
#' @keywords datasets
"air_quality"
