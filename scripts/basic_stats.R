
#' Extracts the time matched concntration vs effect data

#' @param data A data frame.  The data frame must conform to Regeneron standards and contain the variables: `week`, `ndy`, `concn`, `concni`, `pchg`, `actarmcd`, `ada`
#' @return A dataframe with columns `ndy`, `concn`, `effect`.
#' @export
#' @importFrom dplyr left_join group_by summarise
#' @examples
#'   inFile <- system.file("extdata/pk_pd_hysteresis.xls", package = "pkGraph")
#'   library(readxl)
#'   theData <-read_excel(inFile)
#'   timeMatchedPK(data = theData)
fun.quantile <- function(x) quantile(x,na.rm=TRUE, type=2 )    # consistent with SAS definition.
fun.mean <- function(x) mean(x,na.rm=TRUE)
fun.median <- function(x) median(x,na.rm=TRUE)
fun.range <- function(x) range(x,na.rm=TRUE)
fun.sum <- function(x) sum(x,na.rm=TRUE)

fun.min <- function(x) min(x,na.rm=TRUE)
fun.max <- function(x) max(x,na.rm=TRUE)
fun.SD <- function(x) sd(x,na.rm=TRUE)
fun.SE <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))
fun.CV <-  function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100
fun.N <- function(x) length(x)
fun.uniqN= function(x) length(unique(x))
fun.x <- function(x) x

fun.Q1 <- function(x) quantile(x, c(0.25), type=2, na.rm=TRUE)   # consistent with SAS definition.
fun.Q2 <- function(x) quantile(x, c(0.5), type=2, na.rm=TRUE)  # consistent with SAS definition.
fun.Q3 <- function(x) quantile(x, c(0.75), type=2, na.rm=TRUE)  # consistent with SAS definition.

fun.pct50 <- function(x) quantile(x, c(0.5), type=2, na.rm=TRUE)  # consistent with SAS definition.

fun.pct97p5 <- function(x) quantile(x, c(0.975), type=2, na.rm=TRUE)  # consistent with SAS definition.
fun.pct2p5 <- function(x) quantile(x, c(0.025), type=2, na.rm=TRUE)  # consistent with SAS definition.

fun.pct95 <- function(x) quantile(x, c(0.95), type=2, na.rm=TRUE)  # consistent with SAS definition.
fun.pct5 <- function(x) quantile(x, c(0.05), type=2, na.rm=TRUE)  # consistent with SAS definition.

fun.pct90 <- function(x) quantile(x, c(0.9), type=2, na.rm=TRUE)  # consistent with SAS definition.
fun.pct10 <- function(x) quantile(x, c(0.1), type=2, na.rm=TRUE)  # consistent with SAS definition.

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


#https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

# gm_mean <- function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }
# 
# #For example, to compute per-group geometric means of data$value
# gm_mean <- function(x, na.rm=TRUE){
#   #exp(tapply(log(data$value), data$group, mean))     # exp(mean(log(x)))  ==  prod(a)^(1/length(a)) 
#   exp(mean(log(x[is.finite(log(x))]),na.rm=na.rm))
#   
# }

#
#In case there is missing values in your data, this is not a rare case. you need to add one more argument. 

# CI(95%) 
library(pryr)
q97p5 <- pryr::partial(quantile, probs = 0.975, type=2, na.rm=TRUE)  # consistent with SAS definition.
q2p5  <- pryr::partial(quantile, probs = 0.025, type=2, na.rm=TRUE)  # consistent with SAS definition.

# CI(90%) 
q95 <- pryr::partial(quantile, probs = 0.95, type=2, na.rm=TRUE)  # consistent with SAS definition.
q05 <- pryr::partial(quantile, probs = 0.05, type=2, na.rm=TRUE)  # consistent with SAS definition.

# CI(80%) 
q90 <- pryr::partial(quantile, probs = 0.9, type=2, na.rm=TRUE)  # consistent with SAS definition.
q10 <- pryr::partial(quantile, probs = 0.1, type=2, na.rm=TRUE)  # consistent with SAS definition.

# 0.25
q25 <- pryr::partial(quantile, probs = 0.25, type=2, na.rm=TRUE)  # consistent with SAS definition.

# 0.5
q50 <- pryr::partial(quantile, probs = 0.50, type=2, na.rm=TRUE)  # consistent with SAS definition.

# 0.75
q75 <- pryr::partial(quantile, probs = 0.75, type=2, na.rm=TRUE)  # consistent with SAS definition.

#sd <- pryr::partial(sd, na.rm=TRUE)

#sd <- function(x) sd(x,na.rm=TRUE)
#se <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))
#cv <-  function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100




geomean <- function(x, na.rm = FALSE, trim = 0, ...)
{
  exp(mean(log(x[is.finite(log(x))], ...), na.rm = na.rm, trim = trim, ...))
}

geosd <- function(x, na.rm = FALSE, ...)
{
  exp(sd(log(x[is.finite(log(x))], ...), na.rm = na.rm, ...))
}


