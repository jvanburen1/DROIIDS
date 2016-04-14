#' Rates of change function for the DROIIDS Model
#'
#' @description Using the predicted estimates that were created through the \code{prediction} function in the DROIIDS
#' package, this function calculates the rates of change in time across individuals and locations.
#'
#' @param analysis.name A character string giving the name of the analysis that was used in the \code{DROIIDS} function.
#' @param save An indicator (\code{TRUE}/\code{FALSE}) of whether or not a CSV file should be created in the current working directory
#' recording rates of change for each individual. By default, the rates are saved in the global environment as
#' \code{paste("DROIIDS.", <analysis.name>, ".rates", sep="")}.
#' @param return.rate An indicator (\code{TRUE}/\code{FALSE}) of whether or not the rates of change should be returned
#' and printed.
#' @details
#'    A linear regression line is fit to the predicted data using time as the explanatory variable for every location
#'    across all individuals.  This provides an estimate of the linear rate of change.  Using the predicted data instead
#'    of the raw data accounts for measurement variability in the location values.  The \code{prediction} function
#'    will automatically be performed to produce the predicted data if the \code{paste("DROIIDS.", <analysis.name>, ".predicted"
#'    , sep="")} dataset does not already exist in the global environment.
#'
#'    The \code{analysis.name} is the character string that was specified in the \code{DROIIDS} function.  This analysis name will
#'    link the \code{rates.of.change} function to the iterations and output that were saved both in the global environments
#'    and in .txt files in your current working directory.
#'
#'    \code{save = TRUE} will produce a CSV file in the current working directory with rates of change for individuals
#'    at all locations. This provides rate estimates regardless of whether there are missing values in the raw data. The output file will
#'    be named \code{paste("DROIIDS_", <analysis.name>, "_Rates_of_Change.csv", sep="")}.
#'
#' @export


rates.of.change <- function(analysis.name=NULL, save=FALSE, return.rate=TRUE) {

  if(!exists(paste("DROIIDS.", analysis.name, ".predicted", sep=""))) {
    prediction(analysis.name=analysis.name, save=FALSE, return.est=FALSE)
  }

  DROIIDS.rate.info <- get(paste("DROIIDS.", analysis.name, ".info", sep=""))

  DROIIDS.predict.rate.data <- get(paste("DROIIDS.", analysis.name, ".predicted", sep=""))

  # Set up Time as a predictor
  time.predict <- unique(DROIIDS.predict.rate.data[,"time"])

  DROIIDS.ID.nums.rate <- matrix(unique(DROIIDS.predict.rate.data[,"ID"]), ncol=1)

  Rates <- cbind(DROIIDS.ID.nums.rate, t(apply(DROIIDS.ID.nums.rate, 1,
                  function(ID.rate) {
                        DROIIDS.ind.rate.data <- data.frame(DROIIDS.predict.rate.data[(DROIIDS.predict.rate.data[,"ID"] == ID.rate),3:dim(DROIIDS.predict.rate.data)[2]])

                        slope.ests <- apply(DROIIDS.ind.rate.data, 2, function(x) {lm(x ~ time.predict, DROIIDS.ind.rate.data)[["coefficients"]][2]})

                  }
                )))

  assign(paste("DROIIDS.", analysis.name, ".rates", sep=""), Rates, envir=globalenv())

  if(save==TRUE) {
    write.csv(Rates, file=paste(".//DROIIDS_", analysis.name, "_Rates_of_Change.csv", sep=""), row.names=FALSE)
  }

  if(return.rate == TRUE) {
    return(Rates)
  }

}

