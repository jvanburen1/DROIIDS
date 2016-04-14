#' Prediction Function for DROIIDS Data
#'
#' @description This function produces predicted values of data analyzed from the \code{DROIIDS} function.
#'
#' @param analysis.name A character string giving the name of the analysis that was used in the \code{DROIIDS} function.
#' @param save An indicator (\code{TRUE}/\code{FALSE}) of whether or not a CSV file should be created in the current working directory
#' recording the predicted estimates. By default, the predicted estimates are saved in the global environment as \code{paste("DROIIDS.",
#' <analysis.name>, ".predicted", sep="")}.
#' @param return.predict An indicator (\code{TRUE}/\code{FALSE}) of whether or not the predicted estimates should be returned
#' and printed.
#'
#' @details
#'    Posterior means of the latent variables and covariates (if applicable) are used to find the predicted
#'    values using the data model in the DROIIDS equation:
#'
#'      Data Model: \eqn{z_{i,t}=\Phi a_{i,t} + 1_L X_{i} \beta}
#'
#'    Using the posterior means, the vector of estimated latent variables (\eqn{a_{i,t}}) at each time point within an individual is
#'    multiplied through the dimension reduction matrix (\eqn{\Phi}) to produce a single estimate for each location.
#'    If baseline covariates (\eqn{X_{i}}) are included in the analysis, the \eqn{X_{i} \beta} product will be added to each
#'    location's estimate. This produces a "fitted" value across individual, time, and location.
#'
#'    The \code{analysis.name} is the character string that was specified in the \code{DROIIDS} function.  This analysis name will
#'    link the \code{prediction} function to the iterations and output that were saved both in the global environments and in .txt files in
#'    your current working directory.  If the iteration data are not currently in the global environment, this function
#'    will automatically read in all the necessary information before calculating the predicted values.
#'
#'    \code{save = TRUE} will produce a CSV file in the current working directory with predicted estimates for individuals
#'    at all locations and time points. This provides predictions regardless of whether there are missing values in the raw
#'    data. The output file will be named \code{paste("DROIIDS_", <analysis.name>, "_Predicted_Values.csv", sep="")}.
#'
#' @export

prediction <- function(analysis.name=NULL, save=FALSE, return.predict=FALSE) {


  if(!(exists(paste("DROIIDS.", analysis.name, ".info", sep="")) & exists(paste("DROIIDS.", analysis.name, ".z", sep="")))) {
    DROIIDS.Info(analysis.name)
  }

  DROIIDS.Predict.Info <- get(paste("DROIIDS.", analysis.name, ".info", sep=""))

  DROIIDS.Z.Predict <- get(paste("DROIIDS.", analysis.name, ".z", sep=""))
  DROIIDS.Phi.Predict <- get(paste("DROIIDS.", analysis.name, ".Phi", sep=""))
  DROIIDS.posterior.est <- get(paste("DROIIDS.", analysis.name, ".post.est", sep=""))

  a.i.post.means <- DROIIDS.posterior.est[grep("a.i", rownames(DROIIDS.posterior.est)),"Mean"]
  a.matrix.est <- array(a.i.post.means, dim=c(DROIIDS.Predict.Info[["n.Phi.vectors"]],
                                              (DROIIDS.Predict.Info[["time.points"]]+1),
                                              DROIIDS.Predict.Info[["n.people"]]))

  if(DROIIDS.Predict.Info[["n.covariates"]] > 0) {
    beta.matrix.est <- DROIIDS.posterior.est[grep("beta", rownames(DROIIDS.posterior.est)),"Mean"]
    DROIIDS.x.Predict <- get(paste("DROIIDS.", analysis.name, ".x", sep=""))[,-1]
  } else{
    beta.matrix.est <- matrix(0, ncol=1, nrow=1)
    DROIIDS.x.Predict <- matrix(0, ncol=1, nrow=DROIIDS.Predict.Info[["n.people"]])
  }


  Predict.Function.Internal <- function(ID.Number.Predict) {
    a.predict.data <- a.matrix.est[,,ID.Number.Predict]
    x.predict.data <- DROIIDS.x.Predict[ID.Number.Predict,]
    t(DROIIDS.Phi.Predict %*% a.predict.data[,2:(DROIIDS.Predict.Info[["time.points"]]+1)]) + sum(x.predict.data*beta.matrix.est)
  }

  DROIIDS.Z.predict.values <- array(NA, dim=c(DROIIDS.Predict.Info[["time.points"]],
                                      DROIIDS.Predict.Info[["n.locations"]],
                                      DROIIDS.Predict.Info[["n.people"]]))
  for(i in 1:DROIIDS.Predict.Info[["n.people"]]) {
    DROIIDS.Z.predict.values[,,i] <- Predict.Function.Internal(i)
  }

  DROIIDS.Z.predict.matrix <- matrix(apply(DROIIDS.Z.predict.values, 3,
                                            function(x){t(matrix(x, nrow = dim(DROIIDS.Z.predict.values)[1],
                                                                 ncol = dim(DROIIDS.Z.predict.values)[2],
                                                                 byrow=FALSE))
                                            }), nrow = dim(DROIIDS.Z.predict.values)[1]*dim(DROIIDS.Z.predict.values)[3],
                                          ncol = dim(DROIIDS.Z.predict.values)[2], byrow=TRUE)

  DROIIDS.Z.predict.matrix <- cbind(DROIIDS.Z.Predict[,1:2], DROIIDS.Z.predict.matrix)

  assign(paste("DROIIDS.", analysis.name, ".predicted", sep=""), DROIIDS.Z.predict.matrix, envir=globalenv())

  if(save==TRUE) {
    write.csv(DROIIDS.Z.predict.matrix, file=paste(".//DROIIDS_", analysis.name, "_Predicted_Values.csv", sep=""), row.names=FALSE)
  }

  if(return.predict == TRUE) {
    return(DROIIDS.Z.predict.matrix)
  }
}
