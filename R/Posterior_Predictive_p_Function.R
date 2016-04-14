#' Posterior predictive p-value function for the DROIIDS Model
#'
#' @description This function produces a posterior predictive p-value assessing goodness of fit
#' using the posterior estimates from data analyzed with the \code{DROIIDS} function.
#'
#' @param analysis.name A character string giving the name of the analysis that was used in the \code{DROIIDS} function.
#' @param seed The specified seed for generating random numbers.
#'
#'
#' @details
#'    The posterior iterations (\code{paste("DROIIDS.", <analysis.name>, ".Iterations", sep="")}) in the global environment is
#'    used to calculate a posterior predictive p-value.  The \code{posterior.p.value} function will automatically produce this
#'    dataset if it does not already exist in the global environment.
#'
#'    For each iteration in the \code{DROIIDS} function output, predicted mean outcome values are calculated using the
#'    observed iteration estimates for \eqn{a_{i,t}} and \eqn{\beta} with the data model (\eqn{z.mean_{i,t,l}=\Phi[l,] a_{i,t} +
#'    X_{i} \beta}). Potential outcome values are generated from a normal distribution with the
#'    predicted mean outcome values as the mean and respective iteration location variance, \eqn{\sigma^2_{\epsilon,l}}.
#'    That is, \eqn{z.gen_{i,t,l} ~ N(z.mean_{i,t,l}, \sigma^2_{\epsilon,l})}. This occurs for every iteration across all chains.
#'
#'    A goodness of fit calculation is calculated for the generated values, \eqn{z.gen_{i,t,l}}, and the observed values,
#'     \eqn{z.obs_{i,t,l}}, using the following equations for person \eqn{i}, time \eqn{t}, and location \eqn{l}:
#'    \deqn{fit.z.gen = \Sigma_{i} \Sigma_{t} \Sigma_{l} ((z.gen_{i,t,l} - z.mean_{i,t,l})^2 / \sigma^2_{\epsilon,l})}
#'    \deqn{fit.z.obs = \Sigma_{i} \Sigma_{t} \Sigma_{l} ((z.obs_{i,t,l} - z.mean_{i,t,l})^2 / \sigma^2_{\epsilon,l})}
#'
#'    Using these fit values, the proportion of times (fit.z.gen > fit.z) was calculated across all iterations and chains.
#'    When the model fits the data well, this proportion should be near 0.5.  Since this calculation involves observed data,
#'    the posterior predictive p-value cannot be computed when data are imputed.
#'
#'    This function will return the posterior predictive p-value in the global environment with the name
#'    \code{paste("DROIIDS.", <analysis.name>, ".posterior.p", sep="")}.
#'
#'    The \code{analysis.name} is the character string that was specified in the \code{DROIIDS} function.  This analysis name will
#'    link the \code{posterior.p.value} function to the iterations that were saved both in the global environments
#'    and in .txt files in your current working directory.  If the iteration data are not currently in the global
#'    environment, this function will automatically read in all the necessary information before calculating the posterior predictive
#'    p-value.
#'
#' @export

posterior.p.value <- function(analysis.name=NA, seed=NA) {

  if(!(exists(paste("DROIIDS.", analysis.name, ".Iterations", sep="")) & exists(paste("DROIIDS.", analysis.name, ".info", sep="")))) {
    DROIIDS.Info(analysis.name=analysis.name)
  }

  if(!is.na(seed)) {
    if(is.numeric(seed)) {
      if((floor(seed) == seed) & seed > 0) {set.seed(seed)}
      else(stop("Not a valid seed"))
    }
    else(stop("Not a valid seed"))
  }

  post.iterations <- get(paste("DROIIDS.", analysis.name, ".Iterations", sep=""))
  Phi <- as.matrix(get(paste("DROIIDS.", analysis.name, ".Phi", sep="")))
  DROIIDS.post.p.info <- get(paste("DROIIDS.", analysis.name, ".info", sep=""))
  n.Phi.vectors <- DROIIDS.post.p.info[["n.Phi.vectors"]]
  time.points <- DROIIDS.post.p.info[["time.points"]]
  time.points.plus.1 <- time.points+1
  n.people <- DROIIDS.post.p.info[["n.people"]]
  n.covariates <- DROIIDS.post.p.info[["n.covariates"]]
  n.locations <- DROIIDS.post.p.info[["n.locations"]]
  n.N.T.L <- n.people*time.points*n.locations
  z.obs <- matrix(t(get(paste("DROIIDS.", analysis.name, ".z", sep=""))[,3:(n.locations+2)]))

  if(DROIIDS.post.p.info[["impute"]] == FALSE) {

    discrims <- matrix(NA, nrow=dim(post.iterations)[1], ncol=2)
    for(i in 1:dim(post.iterations)[1]) {
      current.iter.info <- post.iterations[i,]

      sigma.eps.curr <- current.iter.info[grep("sigma.eps", colnames(current.iter.info))]
      a.i.curr <- as.numeric(current.iter.info[grep("a.i", colnames(current.iter.info), fixed=TRUE)])

      a.matrix.curr <- array(a.i.curr, dim=c(n.Phi.vectors, time.points.plus.1, n.people))

      if(n.covariates > 0) {
        beta.curr <- current.iter.info[grep("beta", colnames(current.iter.info))]
        DROIIDS.x.post.p <- get(paste("DROIIDS.", analysis.name, ".x", sep=""))[,-1]
      } else{
        beta.curr <- matrix(0, ncol=1, nrow=1)
        DROIIDS.x.post.p <- matrix(0, ncol=1, nrow=DROIIDS.post.p.info[["n.people"]])
      }

      curr.est <- array(NA, dim=c(n.locations, time.points, n.people))
      for(j in 1:n.people) {
        a.matrix.ind <- a.matrix.curr[,,j]
        curr.est[,,j] <- Phi %*% a.matrix.ind[,2:time.points.plus.1] + sum(DROIIDS.x.post.p[j,] * beta.curr)
      }
      curr.est.vec <- as.numeric(curr.est)
      curr.var.vec <- as.numeric(rep(sigma.eps.curr, (n.people*time.points)))
      zhat <- rnorm(n.N.T.L, curr.est.vec, sqrt(curr.var.vec))

      disc.zhat <- sum((zhat-curr.est.vec)^2/curr.var.vec, na.rm=TRUE)
      disc.z <- sum((z.obs-curr.est.vec)^2/curr.var.vec, na.rm=TRUE)

      discrims[i,] <- c(disc.zhat, disc.z)
    }

    posterior.p <- sum(discrims[,1] > discrims[,2])/dim(discrims)[1]

    assign(paste("DROIIDS.", analysis.name, ".posterior.p", sep=""), posterior.p, envir=globalenv())
    return(posterior.p)
  } else {
    stop("Error: Posterior predictive p-value cannot be calculated when data were imputed")
  }

}
