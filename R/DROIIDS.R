#' Dimension Reduction On Independent Infectious Disease Systems
#'
#' @description This function uses a Bayesian hierarchical model to analyze and pool inferences on a set of data
#'    containing correlated observations collected across time on independent subjects. This model provides predicted estimates on
#'    the data accounting for variability in the measurements. The data and process models are:
#'
#'    Data Model: \eqn{z_{i,t}=\Phi a_{i,t} + 1_L X_{i} \beta + \epsilon_{i,t}}
#'
#'    Process Model: \eqn{a_{i,t} = a_{i,t-1} + \eta_{i,t}}
#'
#'    where \eqn{z_{i,t}} is a vector of \code{L} normally distributed correlated outcome variables for individual \code{i} at
#'    time \code{t}. \eqn{\Phi} is a dimension reduction matrix created from the data and \eqn{a_{i,t}} is a vector of latent
#'    variables corresponding to the dimension reduction matrix.  \eqn{X_{i}} and \eqn{\beta} are baseline covariate data
#'    and corresponding population parameters, respectively.  Location errors and latent errors are represented by
#'    \eqn{\epsilon_{i,t}} and \eqn{\eta_{i,t}}, respectively.
#'
#' @param analysis.name A character string giving the name of the analysis.
#' @param z The set of outcome data containing ID, time, and observations. See details below for exact setup of \code{z}.
#' @param x An optional matrix of baseline predictor variables.
#' @param Phi A dimension reduction matrix created from the data.
#' @param n.Phi.vectors The number of vectors in the dimension reduction matrix to be created if Phi is not provided.
#' @param Phi.method The method to create Phi if Phi is not provided.  Options include factor analysis (\code{FA})
#' or empirical orthogonal functions (\code{EOF}).  The default method is factor analysis.
#' @param n.cores The number of computer cores to be used during the analysis.
#' @param n.chains The number of chains to be used in analysis.
#' @param seed The specified seed for generating random numbers.
#' @param burnin.iters The number of iterations to burn-in before determining the autocorrelation thinning value.
#' @param thinning.value Specifies the thinning value to use. This will overwrite the autocorrelation procedure.
#' @param autocorr.thin.level The minimum acceptable autocorrelation between the location variances, latent variances,
#' and coefficient parameters when determining a thinning value.
#' @param impute Indicator (\code{TRUE/FALSE}) of whether missing data in \code{z} should be imputed.
#' @param n.iter.chain The total number of iterations output per chain after thinning.
#' @param converge.check An indicator of whether the Gelman and Rubin univariate convergence tests should be used to
#' test for convergence between chains during the analysis.
#' @param num.iter.loops The number of iteration cycles that should be performed before forced termination.
#'
#' @details
#'    The \code{analysis.name} is the character string that will be used in all outfile names.
#'
#'    Each independent spatio-temporal dataset is stacked to form one overall data frame (\code{z}).  This structure needs to
#'    either be a "matrix" or "data.frame".  Any missed visit in \code{z} needs to be accounted for with indicators of missing
#'    values (\code{NA}).  All ID numbers should be represented the same number of times in \code{z}. The structure for the
#'    outcome data frame is as follows:
#'    \itemize{
#'      \item The first column ("\code{ID}") contains the numeric ID numbers of all subjects.
#'      \item The second column ("\code{time}") contains the specific time points for all subjects.
#'      \item The third through remaining columns of \code{z} contain the location observations.
#'    }
#'
#'    \code{x} is an optional data frame or matrix containing baseline covariates for each individual.  The first column
#'    should be the individual's ID number ("\code{ID}") and correspond to the same order of IDs as the observed data (\code{z}). By default,
#'    the DROIIDS function sets covariate usage to \code{NULL} allowing for usage of the DROIIDS model without covariates.
#'
#'    \code{Phi} is a dimension reduction matrix (a data frame or a matrix). This matrix should be created to account for variability
#'    in the data. Each row in Phi, \eqn{\Phi}\code{_l},  contains the multiplicative values for location \code{l} that will be combined with
#'    the latent variable estimates.  This matrix should have dimension (\code{L x p}) where \code{p} is both the number of dimension reduction matrix vectors
#'    and the number of latent variables estimated per person at each time point.  If the user does not provide a dimension reduction matrix,
#'    the DROIIDS function will create one using the specified number of vectors (\code{n.Phi.vectors}) and the creation method (\code{Phi.method}).
#'    By default, a factor analysis (\code{FA}) with a Varimax rotation on the first \code{n.Phi.vectors} empirical orthogonal functions is
#'    used to create Phi.  Empirical Orthogonal Functions (\code{EOF}) is a second option that can be used to create Phi.  The maximum number of
#'    vectors allowed is half of the number of locations in the analysis. Researchers recommend using no more than 0.1*L vectors in the creation
#'    of Phi where L is the number of locations.
#'
#'    The vector of location errors (\eqn{\epsilon_{i,t} = (\epsilon_{i,t,1}, ..., \epsilon_{i,t,l}, ..., \epsilon_{i,t,L})}) are assumed to be
#'    drawn from independent distributions.  Location error \eqn{\epsilon_{i,t,l}} for individual \code{i}, at time \code{t} and location \code{l}
#'    are assumed to be drawn from a normal distribution with variance \eqn{\sigma^2_{\epsilon,l}}. That is, \eqn{\epsilon_{i,t,l} ~ N(0,\sigma^2_{\epsilon,l})}
#'    for all \code{i} and \code{t}. Similarly, the vector of latent variances (\eqn{\eta_{i,t} = (\eta_{i,t,1}, ..., \eta_{i,t,k}, ..., \eta_{i,t,p})}) are
#'    assumed to be drawn from independent distributions.  The error \eqn{\eta_{i,t,k}} for the \code{k}th vector among the \code{n.Phi.vectors}
#'    latent variables is assumed to follow a normal distribution with variance \eqn{\sigma^2_{\eta,k}}. That is, \eqn{\eta_{i,t,k} ~ N(0,\sigma^2_{\eta,k})}
#'    for all \code{i} and \code{t}.
#'
#'    \code{n.cores} is the number of computer cores that will be used when finding the posterior distribution on the current machine.
#'    This option allows for parallel computing when working with multiple chains.  The default number of cores used is 3.
#'
#'    \code{n.chains} is the number of chains in the analysis used to estimate the posterior distributions. The default
#'    number of chains is 3.
#'
#'    \code{seed} specifies the random number generator.  If no seed is specified, one will be randomly chosen and saved in the output.
#'
#'    \code{burnin.iters} specifies how many iterations to burn prior to determining the thinning level.
#'
#'    After the specified burn-in, 1000 iterations are output and the autocorrelations within and between the location variances, latent variances,
#'    and covariates (if included) are calculated for every lag up to 100 on these 1000 iterations.  The maximum acceptable autocorrelation
#'    between the location variances (\eqn{\sigma^2_\epsilon}), latent variances (\eqn{\sigma^2_\eta}), and covariates (\eqn{\beta}) (if included) when determining the thinning value can be specified with
#'    \code{autocorr.thin.level}.  The lag that has every autocorrelation less than pre-specified acceptance level (\code{autocorr.thin.level})
#'    is chosen. The next increment of 5 from the chosen lag is used as the thinning value. For example, if the lag that has every autocorrelation
#'    less than \code{autocorr.thin.level} is 13, then 15 would be used as the thinning value.
#'
#'    After completion of the first set of \code{n.iter.chain} iterations after the burn-in, the iterations from each chain are imported and convergence is assessed using the Gelman and Rubin
#'    convergence diagnostic (potential scale reduction factor – PSRF) for all location variances, latent variances, and covariate parameters if \code{converge.check = TRUE}.  If at least one parameter in this group
#'    did not converge, the iteration process outputting a new set of \code{n.iter.chain} iterations is repeated using the last observed iteration values of the parameters as
#'    starting values. This iterative cycle continues until there is convergence between the chains or until it looped through a maximum of \code{num.iter.loops}
#'    times, whichever happens first. Convergence tests are not performed if the DROIIDS function is imputing data.
#'
#'    \code{impute} is an indicator (\code{TRUE}/\code{FALSE}) of whether imputation should be performed on missing data.  If there are
#'    missing data and impute = \code{FALSE}, only individuals with a complete set of data will remain in the analysis.  If \code{impute = TRUE}, the
#'    missing data are substituted using the mean of the quadratic fit for that individuals location over time.  After the burn-in and autocorrelation sections, the last observed
#'    latent variables and covariate parameter estimates are used to impute missing data using the data model. These imputed values are then fixed
#'    until a new iteration cycle occurs.
#'
#'    Vague, conjugate priors are used in the posterior estimation process. The following priors are used by the \code{DROIIDS} function:
#'    \itemize{
#'      \item Latent variable \eqn{k, k=1,..p} for individual \eqn{i} at time 0: (\eqn{a_{i,k,t=0} ~ N(0,(\sqrt{1000})^2)}
#'      \item Covariate \eqn{b, b=1,..,B}: \eqn{\beta_b ~ N(0,(\sqrt{1000})^2)}
#'      \item Location variance \eqn{l, l=1,...,L}: \eqn{\sigma^2_\epsilon ~ InverseGamma(0.001, 0.001)}
#'      \item Latent variable variance \eqn{k, k=1,...,p}: \eqn{\sigma^2_\eta ~ InverseGamma(0.001, 0.001)}
#'      \item Missing data \eqn{z_{i,t,l}} for individual \eqn{i}, time point \eqn{t}, and location \eqn{l}: \eqn{z_{i,t,l} ~ N(0,(\sqrt{1000})^2)}
#'    }
#'
#' @return A list containing information of the DROIIDS analysis is produced from this function.  Information in this list is used to create summary statistics,
#'    calculate predicted estimates and rates of changes, and plot the data.
#' @export

DROIIDS <- function(analysis.name=NULL, z=NA, x=NULL, Phi=NULL, n.Phi.vectors=NA, Phi.method="FA", n.cores=3, n.chains=3, seed=NA,
                                burnin.iters=1000, autocorr.thin.level=0.5, thinning.value=NULL, impute=FALSE, n.iter.chain=200,
                                converge.check=TRUE, num.iter.loops=5) {

  start.time <- proc.time() # Record the start time of the analysis

  library(MASS)
  library(parallel)
  library(coda)
  library(compiler)

  dir.create(".//DROIIDS_Iterations", showWarnings = FALSE)

  ##### Read in Data and Set Up Output #####
  if(!(class(z) == "data.frame")) { # Turn z into a data frame
    colnames.z <- colnames(z)
    z <- data.frame(z)
    colnames(z) <- colnames.z
    rm(colnames.z)
  }

  # Creates a name for the analysis if one is not specified
  if(!is.null(analysis.name)) {
    if(!(class(analysis.name) == "character")) {
      stop("Error: Not a valid analysis.name.  Must be a 'character'")
    }
  } else {
    analysis.name.date <- Sys.Date()
    analysis.name.time <- Sys.time()
    analysis.name.time <- paste(noquote(substr(analysis.name.time, 12,13)), noquote(substr(analysis.name.time, 15,16)), noquote(substr(analysis.name.time, 18,19)), sep="-")
    analysis.name <- paste(analysis.name.date, analysis.name.time, sep="_")
  }

  # Keeps only seeds that are positive integers
  if(!is.na(seed)) {
    if(is.numeric(seed)) {
      if((floor(seed) == seed) & seed > 0) {set.seed(seed)}
      else(stop("Not a valid seed"))
    }
    else(stop("Not a valid seed"))
  } else {
    seed <- sample(1:90000000, 1)
    set.seed(seed)
  }

  n.people <- length(unique(z$ID))
  time.points <- length(unique(z$time))
  IDs <- unique(z$ID)
  n.locations <- dim(z)[2] - 2 # Subtract off the columns for IDs and Times
  z.var.names <- names(z)[3:(n.locations+2)]


  # If the user specifies no imputation but there are missing data, only the complete cases are kept
  if((impute == FALSE) & (sum(is.na(z)) > 0)) {
    complete.IDs <- unique(z$ID)[table(z[complete.cases(z),1]) == time.points]
    z <- z[(z$ID %in% complete.IDs),]
    n.people <- length(complete.IDs)
    IDs <- complete.IDs
    if(!is.null(x)) {
      x <- x[(x$ID %in% complete.IDs),]
    }
  }

  if((impute == TRUE) & (sum(is.na(z)) == 0)) {
    impute <- FALSE # There are a complete set of data. No need to impute
  }

  if(!is.null(x)) {
    if(sum(is.na(x)) > 0) {
      stop("Error: There are missing data in X")
    }

    if(!all(unique(z[,1]) == unique(x[,1]))) {
      error.x.z <-cbind(unique(z[,1]), unique(x[,1]), unique(z[,1]) == unique(x[,1]))
      colnames(error.x.z) <- c("ID.z.order", "ID.x.order", "match")
      print(error.x.z)
      stop("Error: The IDs in z do not match the IDs in x")
    }

  }

  if(!(sum(table(z[,1]) == time.points) == n.people)) {
    stop("Error: There are a different number of time points per person. Make sure missing data are incorporated in z")
  }

  if(time.points < 5) {
    stop("Error: Not enough time points")
  } else if(time.points %in% 5:8) {
    short.time.indicator <- readline("Results may be variable with few time points. Do you want to continue? \nType 'y' for yes or 'n' to break. ")

    if(short.time.indicator == "n") {
      stop("User stopped due to few time points")
    } else if(!(short.time.indicator == "y")) {
      stop("Not a valid response")
    }
  }

  if(!is.null(Phi)) {
    n.Phi.vectors <- dim(Phi)[2]
    Phi <- as.matrix(Phi)
    Phi.creation <- "User_Specified"
  } else {
    if(is.numeric(n.Phi.vectors) & ((n.Phi.vectors>0) & (n.Phi.vectors<(n.locations/2)))) {
      if(ceiling(n.Phi.vectors) == n.Phi.vectors) {
        if(Phi.method == "FA") {
          Phi <- as.matrix(factanal(z[,3:dim(z)[2]], factors=n.Phi.vectors, rotation="varimax")[["loadings"]][1:n.locations,])
          Phi.creation <- "Created_With_Factor_Analysis"
        } else if(Phi.method == "EOF") {
          Phi <- as.matrix(princomp(z[,3:dim(z)[2]], cor=TRUE)[["loadings"]][,1:n.Phi.vectors])
          Phi.creation <- "Created_With_Empirical_Orthogonal_Functions"
        } else {
          stop("Not a valide Phi.method")
        }
      } else {
        stop("n.Phi.vectors is not a valid integer between 1 and (n.locations/2)")
      }
    } else {
      stop("n.Phi.vectors is not a valid integer between 1 and (n.locations/2)")
    }
  }


  # Turn outcome data into an array
  Z.matrix <- z[,3:(n.locations+2)] # Starts at the first location and goes to the end of the dataset (2 accounts for ID and time)
  Z.raw <- array(data=NA, dim=c(time.points, n.locations, n.people))
  for(i in 1:n.people) {
      Z.raw[,,i] <- matrix(unlist(Z.matrix[(z$ID == IDs[i]),]), nrow=time.points, ncol=n.locations, byrow=FALSE)
  }


  if(!is.null(x)) {
      X <- as.matrix(x[,-1], ncol=B)
      B <- dim(x)[2] - 1 # Subtract off the columns for IDs and Times
      x.var.names <- names(x)[2:(B+1)]
      covariates <- TRUE
      n.covariates <- B
  } else {
      X <- matrix(0, ncol=1, nrow=n.people)
      B <- 1
      covariates <- FALSE
      n.covariates <- 0
  }

  # Chooses seeds for three chains (setting a prior seed will ensure the same chain seeds are chosen)
  chain.set.seeds <- sample(1:90000000, n.chains, replace=FALSE)

  # Outputs Analysis Info
  cat(paste("analysis.name", "n.locations", "n.people", "time.points", "n.Phi.vectors", "Phi.creation",
            "n.covariates", "n.cores", "n.chains", "n.iter.chain", "seed", "impute", "burnin.iters", "autocorr.thin.level", "thinning.value",
            "n.total.iterations", "analysis.time", sep=","), "\n", file=paste(".//DROIIDS_", analysis.name, "_Info.txt", sep=""))
  cat(paste(analysis.name, n.locations, n.people, time.points, n.Phi.vectors, Phi.creation,
            n.covariates, n.cores, n.chains, n.iter.chain, seed, impute, burnin.iters, autocorr.thin.level,
            sep=","), file =paste(".//DROIIDS_", analysis.name, "_Info.txt", sep=""), append=TRUE)

  # Determine the missing value indices if we request imputation
  if(impute==TRUE) {
    Missing.Indices <- rep(NA, 3)
    for(i in 1:n.people) {
      for(l in 1:n.locations) {
        for(t in 1:time.points) {
          if(is.na(Z.raw[t,l,i])) {
            Missing.Indices <- rbind(Missing.Indices, c(t, l, i))
          }
        }
      }
    }
    Missing.Indices <- Missing.Indices[-1,]
    colnames(Missing.Indices) <- c("time", "location", "person")

  } else {Missing.Indices <- matrix(0,ncol=3,nrow=1)}




  ##### Determine autocorrelations of sigma.eps, sigma.eta, and beta #####
  burnin.autocorr <- 1
  thin.value.5 <- 1
  number.iteration.loops <- 0
  cl <- makeCluster(n.cores)
  clusterExport(cl, varlist=c("Iteration.Function", "Phi", "analysis.name", "n.people", "B", "n.locations", "Z.raw", "X",
                              "time.points", "n.Phi.vectors", "n.cores", "n.chains", "thin.value.5", "Missing.Indices", "impute", "n.iter.chain",
                              "covariates"), envir = environment())
  chain.list <- matrix(cbind(c(1:n.chains), chain.set.seeds, burnin.autocorr, burnin.iters, number.iteration.loops), ncol=5)
  parApply(cl, chain.list, 1, Iteration.Function)
  stopCluster(cl)
  burnin.autocorr <- 0

  if(is.null(thinning.value)) {
      # Reading in the First Chain to Determine Autocorrelations
      Chain1 <- read.csv(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", 1, "_Variable_Iteration_Updates.txt", sep=""), header=TRUE)

      if(covariates == TRUE) {
          Chain1 <- Chain1[-1, 2:(n.locations+n.Phi.vectors+B+1)] # Takes out the initial starting values row and keeps sigma.eps, sigma.eta, and beta. Only Chain1
          acf.values <- acf(Chain1, plot=FALSE, lag.max=100)
      } else {
        Chain1 <- Chain1[-1, 2:(n.locations+n.Phi.vectors+1)] # Takes out the initial starting values row and keeps sigma.eps, sigma.eta, and beta. Only Chain1
        acf.values <- acf(Chain1, plot=FALSE, lag.max=100)
      }

      thin.value.actual <- min(min(which((apply((abs(acf.values[["acf"]]) > autocorr.thin.level), 1, sum) == 0) == TRUE)), 100)
      thin.value.5 <- ceiling(thin.value.actual/5) * 5 # Rounds the lag up to the nearest 5.

      # Remove material read to save memory space
      rm(acf.values) ; rm(Chain1) ;
  } else {
    thin.value.actual <- thinning.value
    thin.value.5 <- thinning.value
  }

  ##### Perform Gibbs Sampling and Save Thinned Output #####
  number.iteration.loops <- 0
  num.non.convergence.1.20 <- 1 # Assign a starting value of non-convergence
  num.non.convergence.1.20.all.vars <- 1 # Assign a starting value of non-convergence

  # The loop stops for whichever event comes first: # loops > 10 OR everything converges
  while(sum((number.iteration.loops >= num.iter.loops), (num.non.convergence.1.20.all.vars == 0)) == 0) {
    number.iteration.loops <- number.iteration.loops + 1

    cl <- makeCluster(n.cores)
    # Iteration.Function reads in the last observed iteration as the starting point so we do not start from scratch
    clusterExport(cl, varlist=c("Iteration.Function", "Phi", "analysis.name", "n.people", "B", "n.locations", "Z.raw", "X",
                                "time.points", "n.Phi.vectors", "thin.value.5", "n.cores", "n.chains", "Missing.Indices", "impute", "n.iter.chain",
                                "covariates"), envir = environment())
    chain.list <- matrix(cbind(c(1:n.chains), chain.set.seeds, burnin.autocorr, 0, number.iteration.loops), ncol=5)
    parApply(cl, chain.list, 1, Iteration.Function)
    stopCluster(cl)

    if((impute == FALSE) & (converge.check == TRUE)) {
      # Read in the covariates for Sigma.Epsilon, Sigma.Eta, and Beta
      posterior.data <- rep(NA, (n.people*n.Phi.vectors*(time.points+1) + n.Phi.vectors + n.locations + 1))
      for(i in 1:n.chains) {
        read.data <- read.csv(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", i, "_Variable_Iteration_Updates.txt", sep=""), header=TRUE)
        data.combined <- read.data

        chain <- rep(i, (dim(data.combined)[1]-1))
        data.final <- cbind(chain, data.combined[-1,])
        posterior.data <- rbind(posterior.data, data.final)
      }
      posterior.data <- posterior.data[-1,]

      if(covariates == FALSE) {
          posterior.data <- posterior.data[,-grep("beta", names(posterior.data))]
      }

      # Converting Data to mcmcs
      mcmc.lists.chains <- lapply(matrix(1:n.chains,ncol=1), function(x) {as.mcmc(posterior.data[posterior.data$chain == x,c(3:dim(posterior.data)[2])])})
      mcmc.Combined.Chains <- mcmc.list(mcmc.lists.chains)

      # Testing Convergence on Sigma.Epsilon, Sigma.Eta, and Beta
      gelman.diagnostics.list <- gelman.diag(mcmc.Combined.Chains, autoburnin=FALSE, multivariate=FALSE)
      gelman.diagnostics.df <- data.frame(gelman.diagnostics.list[["psrf"]])
      num.non.convergence.1.20 <- sum(gelman.diagnostics.df[,2] > 1.20)

      num.non.convergence.1.20.all.vars <- sum(gelman.diagnostics.df[,2] > 1.20)

    }


    # If any of the Sigma.Epsilon, Sigma.Eta, and Beta fail to converge, we automatically force it to restart
    if(num.non.convergence.1.20 > 0) {
      num.non.convergence.1.20.all.vars <- num.non.convergence.1.20
    }

  }

  n.total.iterations <- burnin.iters + 1000 + number.iteration.loops * thin.value.5 * n.iter.chain

  analysis.time <- proc.time() - start.time # Determine how long all the analysis took

  cat(paste(",", thin.value.5, ",", n.total.iterations, ",", analysis.time[3], sep=""), "\n", file=paste(".//DROIIDS_", analysis.name, "_Info.txt", sep=""), append=TRUE)

  DROIIDS.Output <- list(analysis.name=analysis.name, n.locations=n.locations, n.people=n.people, time.points=time.points,
                         n.Phi.vectors=n.Phi.vectors, Phi.creation=Phi.creation, n.covariates=n.covariates, n.cores=n.cores,
                         n.chains=n.chains, n.iter.chain=n.iter.chain, seed=seed, impute=impute, burnin.iters=burnin.iters,
                         autocorr.thin.level=autocorr.thin.level, thinning.value=thin.value.5, n.total.iterations=n.total.iterations,
                         analysis.time=(analysis.time[3]*1))

  assign(paste("DROIIDS.", analysis.name, ".info", sep=""),
         list(analysis.name=analysis.name, n.locations=n.locations, n.people=n.people, time.points=time.points,
              n.Phi.vectors=n.Phi.vectors, Phi.creation=Phi.creation, n.covariates=n.covariates, n.cores=n.cores,
              n.chains=n.chains, n.iter.chain=n.iter.chain, seed=seed, impute=impute, burnin.iters=burnin.iters,
              autocorr.thin.level=autocorr.thin.level, thinning.value=thin.value.5, n.total.iterations=n.total.iterations,
              analysis.time=(analysis.time[3]*1)),
         envir=globalenv())


  # Save z, Phi, and x values into a CSV file
  write.csv(z, file=paste(".//DROIIDS_", analysis.name, "_Outcome_Data.csv", sep=""), row.names=FALSE)
  write.csv(Phi, file=paste(".//DROIIDS_", analysis.name, "_Phi.csv", sep=""), row.names=FALSE)
  if(n.covariates > 0) {
    write.csv(x, file=paste(".//DROIIDS_", analysis.name, "_Covariates.csv", sep=""), row.names=FALSE)
  }

  DROIIDS.Info(analysis.name)

  return(DROIIDS.Output)
}

