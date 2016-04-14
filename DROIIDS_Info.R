#' @export

DROIIDS.Info <- function(analysis.name=NA) {
  DROIIDS.Readin.Temp <- read.csv(paste(".//DROIIDS_", analysis.name, "_Info.txt", sep=""), header=TRUE)
  assign(paste("DROIIDS.", analysis.name, ".info", sep=""),
         list(analysis.name=analysis.name,
              n.locations=DROIIDS.Readin.Temp$n.locations,
              n.people=DROIIDS.Readin.Temp$n.people,
              time.points=DROIIDS.Readin.Temp$time.points,
              n.Phi.vectors=DROIIDS.Readin.Temp$n.Phi.vectors,
              Phi.creation=DROIIDS.Readin.Temp$Phi.creation,
              n.covariates=DROIIDS.Readin.Temp$n.covariates,
              n.cores=DROIIDS.Readin.Temp$n.cores,
              n.chains=DROIIDS.Readin.Temp$n.chains,
              n.iter.chain=DROIIDS.Readin.Temp$n.iter.chain,
              seed=DROIIDS.Readin.Temp$seed,
              impute=DROIIDS.Readin.Temp$impute,
              burnin.iters=DROIIDS.Readin.Temp$burnin.iters,
              autocorr.thin.level=DROIIDS.Readin.Temp$autocorr.thin.level,
              thinning.value=DROIIDS.Readin.Temp$thinning.value,
              n.total.iterations=DROIIDS.Readin.Temp$n.total.iterations,
              analysis.time=DROIIDS.Readin.Temp$analysis.time),
         envir=globalenv())



  posterior.data <- rep(NA, (DROIIDS.Readin.Temp$n.people*DROIIDS.Readin.Temp$n.Phi.vectors*(DROIIDS.Readin.Temp$time.points+1) +
                               DROIIDS.Readin.Temp$n.Phi.vectors + DROIIDS.Readin.Temp$n.covariates + DROIIDS.Readin.Temp$n.locations + 1))
  for(i in 1:DROIIDS.Readin.Temp$n.chains) {
    read.data <- read.csv(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", i, "_Variable_Iteration_Updates.txt", sep=""), header=TRUE)
    data.combined <- read.data

    for(j in 1:(ceiling(DROIIDS.Readin.Temp$n.people/10))) {
      read.data <- read.csv(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", i, "_a_i_Set_", j, "_Iteration_Updates.txt", sep=""), header=TRUE)
      data.combined <- cbind(data.combined, read.data)
    }

    chain <- rep(i, (dim(data.combined)[1]-1))
    data.final <- cbind(chain, data.combined[-1,])
    posterior.data <- rbind(posterior.data, data.final)
  }
  posterior.data <- posterior.data[-1,]
  posterior.data.no.latent <- posterior.data[,-c(grep("a.i.", colnames(posterior.data)))]

  mcmc.lists.chains <- lapply(matrix(1:DROIIDS.Readin.Temp$n.chains,ncol=1), function(x) {as.mcmc(posterior.data[posterior.data$chain == x,c(3:dim(posterior.data)[2])])})
  mcmc.Combined.Chains <- mcmc.list(mcmc.lists.chains) ;

  mcmc.lists.chains.no.latent <- lapply(matrix(1:DROIIDS.Readin.Temp$n.chains,ncol=1), function(x) {as.mcmc(posterior.data.no.latent[posterior.data.no.latent$chain == x,c(3:dim(posterior.data.no.latent)[2])])})
  mcmc.Combined.Chains.no.latent <- mcmc.list(mcmc.lists.chains.no.latent) ;

  post.summary.estimates.list <- summary(mcmc.Combined.Chains, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975))
  post.summary.estimates.df <- data.frame(cbind(post.summary.estimates.list[["statistics"]], post.summary.estimates.list[["quantiles"]],
                                                post.summary.estimates.list[["quantiles"]][,"97.5%"] - post.summary.estimates.list[["quantiles"]][,"2.5%"]))

  assign(paste("DROIIDS.", analysis.name, ".Iterations", sep=""), posterior.data, envir=globalenv())
  assign(paste("DROIIDS.", analysis.name, ".mcmc", sep=""), mcmc.Combined.Chains, envir=globalenv())
  assign(paste("DROIIDS.", analysis.name, ".mcmc.no.latent", sep=""), mcmc.Combined.Chains.no.latent, envir=globalenv())
  assign(paste("DROIIDS.", analysis.name, ".post.est", sep=""), post.summary.estimates.df, envir=globalenv())


  DROIIDS.z <- read.csv(paste(".//DROIIDS_", analysis.name, "_Outcome_Data.csv", sep=""), header=TRUE)
  DROIIDS.Phi <- as.matrix(read.csv(paste(".//DROIIDS_", analysis.name, "_Phi.csv", sep=""), header=TRUE), ncol=dim.Phi)
  if(DROIIDS.Readin.Temp$n.covariates > 0) {
    DROIIDS.x <- read.csv(file=paste(".//DROIIDS_", analysis.name, "_Covariates.csv", sep=""), header=TRUE)
    assign(paste("DROIIDS.", analysis.name, ".x", sep=""), DROIIDS.x, envir=globalenv())
  }
  assign(paste("DROIIDS.", analysis.name, ".z", sep=""), DROIIDS.z, envir=globalenv())
  assign(paste("DROIIDS.", analysis.name, ".Phi", sep=""), DROIIDS.Phi, envir=globalenv())


}
