#' @export

Iteration.Function <- function(Parameters) {

  chain.iters <- Parameters[1]
  set.seed(Parameters[2]+10*Parameters[5])
  burnin.autocorr <- Parameters[3]
  burnin.iters <- Parameters[4]

  library(MASS)
  library(parallel)
  library(coda)
  library(compiler)


  ### Specify Hyperparameters ###
    # Hyperparameters for normal a at time 0 #
    mu.prior.a0 <- rep(0,n.Phi.vectors)
    sigma.prior.a0 <- rep(1000, n.Phi.vectors)
    sigma.prior.a0.i <- 1/sigma.prior.a0
    sigma.mu.prior.a0 <- diag(sigma.prior.a0.i, ncol=length(sigma.prior.a0.i)) %*% mu.prior.a0
    # Hyperparameters for inverse gamma sigma.eps (z variance matrix) #
    a.prior.sigma.eps <- 0.001
    b.prior.sigma.eps <- 0.001
    # Hyperparameters for inverse gamma sigma.eta (a variance matrix) #
    a.prior.sigma.eta <- 0.001
    b.prior.sigma.eta <- 0.001
    # Hyperparameters for beta
    mu.prior.beta <- c(rep(0, B))
    sigma.prior.beta <- rep(1000, B)
    sigma.prior.beta.i <- diag(1/sigma.prior.beta, ncol=length(sigma.prior.beta))
    sigma.mu.prior.beta <- sigma.prior.beta.i %*% mu.prior.beta
    # Hyperparameters for Z missing
    mu.Z.miss <- 0
    Sigma.Z.miss <- 1000
    Sigma.Z.miss.i <- 1/Sigma.Z.miss
    Sig.mu.Z.miss <- solve(Sigma.Z.miss) * mu.Z.miss
    Sig.mu.Z.miss.sqrt <- sqrt(solve(Sigma.Z.miss) * mu.Z.miss)


  ### Setup Parameter Matrices ###

    # a.i. matrix goes a.i.(person.number).(row.number).(column.number)
    a.i.var.list.setup <- expand.grid(c(1:n.Phi.vectors), c(1:(time.points+1)), c(1:n.people))
    a.i.var.list.setup <- noquote(paste0("a.i.", a.i.var.list.setup[,3], ".", a.i.var.list.setup[,1], ".", a.i.var.list.setup[,2]))
    a.i.var.list <- noquote(paste(a.i.var.list.setup, collapse=", "))

    sigma.eps.var.list.setup <- noquote(paste0("sigma.eps.", 1:n.locations))
    sigma.eps.var.list <- noquote(paste(sigma.eps.var.list.setup, collapse=", "))

    sigma.eta.var.list.setup <- noquote(paste0("sigma.eta.", 1:n.Phi.vectors))
    sigma.eta.var.list <- noquote(paste(sigma.eta.var.list.setup, collapse=", "))

    beta.var.list.setup <- paste0("beta.", 1:B)
    beta.var.list <- noquote(paste(beta.var.list.setup, collapse=", "))


  # Determining Starting Values for Iterations
    if(burnin.autocorr == 1) {
      # Specifies the first initial matrix for each parameter set
      a.i.old <- array(rnorm((n.Phi.vectors*(time.points+1)*n.people), mean=0, sd=5), dim=c(n.Phi.vectors, (time.points+1), n.people))
      sigma.eps.old <- runif(n.locations, min=0.01, max=30)
      sigma.eta.old <- runif(n.Phi.vectors, min=0.01, max=30)
      beta.old <- matrix(rnorm(B, mean=0, sd=5), ncol=1)
    }

    if(burnin.autocorr == 0) {
      # Read in last variable iteration estimates for starting values
      n.rows.txt <- dim(read.table(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_Variable_Iteration_Updates.txt", sep=""),
                                   sep=",", colClasses = c("numeric", rep("NULL", (n.locations+n.Phi.vectors+B))), header=TRUE))[1]

      current.variable.est <- read.table(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_Variable_Iteration_Updates.txt", sep=""), skip=n.rows.txt,
                                         header=FALSE, sep=",")[-1]

      for(j in 1:(ceiling(n.people/10))) {
        read.data <- as.vector(read.table(paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", j, "_Iteration_Updates.txt", sep=""), skip=n.rows.txt,
                                          header=FALSE, sep=","))
        current.variable.est <- cbind(current.variable.est, read.data)
      }

      current.variable.est <- as.matrix(current.variable.est, nrow=1)

      # Specifies the first initial matrix for each parameter set based on previous interation end point
      sigma.eps.old <- current.variable.est[1:n.locations]
      sigma.eta.old <- current.variable.est[(n.locations+1):(n.locations+n.Phi.vectors)]
      beta.old <- current.variable.est[(n.locations+n.Phi.vectors+1):(n.locations+n.Phi.vectors+B)]
      a.i.old <- array(as.vector(current.variable.est[(n.locations+n.Phi.vectors+B+1):length(current.variable.est)]), dim=c(n.Phi.vectors, (time.points+1), n.people))

    }

    a.i.new <- array(data=NA, dim=c(n.Phi.vectors, (time.points+1), n.people))
    sigma.eps.new <- NA
    sigma.eta.new <- NA
    beta.new <- NA


  # Functions to be used in estimation
    a.i.t0.func <- function(a.i.t0.data) {
      mu.a.t0 <- sigma.a.t0 %*% (sigma.eta.current.inv * (a.i.t0.data[(1:n.Phi.vectors)]) + sigma.mu.prior.a0)
      mvrnorm(1, mu=mu.a.t0, Sigma=sigma.a.t0)
    }

    a.i.t1.tT_1.func <- function(a.i.t1.tT_1.data) {
      mu.a.t <- sigma.a.t %*% (t(Phi) %*% diag(sigma.eps.current.inv) %*% (a.i.t1.tT_1.data[(1:n.locations)] -
                                                                             vec.1.z %*% a.i.t1.tT_1.data[(n.locations+n.Phi.vectors+n.Phi.vectors+1):(n.locations+n.Phi.vectors+n.Phi.vectors+B)] %*% beta.current) +
                               diag(sigma.eta.current.inv, ncol=n.Phi.vectors, nrow=n.Phi.vectors) %*% (a.i.t1.tT_1.data[(n.locations+1):(n.locations+n.Phi.vectors)]) +
                               diag(sigma.eta.current.inv, ncol=n.Phi.vectors, nrow=n.Phi.vectors) %*% (a.i.t1.tT_1.data[(n.locations+n.Phi.vectors+1):(n.locations+n.Phi.vectors+n.Phi.vectors)]))

      mvrnorm(1, mu=mu.a.t, Sigma=sigma.a.t)
    }

    a.i.tT.func <- function(a.i.tT.data) {
      mu.a.T <- sigma.a.T %*% (t(Phi) %*% diag(sigma.eps.current.inv) %*% (a.i.tT.data[(1:n.locations)] - vec.1.z %*% a.i.tT.data[(n.locations+n.Phi.vectors+1):(n.locations+n.Phi.vectors+B)] %*% beta.current) +
                               diag(sigma.eta.current.inv, ncol=n.Phi.vectors, nrow=n.Phi.vectors) %*% (a.i.tT.data[(n.locations+1):(n.locations+n.Phi.vectors)]))
      mvrnorm(1,mu=mu.a.T, Sigma=sigma.a.T)
    }

    beta.func <- function(){
      apply(sapply(1:n.people, function(i){X[i,] * sum(sigma.eps.current.inv * (t(Z.iter[,,i]) - Phi %*% a.i.new[,2:(time.points+1),i]))}), 1, sum)
    }

    Z.miss.func <- function(z.time, z.location, z.person, sigma.eps.current.inv, a.i.old) {
      sig.Z.miss <- 1/(sigma.eps.current.inv[z.location] + Sigma.Z.miss.i)
      mu.Z.miss <- sig.Z.miss * (sigma.eps.current.inv[z.location] %*% (matrix(Phi[z.location,], nrow=1) %*% matrix(a.i.old[,(z.time+1), z.person], ncol=1) + X[z.person,] %*% beta.current) + Sig.mu.Z.miss)
      rnorm(1, mean=mu.Z.miss, sd=sqrt(sig.Z.miss))
    }

  # Sets up matrices with data that will be used for the a.i functions
    vec.1.p <- matrix(1, ncol=1, nrow=n.Phi.vectors)
    vec.1.z <- matrix(1, ncol=1, nrow=n.locations)
    vec.1.T <- matrix(1, ncol=1, nrow=time.points)
    matrix.1.z.T <- matrix(1, ncol=time.points, nrow=n.locations)

    a.i.t0.data <- matrix(NA, ncol=n.Phi.vectors, nrow=n.people)
    a.i.t1.tT_1.data <- matrix(NA, ncol=(n.locations + n.Phi.vectors + n.Phi.vectors + B), nrow=n.people)
    a.i.tT.data <- matrix(NA, ncol=(n.locations + n.Phi.vectors + B), nrow=n.people)


  # Prints Initial Starting Values
    num.a.i.files <- ceiling(n.people / 10)

    # sigma.eps, sigma.eta, and beta
    cat(paste("Iteration", sigma.eps.var.list, sigma.eta.var.list, beta.var.list, sep=", "), "\n",
        file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_Variable_Iteration_Updates.txt", sep=""))
    cat(paste(0, noquote(paste(as.vector(sigma.eps.old), collapse=",")), noquote(paste(as.vector(sigma.eta.old), collapse=",")),
              noquote(paste(as.vector(beta.old), collapse=",")), sep=","), "\n",
        file =paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_Variable_Iteration_Updates.txt", sep=""), append=TRUE)


    # a.i values in increments of 10 to reduce file sizes
      n.people.per.file <- 10
      n.a.i.printout.per.file <- n.people.per.file*n.Phi.vectors*(time.points+1)
      for(a.i.iters in 1:(num.a.i.files-1)) {
      cat(paste(noquote(paste(a.i.var.list.setup[((a.i.iters-1)*n.a.i.printout.per.file + 1) : ((a.i.iters)*n.a.i.printout.per.file)], collapse=", "))), "\n",
          file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", a.i.iters, "_Iteration_Updates.txt", sep=""))
      cat(paste(noquote(paste(as.vector(a.i.old[,,((a.i.iters-1)*n.people.per.file + 1:n.people.per.file)]), collapse=", ")), sep=","),
          file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", a.i.iters, "_Iteration_Updates.txt", sep=""),
          append=TRUE)
      cat("\n", file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", a.i.iters, "_Iteration_Updates.txt", sep=""),
          append=TRUE)
      }

    # remaining initial a.i. values
      cat(paste(noquote(paste(a.i.var.list.setup[((num.a.i.files-1)*n.a.i.printout.per.file + 1) : length(a.i.var.list.setup)], collapse=", "))), "\n",
          file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", num.a.i.files, "_Iteration_Updates.txt", sep=""))
      cat(paste(noquote(paste(as.vector(a.i.old[,,((num.a.i.files-1)*n.people.per.file + 1:(n.people - (num.a.i.files-1)*n.people.per.file))]), collapse=", ")), sep=","),
          file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", num.a.i.files, "_Iteration_Updates.txt", sep=""),
          append=TRUE)
      cat("\n", file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", num.a.i.files, "_Iteration_Updates.txt", sep=""),
          append=TRUE)

  # Determine Z for set of iterations
  # Determine the number of iterations to run
    if(burnin.autocorr==1) {
        n.iterations <- burnin.iters + 1000

        Z <- Z.raw
        if(impute==TRUE) {
            for(miss.indice in 1:dim(Missing.Indices)[1]) {
              miss.person <- Missing.Indices[miss.indice,3]
              miss.time <- Missing.Indices[miss.indice,1]
              miss.location <- Missing.Indices[miss.indice,2]

              time.impute <- 1:time.points
              time.impute2 <- time.impute^2
              data.impute <- Z.raw[,miss.location,miss.person]

              fit.impute <- lm(data.impute~time.impute+time.impute2)
              mean.impute <- predict(fit.impute, newdata=data.frame(time.impute=miss.time,time.impute2=miss.time^2), type="response", interval="confidence", level=0.95, se.fit=TRUE)[["fit"]][1]
              Z[miss.time,miss.location,miss.person] <- mean.impute
            }
        }
        Z.iter <- Z
    }
    if(burnin.autocorr==0) {
        n.iterations <- n.iter.chain*thin.value.5

        sigma.eta.current.inv <- 1/sigma.eta.old # Assuming diagnol matrices (0 off diagnols)
        sigma.eps.current.inv <- 1/sigma.eps.old # Assuming diagnol matrices (0 off diagnols)
        beta.current <- beta.old

        Z.iter <- Z.raw
        if(impute==TRUE) {
            n.miss.obs <- dim(Missing.Indices)[1]

            for(miss.indices in 1:n.miss.obs) {
              miss.values <- Missing.Indices[miss.indices,]
              Z.iter[miss.values[1], miss.values[2], miss.values[3]] <- Z.miss.func(miss.values[1], miss.values[2], miss.values[3], sigma.eps.current.inv, a.i.old)
            }
        }
    }


  for(j in 1:n.iterations) {
      sigma.eta.current.inv <- 1/sigma.eta.old # Assuming diagnol matrices (0 off diagnols)
      sigma.eps.current.inv <- 1/sigma.eps.old # Assuming diagnol matrices (0 off diagnols)
      beta.current <- beta.old


      # Update a.i at t=0
      a.i.1.current <- a.i.old[,2,] # 2nd col is a=1. 1st col is a=0
      sigma.a.t0 <- diag(1/(sigma.eta.current.inv + sigma.prior.a0.i), ncol=n.Phi.vectors, nrow=n.Phi.vectors)

      a.i.t0.data[1:n.people, 1:n.Phi.vectors] <- t(a.i.1.current)
      a.i.new[,1,] <- apply(a.i.t0.data, 1, a.i.t0.func)


      # Update a.i at t for t=1:T-1
      sigma.a.t <- solve((t(Phi) %*% diag(sigma.eps.current.inv) %*% Phi + diag(sigma.eta.current.inv, ncol=n.Phi.vectors, nrow=n.Phi.vectors) + diag(sigma.eta.current.inv, ncol=n.Phi.vectors, nrow=n.Phi.vectors)))
      for(t in 1:(time.points-1)) {
          a.i.t1.tT_1.data[1:n.people, 1:n.locations] <- t(Z.iter[t,,])
          a.i.t1.tT_1.data[1:n.people, (n.locations+1):(n.locations+n.Phi.vectors)] <- t(a.i.old[,(t+2),])
          a.i.t1.tT_1.data[1:n.people, (n.locations+n.Phi.vectors+1):(n.locations+n.Phi.vectors+n.Phi.vectors)] <- t(a.i.new[,t,])
          a.i.t1.tT_1.data[1:n.people, (n.locations+n.Phi.vectors+n.Phi.vectors+1):(n.locations+n.Phi.vectors+n.Phi.vectors+B)] <- X

          a.i.new[,(t+1),] <- apply(a.i.t1.tT_1.data, 1, a.i.t1.tT_1.func)
      }


      # Update a.i at t=T
      a.i.tT.data[1:n.people, 1:n.locations] <- t(Z.iter[time.points,,])
      a.i.tT.data[1:n.people, (n.locations+1):(n.locations+n.Phi.vectors)] <- t(a.i.new[,time.points,])
      a.i.tT.data[1:n.people, (n.locations+n.Phi.vectors+1):(n.locations+n.Phi.vectors+B)] <- X

      sigma.a.T <- solve(t(Phi) %*% diag(sigma.eps.current.inv) %*% Phi + diag(sigma.eta.current.inv, ncol=n.Phi.vectors, nrow=n.Phi.vectors))
      a.i.new[,(time.points+1),] <- apply(a.i.tT.data, 1, a.i.tT.func)


      # Update beta
      if(covariates == TRUE) {
          sigma.beta <- time.points * t(X) %*% diag(rep(sum(sigma.eps.current.inv), n.people)) %*% X
          sigma.beta <- solve(sigma.beta + sigma.prior.beta.i)

          mu.beta <- beta.func()
          mu.beta <- sigma.beta %*% (mu.beta + sigma.mu.prior.beta)

          beta.new <- matrix(mvrnorm(1, mu=mu.beta, Sigma=sigma.beta), ncol=1)
      } else {
        beta.new <- matrix(0, ncol=1, nrow=1)
      }

      # Update sigma.eta
      eta.shape <- time.points*n.people/2 + a.prior.sigma.eta
      sigma.eta.new <- sapply(1:n.Phi.vectors, function(w){
          eta.rate <- sum(sapply(1:n.people, function(i){(a.i.new[w,2:(time.points+1),i] - a.i.new[w,1:time.points,i])^2})) / 2 + b.prior.sigma.eta
          1/rgamma(1, shape=eta.shape, rate=eta.rate)
      })


      # Update sigma.eps
      eps.shape <- time.points*n.people/2 + a.prior.sigma.eps
      sigma.eps.new <- sapply(1:n.locations, function(w){
          eps.rate <- sum(sapply(1:n.people, function(i){(Z.iter[,w,i] - t(Phi[w,] %*% a.i.new[,2:(time.points+1),i]) - vec.1.T %*% X[i,] %*% beta.new)^2})) / 2 + b.prior.sigma.eps
          1/rgamma(1, shape=eps.shape, rate=eps.rate)
      })




      # Thinning the data within the iterations
      if((j %% thin.value.5 == 0) & (j > burnin.iters)) {

          # Prints out updated values for sigma.eps, sigma.eta, and beta
          cat(paste(j, noquote(paste(as.vector(sigma.eps.new), collapse=",")), noquote(paste(as.vector(sigma.eta.new), collapse=",")),
                    noquote(paste(as.vector(beta.new), collapse=",")), sep=","),
              file =paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_Variable_Iteration_Updates.txt", sep=""), "\n",
              append=TRUE)

          # Prints out updated a.i values in increments of 10 to reduce file sizes
          for(a.i.iters in 1:(num.a.i.files-1)) {
            cat(paste(noquote(paste(as.vector(a.i.new[,,((a.i.iters-1)*n.people.per.file + 1:n.people.per.file)]), collapse=", ")), sep=","),
                file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", a.i.iters, "_Iteration_Updates.txt", sep=""), "\n",
                append=TRUE)
          }

          # Prints out remaining updated a.i. values
          cat(paste(noquote(paste(as.vector(a.i.new[,,((num.a.i.files-1)*n.people.per.file + 1:(n.people - (num.a.i.files-1)*n.people.per.file))]), collapse=", ")), sep=","),
              file=paste(".//DROIIDS_Iterations//DROIIDS_", analysis.name, "_Chain_", chain.iters, "_a_i_Set_", num.a.i.files, "_Iteration_Updates.txt", sep=""), "\n",
              append=TRUE)
      }

      a.i.old <- a.i.new
      sigma.eps.old <- sigma.eps.new
      sigma.eta.old <- sigma.eta.new
      beta.old <- beta.new

  } # End Gibbs sampling for a single chain

}
