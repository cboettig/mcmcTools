step_fn <- function(pars, stepsizes = .02){
# Sequential random updating 
  j <- sample(1:length(pars), 1)
  pars[j] <- rnorm(1, pars[j], stepsizes)
  pars
}

## Proposal density is symmetric, so we won't need Q
Q <- function(pars, proposed){
  dnorm(pars, proposed, stepsizes=1, log=TRUE)
}

mcmc_fn <- function(pars, loglik, prior, MaxTime=1e3, stepsizes=.02, ...){
  history <- matrix(NA, nrow=MaxTime, ncol=(1+length(pars)))
  for(t in 1:MaxTime){
    Pi <- loglik(pars)+prior(pars)
    history[t,] <- c(Pi, pars)
    proposed <- step_fn(pars, stepsizes)
    alpha <- exp(loglik(proposed)+prior(proposed) - Pi)
  if (alpha > runif(1) )
      pars <- proposed
  }
  history
}


beta <- function(i, Delta_T=1){
## Heating fraction for heated chains
## Args
##  i: chain number.  i=1 corresponds to the cold chain
##  Delta_T: temperature scaling.  Defaults to 1.  Higher values explore more liberally
  1/(1+Delta_T*(i-1))
}


mcmcmc_fn <- function(pars, loglik, prior, MaxTime=1e3, indep=100, stepsizes=.02, ...){
# Metropolis Coupled Markov Chain Monte Carlo
# Args:
#   pars: a list of length n_chains, with numerics pars[[i]] that can be passed to loglik
#   loglik: a function to calculate the log-likelihood of chain i at pars[[i]], 
#   prior: a function to calculate the prior density
#   MaxTime: length of time to run the chain
#   indep: period of time for which chains should wander independently
#   step sizes of proposal distribution (can be numeric of length 1 or length pars)
# Returns:
#   chains: list containing matrix for each chain, first col is loglik + log prior prob,
#           remaining columns are fn parameters in order given in the pars[[i]]
  n_chains <- length(pars)
  n_pars <- length(pars[[1]])

  # in case we want to store the complete history.  Should have the option
  # of writing this to a file for speed? Or do that all in C...
  chains <- lapply(1:n_chains, function(i)  matrix(NA, nrow=MaxTime, ncol=(1+n_pars)) )

  # The independent intervals, lets us run chains in parallel during these periods
  Interval <- matrix(1:MaxTime, nrow=indep)

  # Note the outer time loop over intervals, and 
  # an inner time loop that can be parallelized over chains
  for(s in 1:(MaxTime/indep)){
  # Evolve chains independently for "indep" time steps
    out <- sfLapply(1:n_chains, 
            function(i){
              out <- matrix(NA, ncol=n_pars+1, nrow=indep)
              # Inner time loop
              for(t in 1:indep){
                Pi <- loglik(pars[[i]]) + prior(pars[[i]])
                out[t,] <- c(Pi, pars[[i]]) # more simply could print this to file, save mem
                proposed <- step_fn(pars[[i]], stepsizes)
                # Normal Hastings ratio weighted by temp fn beta 
                alpha <- exp( beta(i) * ( loglik(proposed)+prior(proposed) - Pi ) )
                if ( alpha  > runif(1) )
                  pars[[i]] <- proposed
              }
              out
            })
    # write to output
    for(i in 1:n_chains){
      chains[[i]][Interval[,s],] <- out[[i]]
      pars[[i]] <- out[[i]][indep,][-1] # copy parameters (but not loglik)
    }
    # This isn't quite precise bc this isn't being counted as a timestep yet!
    # Propose a chain swap every "indep" time steps beween pars_i and j
    pick <- sample(1:length(pars), 2)
    i <- pick[1]; j <- pick[2] # for convience
    R <- 
      beta(i) * ( loglik(pars[[j]]) + prior(pars[[j]]) - loglik(pars[[i]]) - prior(pars[[i]]) ) +
      beta(j) * ( loglik(pars[[i]]) + prior(pars[[i]])  - loglik(pars[[j]]) - prior(pars[[j]]) )
    ## verbose output about swaps
    # print(paste("swap chain", i, "with", j, "proposed, R =", exp(R)))
    if(exp(R) > runif(1)){
      # print("swap accepted")
      pars[[i]] <- pars[[j]] # accept the swap
    }

  }
  # Returns the full history of all chains
  chains 
}





