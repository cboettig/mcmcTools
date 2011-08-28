# montecarlotest.R
## Dependencies: snowfall (parallel computing)


## requries loglik and getParameters to be defined functions. Here we just set up there generics.   

# loglik isn't a standard S3 generic, so we add it:
loglik <- function(x, y, ...) UseMethod("loglik")

# getParameters isn't a standard S3 generic, so we add it:
getParameters <- function(x, y, ...) UseMethod("getParameters")




## Functions compare_models and collect work to break up the 
## function montecarlotest() so that it can be deployed on arbitrary
## parallel fashion/architecture. Typical use looks like this: 
# #Parallelize comparison A,B in any fashion:  
# Asims <- sfLapply(1:nboot, function(i) compare_models(A, B))
# #and maybe some more sims using snow package
# more_Asims <- parLapply(cluster, 1:nboot, function(i) compare_models(A, B)) 
# Bsims <- parLapply(cluster, 1:nboot, function(i) compare_models(B, A))
#  #Concatentation will work if using Lapply loops only!!
# pow <- collect(c(Asims,more_Asims), Bsims, A, B)
# plot(pow)


compare_models <- function(A, B, interval=NA, ...){
# Function simulates data under model A,
#   re-estimates model A and B from this data.
# Args:
#   A: a model with methods update, simulate, loglik, getParameters
#   B: another model with these methods
#   interval: optional arg to simulate(A, interval), which may change 
#     the simulation interval.  NA if not used. 
#   ...: additional options to update function
# Returns the 2*difference in log likelihood, as 
#   well as the estimated parameters of model A 
#   and B from the simulated data.  
		if(is.na(interval)) 
      data <- simulate(A)
    else 
      data <- simulate(A, interval)
		if (is(data,"list")) 
      data <- data$rep.1 # handle data that has replicates
		A <- update(A, data, ...)
		B <- update(B, data, ...)
		lr <- -2*(loglik(A) - loglik(B)) 
		list(lr, getParameters(A), getParameters(B))
}

collect <- function(A_sim, B_sim, A, B){
  if(! is(A_sim, "list"))
    stop(paste("A_sim is in format", class(A_sim), "instead of list"))
  if(! is(B_sim, "list"))
    stop(paste("B_sim is in format", class(B_sim), "instead of list"))

	## grab the distribution of likelihood ratios
	A_dist <- sapply(A_sim, function(x) x[[1]]) 
	B_dist <- sapply(B_sim, function(x) -x[[1]]) 
  nboot <- length(A_dist)
	## grab dist of pars observed for null model under null model sims
	A_pars <- sapply(A_sim, function(x) x[[2]])
	B_pars <- sapply(B_sim, function(x) x[[3]])
	Asim_Bpars <- sapply(A_sim, function(x) x[[3]])
	Bsim_Apars <- sapply(B_sim, function(x) x[[2]])
	## format the output
	output <- list(null=A, test=B, nboot=nboot, 
                 null_dist=A_dist, test_dist=B_dist, 
                 null_par_dist = A_pars,
                 test_par_dist = B_pars, 
                 null_sim_test_pars = Asim_Bpars,
                 test_sim_null_pars = Bsim_Apars)
	class(output) <- "pow"
	output
}

montecarlotest <- function(null, test, nboot = 100, cpu = 2, 
	                         GetParNames=TRUE, times=NA, ...){
# bootstrap the likelihood ratio   
# Args:
#   null: a model with methods "simulate", "update", "loglik", "getParameters"
#   test: another such model
#   nboot: number of replicates
#   cpu: number of cpus available on the machine.  Note: if sfInit has already
#        been called and not stopped, those conditions will be used and this 
#        will be ignored.  
#   GetParNames: optionally attempt to report parameter names in the bootstrap
#   ...: options passed to the optim routine
#   times: simulations can be sampled at specified times.  If NA, will draw
#          sample points in simulating the process that match those of the 
#          original data.  If just a scalar number is given, will use the same 
#          interval, but draw that number of sample points.  If a numeric 
#          vector is given, it will sample at those time intervals.  
#
#          FIXME needs to be extended to sample at times longer/shorter 
#          than original domain. this is really a function of the simulate
#          method for the model objects, the montecarlotest function just 
#          passes this to the simulation function.  
# 
# Returns:
#   A list with the following objects:
#       null, test: the models passed in as null and test, 
#       nboot: number of replicate simulations used,
#       null_dist, test_dist: The distribution of Likelihood ratios observed
#                              when simulating under the respective model
#       null_par_dist: distribution of parameters observed for null model, 
#                      simulating under the null model,
#       test_par_dist: likewise for test model
#       null_sim_test_pars: parameters estimated for the test model over 
#                           simulations for the null model
#       test_sim_null_pars: parameters estimated fro the null model over
#                           simulations for the test model. 
#
	## are we in parallel already? If not, how many cpus do we have?
	if(cpu>1 & !sfIsRunning()){ 	
		sfInit(parallel=TRUE, cpu=cpu)
#    libs <- names(sessionInfo()[[5]])
#    for(i in libs) 
#		 sfLibrary(i) 
		sfExportAll()
	} else if(cpu<2 & !sfIsRunning()){
    sfInit()
	} 

	## generate each distribution
	null_sim <- sfSapply(1:nboot, function(i){
		if(is.na(times)) 
      data <- simulate(null)
    else 
      data <- simulate(null, times=times)
		if (is(data,"list")) data <- data$rep.1 # handle data that has replicates
		null <- update(null, data, ...)
		test <- update(test, data, ...)
		lr <- -2*(loglik(null) - loglik(test)) 
		list(lr, getParameters(null), getParameters(test))
	})
	test_sim <- sfSapply(1:nboot, function(i){
			if(is.na(times)) 
      data <- simulate(test)
    else 
      data <- simulate(test, times=times)
		if (is(data,"list")) data <- data$rep.1 # handle data that has replicates
		null <- update(null, data, ...)
		test <- update(test, data, ...)
		lr <- -2*(loglik(null) - loglik(test))
		list(lr, getParameters(null), getParameters(test))
	})




	## grab the distribution of likelihood ratios
	null_dist <- unlist(null_sim[1,]) 
	test_dist <- unlist(test_sim[1,]) 
	
	## grab dist of pars observed for null model under null model sims
	null_bootstrap_pars <- matrix( unlist(null_sim[2,]), ncol=nboot)
	if(GetParNames) rownames(null_bootstrap_pars) <- names(getParameters(null))

	## grab dist of pars observed for test model under test model sims
	test_bootstrap_pars <- matrix( unlist(test_sim[3,]), ncol=nboot)
	if(GetParNames) rownames(test_bootstrap_pars) <- names(getParameters(test))

	### Grab the cross-case bootstraps too:
	## grab dist of pars observed for test model under null model sims
	null_sim_test_pars <- matrix( unlist(null_sim[3,]), ncol=nboot)
	if(GetParNames) rownames(null_sim_test_pars) <- names(getParameters(test))

	## grab dist of pars observed for null model under test model sims
	test_sim_null_pars <- matrix( unlist(test_sim[2,]), ncol=nboot)
	if(GetParNames) rownames(test_sim_null_pars) <- names(getParameters(null))

	## format the output
	output <- list(null=null, test=test, nboot=nboot, 
                 null_dist=null_dist, test_dist=test_dist, 
                 null_par_dist = null_bootstrap_pars,
                 test_par_dist = test_bootstrap_pars, 
                 null_sim_test_pars = null_sim_test_pars,
                 test_sim_null_pars = test_sim_null_pars)
	class(output) <- "pow"
	output
}

## DOF calculation ## FIXME should be more generic 
dof <- function(object){
  if(is(object, "fitContinuous")) dof<- object[[1]]$k
  else if(is(object, "hansentree")){
    dof <- length(object@sqrt.alpha)+length(object@sigma)+
                  sum(sapply(object@theta,length))
  } else if(is(object, "browntree")) { 
    dof <- length(object@sigma)+sum(sapply(object@theta,length))
  } else {
    dof <- object$k
    if(is.null(object$k)) 
      print(paste("cannot determine degrees of freedom,i
                  please input for model"))
  }
  dof
}



## should also add the functions for summary parameter distributions



## Some generic overlap comparisons.  
## should be replaced with ROC curve tools, (seperate file) 

KL_divergence <- function(null_dist, test_dist){
	nd <- density(null_dist)
	td <- density(test_dist)
	P <- nd$y[nd$y > 0 & td$y > 0]
	Q <- td$y[nd$y > 0 & td$y > 0]
	KL_div <- P %*% log(P/Q)
	KL_div	
}

overlap <- function(pow, bw="nrd0"){
	nd <- density(pow$null_dist, bw=bw)
	td <- density(pow$test_dist, bw=bw)
	td$y %*% nd$y
}


