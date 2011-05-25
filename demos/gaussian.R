require(mcmcTools)

X <- rnorm(20, 2, 5)

loglik <- function(pars){
  sum( dnorm(X, mean=pars[1], sd=pars[2], log=TRUE) )
}

prior <- function(pars){
  1/pars[2]^2
}

# chains can have different or matching starting conditions, will mix anyway
 pars <- list(c(mu=1,sigma=2), c(mu=1,sigma=2))

chains <- mcmcmc_fn(pars, loglik, prior, MaxTime=1e3, indep=100,
                    stepsizes=c(.2, .1))

burnin=1e2


## plot results
layout(t(matrix(c(1,1,2,3), nrow=2)))
par(cex.lab=1.5, cex.axis=1.5)
plot(chains[[2]][-burnin, 1], type='l', col=rgb(0,0,1,.8),
     main="mixing chains", ylab="log prob density")
lines(chains[[1]][-burnin, 1], col="black")
legend("bottomright", c("normal", "heated"), col=c("black", "blue"), lty=1,cex=1.5)
plot(density((chains[[1]][-burnin, 2])), lwd=3, main="mu")
abline(v=2, col="red", lwd=2, lty=2)
plot(density((chains[[1]][-burnin, 3])), lwd=3, main="sigma")
abline(v=5, col="red", lwd=2, lty=2)

