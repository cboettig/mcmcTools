## The nightmare ever-expanding summary/plot function
plot.pow <- function(pow, threshold = .95, legend = FALSE, type = "density", 
                     test_dist = TRUE, shade_power=FALSE, shade_p = FALSE,
                     show_data = TRUE, shade = TRUE, color1 = rgb(0,0,1,.5), 
                     color2 = rgb(1,0,0,.5), border1 = "blue", border2 = "red",
                     shade_aic = FALSE, print_text=TRUE, show_text = c("p"),
                     xlim = NULL, ylim = NULL, null_dist = TRUE, bw = "nrd0",
                     info_criterion = c("aic", "bic", "aicc", "threshold"), 
                     marker = c("line", "triangle"), ...){
# A plotting and analysis function for bootstrap model comparisons
# Args:
#  pow an output object of class pow
#  threshold: used for p/power illustrations and the threshold info_criterion
#  legend: show legend? (logical)
#  type: density plot or histogram?
#  test_dist: show the distribution created by simulation under model II?
#  shade_power: show shading of the test distribution reflecting power (logical)
#  shade_p: show shading for the p value? (logical)
#  show_data: show a mark for the observed lik. ratio of the models? (logical)
#  shade: shade in the distributions? (logical)
#  shade_aic: show the 
#  print_text: print various summary stats to terminal? (logical)
#  show_text
#  marker: show data falls using triangle or line?


  marker <- match.arg(marker)
	## Calculate the densities
	nd <- density(pow$null_dist, bw=bw)
	td <- density(pow$test_dist, bw=bw)
  n_null <- length(pow$null_dist)
  n_test <- length(pow$test_dist)


 n_data_pts <- function(object){
    if(is(object, "fitContinuous")){
      warning("aicc not caculated, returning aic")
      n_data <- Inf
    } else if(is(object, "browntree")){
      n_data <- object@nterm
    } else if(is(object, "hansentree")){
      n_data <- object@nterm
    }
  as.double(n_data)
  }

  ## Select the information criterion
	info_criterion = match.arg(info_criterion)
	if(info_criterion=="aic"){
		threshold_mark <-  2*(dof(pow$test) - dof(pow$null)) 
  } else if(info_criterion=="aicc"){
      k1 <- dof(pow$test)
      k2 <- dof(pow$null)
      n_data <- n_data_pts(pow$null)
      print(paste("AICc correction for", n_data, "tips"))
      threshold_mark <- 2*k1+2*k1*(k1+1)/(n_data-k1-1) - 2*k2+2*k2*(k2+1)/(n_data-k2-1) 

  } else if(info_criterion=="bic"){
      k <- pow$null@nterm
      threshold_mark <- log(k)*(dof(pow$test) - dof(pow$null)) 
	} else if(info_criterion=="threshold") {
		  threshold_tail <- sort(pow$null_dist)[ round(threshold*n_null) ]
		  threshold_mark <- threshold_tail #nd$x[tmp]
		  print(paste("threshold", threshold_mark))
	}

	## Calculate statistics FIXME (should be a seperate function of pow)
	threshold_tail <- sort(pow$null_dist)[ round(threshold*n_null) ]
	power <- sum(pow$test_dist > threshold_tail)/n_test
	lr <- -2*(loglik(pow$null) - loglik(pow$test))
	p <- 1-sum(pow$null_dist < lr)/pow$n_null
#	print(paste("power is ", power, ", p = ", p))
	# Check the reverse test (equivalent to swapping null and test)
	reverse_p <- sum(pow$test_dist < lr)/n_test
	reverse_threshold_tail <- sort(-pow$test_dist)[ round(threshold*n_test) ]
	reverse_power <- sum(-pow$null_dist > reverse_threshold_tail)/n_null
#	print(paste("reverse test power ", reverse_power, 
#              ", reverse test p = ", reverse_p))
	aic_wrong <- sum(pow$null_dist > threshold_mark)/n_null
	aic_power <- sum(pow$test_dist > threshold_mark)/n_test





	## Calculate Axis Limits
	if(is.null(xlim)) 
    xlim <- c(min(td$x, nd$x), max(td$x,nd$x))
	if(is.null(ylim)) 
    ylim <- c(min(td$y, nd$y), max(td$y,nd$y))


  ylim[2] <- 1.1*ylim[2]
  xlim[2] <- 1.1*xlim[2]


	## Density plots
	if(type != "hist"){
        if(shade_aic==TRUE | shade_power==TRUE){ 
            plottype="s"
        } else { 
            plottype ="n"
        }
		## Plot the null distribution with appropriate shading
		if(null_dist){ 
			plot(nd, xlim=xlim, ylim=ylim, type=plottype,
           col=border1, ...) 
			if(shade_p){
				shade_p <- which(nd$x > lr)
				polygon(c(lr,nd$x[shade_p]), c(0,nd$y[shade_p]),
                col=color1, border=border1)
			} else if(shade){
				polygon(nd$x, nd$y, col=color1, border=border1)
			}
		} else { 
			plot(0,0, xlim=xlim, ylim = ylim, xlab=" Likelihood Ratio",
           ...)  ## blank plot
		}

		## Plot the test distribution with appropriate shading
		if(test_dist){
			if(!null_dist) plot(td, xlim=xlim, type=plottype,
                          col=border2, ...)  ## just plot test dist 
			else lines(td, type=plottype, col=border2)
			threshold_tail <- sort(pow$null_dist)[ round(threshold*n_null) ]
			if(shade_power){
				shade_power <- which(td$x > threshold_tail)
				polygon(c(threshold_tail, td$x[shade_power]), c(0,td$y[shade_power]),
                col=color2, border=border2)
			} else if(shade){
				polygon(td$x, td$y, col=color2, border=border2)
			}
		}
		## AIC shading 
		if(shade_aic){
			shade_p <- which(nd$x > threshold_mark)
			polygon(c(threshold_mark,nd$x[shade_p]), c(0,nd$y[shade_p]), 
              col=color1, border=border1)
			if(test_dist){ # If test_dist on, shade the errors under the test too
				shade_p <- which(td$x < threshold_mark)
				polygon(c(threshold_mark,td$x[shade_p]), c(0,td$y[shade_p]), 
                col=color2, border=border2)
			}

		}

	## Plot histograms instead of density plots
	} else {
	hist(pow$null_dist, xlim=xlim, lwd=3, col=color1,
       border="white", ...)
		if(test_dist){
			hist(pow$test_dist, add=T, lwd=0, col=color2,
           border="white", ...)
		}
	}


  

  ##  Add lines to plots 
	if(show_data){ 
    if(marker=="line")
		  abline(v=lr, lwd=3, col="black", lty=2 )
    else if(marker=="triangle")
		points(lr,yshift(1), cex=1.5, col="black", pch=25, fg="black", 
            bg="black")
	}

	## add legend
	if(legend){
        if (shade==TRUE){
		    legend("topright", c("sim under Null", "sim under Test", "observed"), 
               pch=c(15,15,25), fg=c("white", "white", "black"), 
               col=c(color1, color2, "black"))
        }
        else if (shade_aic==TRUE & test_dist==TRUE){
   		    legend("topright", 
                 c( paste("False Pos (", round(aic_wrong*100,2),
                          "%)", sep=""), 
                    paste("False Neg (", round((1-aic_power)*100,2),
                          "%)", sep="")),
                    pch=c(15,15), col=c(color1, color2), bty="n")
        }
        else if (shade_aic==TRUE & test_dist==FALSE){
   		    legend("topright", paste("False Alarms (", round(aic_wrong*100,3),
                 "%)", sep=""), pch=15, col=color1, bty="n")
        }
    }
}

## a quick labling functions, find percent distance along x and y axis, starting at origin
xshift <- function(xsteps){
	deltax <- (par()$xaxp[2]-par()$xaxp[1])/100
	par()$xaxp[1]+xsteps*deltax
}
yshift <- function(ysteps){
	deltay <- (par()$yaxp[2]-par()$yaxp[1])/100
	par()$yaxp[1]+ysteps*deltay
}







