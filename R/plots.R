## The nightmare ever-expanding summary/plot function
plot.pow <- function(pow, threshold=.95, main="", legend=FALSE, type="density", 
                     test_dist=TRUE, shade_power=FALSE, shade_p=FALSE,
                     show_aic=FALSE, show_data=TRUE, shade=TRUE,
                     shade_aic=FALSE, print_text=TRUE, show_text = c("p"),
                     xlim=NULL, ylim=NULL, null_dist=TRUE, bw = "nrd0",
                     info_criterion=c("aic", "bic", "aicc", "threshold"), ...){

	## Calculate the densities
	nd <- density(pow$null_dist, bw=bw)
	td <- density(pow$test_dist, bw=bw)
  n_null <- length(pow$null_dist)
  n_test <- length(pow$test_dist)

  ## Select the information criterion
	info_criterion = match.arg(info_criterion)
	if(info_criterion=="aic"){
		threshold_mark <-  2*(dof(pow$test) - dof(pow$null)) 
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
	print(paste("power is ", power, ", p = ", p))
	# Check the reverse test (equivalent to swapping null and test)
	reverse_p <- sum(pow$test_dist < lr)/n_test
	reverse_threshold_tail <- sort(-pow$test_dist)[ round(threshold*n_test) ]
	reverse_power <- sum(-pow$null_dist > reverse_threshold_tail)/n_null
	print(paste("reverse test power ", reverse_power, 
              ", reverse test p = ", reverse_p))
	aic_wrong <- sum(pow$null_dist > threshold_mark)/n_null
	aic_power <- sum(pow$test_dist > threshold_mark)/n_test





	## Calculate Axis Limits
	if(is.null(xlim)) 
    xlim <- c(min(td$x, nd$x), max(td$x,nd$x))
	if(is.null(ylim)) 
    ylim <- c(min(td$y, nd$y), max(td$y,nd$y))



	## Density plots
	if(type != "hist"){
        if(shade_aic==TRUE | shade_power==TRUE){ 
            plottype="s"
        } else { 
            plottype ="n"
        }
		## Plot the null distribution with appropriate shading
		if(null_dist){ 
			plot(nd, xlim=xlim, ylim=ylim, main=main, type=plottype,
           col=rgb(0,0,1,1), ...) 
			if(shade_p){
				shade_p <- which(nd$x > pow$lr)
				polygon(c(pow$lr,nd$x[shade_p]), c(0,nd$y[shade_p]),
                col=rgb(0,0,1,.5), border=rgb(0,0,1,1))
			} else if(shade){
				polygon(nd$x, nd$y, col=rgb(0,0,1,.5), border=rgb(0,0,1,1))
			}
		} else { 
			plot(0,0, xlim=xlim, ylim = ylim, main=main, xlab=" Likelihood Ratio",
           ...)  ## blank plot
		}

		## Plot the test distribution with appropriate shading
		if(test_dist){
			if(!null_dist) plot(td, xlim=xlim, main=main, type=plottype,
                          col=rgb(1,0,0,1), ...)  ## just plot test dist 
			else lines(td, type=plottype, col=rgb(1,0,0,1))
			threshold_tail <- sort(pow$null_dist)[ round(threshold*n_null) ]
			if(shade_power){
				shade_power <- which(td$x > threshold_tail)
				polygon(c(threshold_tail, td$x[shade_power]), c(0,td$y[shade_power]),
                col=rgb(1,0,0,.5), border=rgb(1,0,0,1))
			} else if(shade){
				polygon(td$x, td$y, col=rgb(1,0,0,.5), border=rgb(1,0,0,1))
			}
		}
		## AIC shading 
		if(shade_aic){
			shade_p <- which(nd$x > threshold_mark)
			polygon(c(threshold_mark,nd$x[shade_p]), c(0,nd$y[shade_p]), 
              col=rgb(1,.5,0,.5), border=rgb(0,0,1,0))
			if(test_dist){ # If test_dist on, shade the errors under the test too
				shade_p <- which(td$x < threshold_mark)
				polygon(c(threshold_mark,td$x[shade_p]), c(0,td$y[shade_p]), 
                col=rgb(1,1,0,.5), border=rgb(1,0,0,0))
			}

		}

	## Plot histograms instead of density plots
	} else {
	hist(pow$null_dist, xlim=xlim, lwd=3, col=rgb(0,0,1,.5),
       border="white", main=main, ...)
		if(test_dist){
			hist(pow$test_dist, add=T, lwd=0, col=rgb(1,0,0,.5),
           border="white", main=main, ...)
		}
	}


  

  ##  Add lines to plots 
	if(show_data){ 
		#abline(v=pow$lr, lwd=3, col="darkred", lty=2 )
		points(lr,yshift(1), cex=1.5, col="black", pch=25, fg="black", 
            bg="black")
	}
	if(show_aic){
    #    abline(v=threshold_mark, lwd=3, col="darkgray", lty=3) 
    	points(threshold_mark,yshift(1), cex=1.5, col="red", pch=25, fg="red", bg="red")
}


	## add legend
	if(legend){
        if (shade==TRUE){
		    legend("topright", c("sim under Null", "sim under Test", "observed"), 
               pch=c(15,15,25), fg=c("white", "white", "black"), 
               col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "black"))
        }
        else if (shade_aic==TRUE & test_dist==TRUE){
   		    legend("topright", 
                 c( paste("False Alarms (", round(aic_wrong*100,3),
                          "%)", sep=""), 
                    paste("Missed Events (", round((1-aic_power)*100,3),
                          "%)", sep="")),
                    pch=c(15,15), col=c(rgb(1,.5,0,.5), rgb(1,1,0,.5)))
        }
        else if (shade_aic==TRUE & test_dist==FALSE){
   		    legend("topright", paste("False Alarms (", round(aic_wrong*100,3),
                 "%)", sep=""), pch=15, col=rgb(1,.5,0,.5))
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




