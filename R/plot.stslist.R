## =============================
## PLot of STS sequence objects
## =============================

plot.stslist <- function(x, tlim=NULL, weighted=TRUE, sortv=NULL, 
	cpal=NULL, missing.color=NULL, ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, cex.plot=1, ...) {

	n <- nrow(x)
	seql <- ncol(x)
	statl <- alphabet(x)
	nr <- attr(x,"nr")

	if (is.null(xtlab))
		xtlab <- colnames(x)

	## Range 
	if (is.null(tlim)) {
		if (n>=10) tlim <- 1:10
		else tlim=1:n
	}
	else if (tlim[1]==0) 
			tlim <- 1:n
	else if (max(tlim) > n) 
			tlim <- 1:n
	
	## Sorting
	if (!is.null(sortv)) {
		if (length(sortv)!=n) {
			stop(call.=FALSE, "sortv must contain one value for each row in the sequence object")
		}

		x <- x[order(sortv),]
		sortlab <- paste(", sorted")
		}
	else sortlab <- NULL

	## 
	if (is.null(cpal))
		cpal <- attr(x,"cpal")

	## Adding an entry for missing in the legend
	if (any(x==nr)) {
		if (is.null(missing.color)) missing.color <- attr(x,"missing.color")
		cpal <- c(cpal, missing.color)
		statl <- c(statl, nr)
	}

	## Storing the optional parameters in a list
	olist <- list(...)

	ssamp <- x[tlim,]
	seqbar <- apply(ssamp, 1, seqgbar, statl=statl, seql=seql)

	## WEIGHTS
	## Weights
	weights <- attr(x, "weights")

	if (!weighted || is.null(weights)) {
		weights <- rep(1.0, nrow(x))
	}
	## Also takes into account that in unweighted sequence objects created with 
	## older TraMineR versions the weights attribute is a vector of 1
	## instead of NULL  
	if (all(weights==1)) 
		weighted <- FALSE

	if (weighted) {wlab <- "weighted "}
	else {wlab <- NULL}

	if (is.null(ylab))
		ylab <- paste(length(tlim)," seq. ", "(", wlab,"n=",sum(weights),")", 
			sortlab, sep="")

	## The PLot
	barplot(seqbar,col=cpal, width=weights,
		ylab=ylab,
		horiz=TRUE,
		yaxt="n",
		axes=FALSE,
		las=1, 
		...
	)

	## Plotting the x axis
	if (xaxis) 
		axis(1, at=1:seql-0.5, labels=xtlab, 
		# mgp=c(3,0.5,0), 
		cex.axis=cex.plot)

	## Plotting the y axis
	if (is.null(yaxis) || yaxis) {
		if ("space" %in% names(olist)) sp <- olist[["space"]]
		else sp <- 0.2

		idxmax <- length(tlim)
		
		if (!weighted) {
			y.lab.pos <- sp+0.5

			if (idxmax>1) {
				for (p in 2:idxmax) {
					y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*sp)+(0.5+sp))
				}
			}
		}
		else {
			y.lab.pos <- (weights[1]/2)+1
			sep <- sp*mean(weights)

			if (idxmax>1) {
				for (p in 2:idxmax)
					y.lab.pos <- c(y.lab.pos, sum(weights[1:p])+(p*sep)-weights[p]/2)
			}
		}

		axis(2, at=y.lab.pos, mgp=c(1.5,0.5,0), labels=tlim, tick=FALSE, cex.axis=cex.plot)
	}

}

