## =============================
## PLot of STS sequence objects
## =============================

plot.stslist <- function(x, tlim=NULL, sortv=NULL, 
	cpal=NULL, missing.color=NULL, ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, cex.plot=1, ...) {

	n <- nrow(x)
	seql <- ncol(x)
	statl <- alphabet(x)
	nr <- attr(x,"nr")

	if (is.null(xtlab))
		xtlab <- colnames(x)

	## Range 
	if (is.null(tlim)) tlim <- 1:10
	else if (tlim[1]==0) 
			tlim <- 1:n
	else if (max(tlim) > n) 
			tlim <- 1:n

	if (!is.null(sortv)) {
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

	if (is.null(ylab))
		ylab <- paste("Seq. ", min(tlim)," to ",max(tlim), " (n=",n,")", sortlab, sep="")

	olist <- list(...)

	ssamp <- x[tlim,]
	seqbar <- apply(ssamp, 1, seqgbar, statl=statl, seql=seql)

	## The PLot
	barplot(seqbar,col=cpal,
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

	if (is.null(yaxis) || yaxis) {
		if ("space" %in% names(olist)) sp <- olist[["space"]]
		else sp <- 0.2

		y.lab.pos <- sp+0.5

		for (p in 2:max(tlim))
			y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*sp)+(0.5+sp))

		axis(2, at=y.lab.pos, mgp=c(1.5,0.5,0), labels=tlim, tick=FALSE, cex.axis=cex.plot)
	}

}

