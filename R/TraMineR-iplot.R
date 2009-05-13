## ====================
## Sequence index plot
## ====================

TraMineR.iplot <- function(data, np, title, tlim, sortv, statl, cpal, ylab, yaxis, axisp, xtlab, cex.plot, ...) {

	n <- nrow(data)
	seql <- ncol(data)

	## Range 
	if (is.null(tlim)) tlim <- 1:10
	else if (tlim[1]==0) 
			tlim <- 1:n
	else if (max(tlim) > n) 
			tlim <- 1:n

	if (!is.null(sortv)) {
		data <- data[order(sortv),]
		sortlab <- paste(", sorted")
		}
	else sortlab <- NULL

	if (missing(ylab))
		ylab <- paste("Seq. ", min(tlim)," to ",max(tlim), " (n=",n,")", sortlab, sep="")

	ssamp <- data[tlim,]
	seqbar <- apply(ssamp, 1, seqgbar, statl=statl, seql=seql)

	## The PLot
	barplot(seqbar,col=cpal,
		ylab=ylab,
		main=title,
		horiz=TRUE,
		yaxt="n",
		axes=FALSE,
		las=1, 
		...
	)

	## Plotting the x axis
	if (any(np==axisp)) 
		axis(1, at=1:seql-0.5, labels=xtlab, 
		# mgp=c(3,0.5,0), 
		cex.axis=cex.plot)

	if (!is.null(yaxis) && yaxis) {
		y.lab.pos <- 0.7
		for (p in 2:max(tlim))
			y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*0.2)+0.7)

		axis(2, at=y.lab.pos, mgp=c(1.5,0.5,0), labels=tlim, tick=FALSE, cex.axis=cex.plot)
	}

}

