## ====================
## Mean times
## ====================

TraMineR.mtplot <- function(data, np, title, cpal, ylab, axisp, xtlab, cex.plot, ylim=NULL, ...) {

	n <- nrow(data)
	seql <- ncol(data)

	istatd <- suppressMessages(seqistatd(data))
	mtime <- apply(istatd,2,mean)

	if (missing(ylab)) 
		ylab <- paste("Mean time (n=",n,")",sep="")

	if (is.null(ylim))
		ylim <- c(0,seql)

	barplot(mtime,
		## mgp=c(2.5,0.6,0),
		cex.names=cex.plot,
		cex.axis=cex.plot,
		col=cpal,
		ylim=ylim,
		main=title,
		ylab=ylab,
		axes=FALSE,
		...)

	## Plotting the axes
	## axis(1, at=1:nbstat, labels=ltext, cex.axis=cex.plot)
	axis(2, at=round(seq(0,seql,length.out=6),0), cex.axis=cex.plot)
}

