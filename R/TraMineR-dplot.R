## ==============================
## PLot of the state distribution
## ==============================

TraMineR.dplot <- function(data, np, title, cpal, ylab, yaxis, axisp, xtlab, cex.plot, ...) {

	n <- nrow(data)
	seql <- ncol(data)

	nseq <- seqstatd(data, digits=NULL)$Frequencies

	if (missing(ylab)) 
		ylab <- paste("Freq. (n=",n,")",sep="")

	space=0

	barplot(nseq,
		space=space,
		axes=FALSE,
		axisnames=FALSE,
		## cex.axis=cex.plot,
		## cex.names=cex.plot,
		col=cpal,
		main=title,
		ylab=ylab,
		...)

	## Plotting the x axis
	if (any(np==axisp))
		axis(1, at=1:seql-0.5, labels=xtlab, pos=-0.02, cex.axis=cex.plot)

	if (is.null(yaxis) || yaxis)	
		axis(2, cex.axis=cex.plot)
}
