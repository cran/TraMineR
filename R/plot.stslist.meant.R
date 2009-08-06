## ====================
## Mean times plot
## ====================

plot.stslist.meant <- function(x, cpal=NULL, ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, cex.plot=1, ylim=NULL, ...) {

	n <- attr(x,"nbseq")
	seql <- length(attr(x,"xtlab"))

	if (is.null(ylab)) 
		ylab <- paste("Mean time (n=",n,")",sep="")

	if (is.null(ylim))
		ylim <- c(0,seql)

	if (is.null(xtlab))
		xtlab <- attr(x,"xtlab")

	if (is.null(cpal))
		cpal <- attr(x,"cpal")

	barplot(as.vector(x),
		## mgp=c(2.5,0.6,0),
		names.arg=if (xaxis) rownames(x) else NULL,
		cex.names=cex.plot,
		cex.axis=cex.plot,
		col=cpal,
		ylim=ylim,
		ylab=ylab,
		axes=FALSE,
		...)

	## Plotting the axes
	## axis(1, at=1:nbstat, labels=ltext, cex.axis=cex.plot)

	if (yaxis)
		axis(2, at=round(seq(0, max(ylim), length.out=6),0), cex.axis=cex.plot)
}


