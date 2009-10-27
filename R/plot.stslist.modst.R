## =============================
## Plot the modal state sequence
## =============================

plot.stslist.modst <- function(x, cpal=NULL,
	ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, cex.plot=1, ...) {

	seql <- ncol(x)
	statl <- attr(x,"alphabet")
	nbstat <- length(statl)
	n <- attr(x, "nbseq")

	if (is.null(cpal)) cpal <- attr(x,"cpal")

	if (is.null(xtlab)) xtlab <- colnames(x)

	if (is.null(ylab)) ylab <- paste("State freq. (n=",n,")",sep="")

	dist <- attr(x,"Distances")

	## ============================
	## Max distance for axis limits
	## ============================
	mod <- attr(x, "Frequencies")
	prof.freq <- matrix(0, nrow=nbstat, ncol=seql)

	## Preparing the matrix for plot
	for (i in 1:seql) {
		smax <- which(statl==x[,i])
		prof.freq[smax,i] <- mod[i]
	}

	## Frequency of the representative sequence
	nbrep <- attr(x,"Occurences")
	ctfreq <- round((nbrep/n)*100,1)
	txt <- paste("Modal state sequence (",nbrep," occurences, freq=", ctfreq ,"%)", sep="")

	barplot(prof.freq,
		space=0,
		## mgp=c(2.5,0.6,0),
		cex.names=cex.plot,
		ylim=c(0,1.2),
		col=cpal,
		## main=title,
		axisnames=FALSE,
		ylab=ylab,
		axes=FALSE,
		...)

	text(seql/2, 1.1, txt, 
		cex=cex.plot)

	## Plotting the x axis
	if (xaxis) 
		axis(1, at=1:seql-0.5, labels=xtlab, pos=-0.02,
		# mgp=c(3,0.5,0), 
		cex.axis=cex.plot)

	## Axis for the state frequencies
	if (yaxis)
		axis(2, at=seq(0,1.0,0.25), labels=c("0","0.25",".5","0.75","1"), 
			las=2, cex.axis=cex.plot)

}
