## ================================
## PLot of the sequences frequency
## ================================

TraMineR.fplot <- function(data, np, title, tlim, cpal, pbarw, ylab, yaxis, axisp, xtlab, cex.plot, ...) {

	n <- nrow(data)
	seql <- ncol(data)
	statl <- attr(data,"alphabet")

	if (missing(ylab)) 
		ylab <- paste("Cum. % freq. (n=",n,")",sep="")

	if (is.null(tlim)) 
		tlim <- 10

	tab <- seqtab(data, format="STS")
	ndseq <- nrow(tab)
	if (ndseq<tlim) tlim <- ndseq
	tab <- tab[1:tlim,]

	seqlist <- suppressMessages(seqdecomp(row.names(tab)))

	seqbar <- apply(seqlist,1, seqgbar, seql=seql, statl=statl)

	if (pbarw==TRUE) barw=tab$Percent 
	else barw=1

	## The plot
	barplot(seqbar,col=cpal, width=barw,
		ylab=ylab,
		main=title,
		horiz=TRUE,
		axes=FALSE,
		axisnames=FALSE,
		...)
	
		## Plotting the x axis
		if (any(np==axisp)) 
			axis(1, at=1:seql-0.5, labels=xtlab, 
				# mgp=c(2.5,0.5,0), 
				cex.axis=cex.plot)

		## Plotting the y axis
		if (!pbarw) {
			## y.lab <- tab$Percent

			## y.lab.pos <- 0.7
			## for (p in 2:length(y.lab))
			##	y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*0.2)+0.7)

			ylab.pos <- c(0.2,(tlim+(tlim*0.2)))
		} 
		else { 
			## y.lab <- tab$Percent[tab$Percent>=1]

			## y.lab.pos <- (tab$Percent[1]/2)+1
			## for (p in 2:length(y.lab))
			##	y.lab.pos <- c(y.lab.pos, sum(y.lab[1:(p-1)]+0.4)+(y.lab[p]/2+0.5))
			ylab.pos <- c(0.2,(sum(barw)+(tlim*0.2*mean(barw))))
		}

		##	axis(2, at=y.lab.pos, 
		##		labels=paste(round(y.lab,1),sep=""), 
		##		tick=FALSE,
		##		mgp=c(1.5,0.5,0), 
		##		las=1, cex.axis=cex.plot)

		if (is.null(yaxis) || yaxis)
			axis(2, at=ylab.pos, 
				labels=paste(c(0, round(sum(tab$Percent),1)),"%",sep=""), 
				tick=TRUE,
				## mgp=c(1.5,1,0), 
				## las=1, 
				cex.axis=cex.plot)
}
	
