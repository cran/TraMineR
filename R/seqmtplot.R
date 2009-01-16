## ============================================
## PLOT OF THE STATES DISTRIBUTION BY TIME UNIT
## ============================================

seqmtplot <- function(seqdata, group=NULL, title=NULL,
	cpal=NULL,  
	ylab, axes="all", cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1,  
	use.layout=(!is.null(group) | withlegend!=FALSE), 
	legend.prop=NA, rows=NA, cols=NA, ...) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one")

	## ====================
	lseq <- seqdim(seqdata)[2]

	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)
	
	if (is.null(cpal)) cpal <- attr(seqdata,"cpal")

	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	if (missing(ylab)) ylab.auto=TRUE else ylab.auto=FALSE

	## ==============================
	## Preparing if group is not null
	## ==============================
	if (!is.null(group)) {
		## Eliminate the unused levels
		group <- factor(group)
		nplot <- length(levels(group))
		gindex <- vector("list",nplot)
				
		for (s in 1:nplot)
			gindex[[s]] <- which(group==levels(group)[s])

		nrplot <- ceiling(nplot/2)
	}
	else {
		nplot <- 1
		gindex <- vector("list",1)
		gindex[[1]] <- 1:seqdim(seqdata)[1]
	}

	## ===================
	## Defining the layout
	## ===================
	if (use.layout | !is.null(group) ) {
		lout <- TraMineR.setlayout(nplot, rows, cols, withlegend, axes, legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)
		axisp <- lout$axisp
		legpos <- lout$legpos
	}
	else {
		axisp <- 1
		legpos <- NULL
	}

	## ========
	## Plotting
	## ========
	for (np in 1:nplot) {
		subdata <- seqdata[gindex[[np]],]
		istatd <- suppressMessages(seqistatd(subdata))
		mtime <- apply(istatd,2,mean)

		if (nplot>1) {
			subtitle <- levels(group)[np]
			if (!is.null(title)) subtitle <- paste(title,"-",subtitle)
		} 
		else subtitle <- title

		n <- seqdim(subdata)[1]

		if (ylab.auto) ylab <- paste("Mean time (n=",n,")",sep="")

		barplot(mtime,
			## mgp=c(2.5,0.6,0),
			cex.names=cex.plot,
			cex.axis=cex.plot,
			col=cpal,
			ylim=c(0,lseq),
			main=subtitle,
			ylab=ylab,
			axes=FALSE,
			...)

		## Plotting the axes
		## axis(1, at=1:nbstat, labels=ltext, cex.axis=cex.plot)
		axis(2, at=round(seq(0,lseq,length.out=6),0), cex.axis=cex.plot)
	
	}

	## Plotting the legend
	if (!is.null(legpos)) TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)
}
