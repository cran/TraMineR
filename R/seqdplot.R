## ============================================
## PLOT OF THE STATES DISTRIBUTION BY TIME UNIT
## ============================================

seqdplot <- function(seqdata, group=NULL, title=NULL,
	cpal=NULL, 
	withborder=TRUE, 
	ylab, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1,
	use.layout=(!is.null(group) | withlegend!=FALSE), 
	legend.prop=NA, rows=NA, cols=NA, ...) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one")

	## ====================
	lseq <- seqdim(seqdata)[2]

	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)

	if (is.null(xtlab)) xtlab <- colnames(seqdata)
	else colnames(seqdata) <- xtlab
	
	if (is.null(cpal)) cpal <- attr(seqdata,"cpal")
	else cpal <- cpal

	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	## Message if obsolete option withborder is specified
	if (!missing(withborder)) 
		warning(" [!] option 'withborder' is obsolete, use 'border=NA' to plot without border")

	if (missing(ylab)) ylab.auto=TRUE 
	else ylab.auto=FALSE

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
		if(axes!=FALSE) axisp <- 1
		else axisp <- 0
		legpos <- NULL
	}

	## ========
	## Plotting
	## ========
	for (np in 1:nplot) {
		subdata <- seqdata[gindex[[np]],]
		nseq <- seqstatd(subdata, digits=NULL)$Frequencies

		if (nplot>1) {
			subtitle <- levels(group)[np]
			if (!is.null(title)) subtitle <- paste(title,"-",subtitle)
		} 
		else subtitle <- title

		n <- seqdim(subdata)[1]

		if (ylab.auto)	ylab <- paste("Freq. (n=",n,")",sep="")

		barplot(nseq,
			space=0,
			axes=FALSE,
			axisnames=FALSE,
			## cex.axis=cex.plot,
			## cex.names=cex.plot,
			col=cpal,
			main=subtitle,
			ylab=ylab,
			...)

		## Plotting the x axis
		if (any(np==axisp))
			axis(1, at=1:lseq-0.5, labels=xtlab, pos=-0.02, cex.axis=cex.plot)

		axis(2, cex.axis=cex.plot)

	}

	if (!is.null(legpos)) TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)
}
