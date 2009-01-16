## =============================
## Plotting individual sequences
## =============================

seqiplot <- function(seqdata, group=NULL, tlim=1:10, sortv=NULL, title=NULL, 
	cpal=NULL, missing.color=NULL,
	withborder=TRUE, 
	ylab, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1,
	use.layout=(!is.null(group) | withlegend!=FALSE), 
	legend.prop=NA, rows=NA, cols=NA, ...) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use 'seqdef' function to create one")
	
	## 
 	statl <- attr(seqdata,"alphabet")
	nr <- attr(seqdata,"nr")
	nbstat <- length(statl)
	seql <- seqdim(seqdata)[2]

	## Message if obsolete option withborder is specified
	if (!missing(withborder)) 
		warning(" [!] option 'withborder' is obsolete, use 'border=NA' to plot without border")
	
	## ================================
	## Setting color palette and labels
	## ================================
	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	if (is.null(missing.color)) missing.color <- attr(seqdata,"missing.color") 

	if (is.null(cpal)) cpal <- c(attr(seqdata,"cpal"))

	if (is.null(xtlab)) xtlab <- colnames(seqdata)

	if (missing(ylab)) ylab.auto=TRUE else ylab.auto=FALSE

	## Adding an entry for missing in the legend
	if (any(seqdata==nr)) {
		cpal <- c(cpal,missing.color)
		ltext <- c(ltext,"missing")
		statl <- c(statl,nr)
		nbstat <- nbstat+1
		}

	## ==============================
	## Preparing if group is not null
	## ==============================
	if (!is.null(group)) {
		## Eliminating the unused levels
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

	## =======
	## Ploting
	## =======
	for (np in 1:nplot) {
		subdata <- seqdata[gindex[[np]],]
		nbsub <- seqdim(subdata)[1]

		if (tlim[1]==0) 
			sublim <- 1:nbsub
		else {
			if (max(tlim) <= nbsub) 
				sublim <- tlim
			else sublim <- 1:nbsub
		}

		if (!is.null(sortv)) {
			subsort <- sortv[gindex[[np]]] 
			subdata <- subdata[order(subsort),]
			sortlab <- paste(", sorted")
		}
		else sortlab <- NULL

		if (ylab.auto)
			ylab <- paste("Seq. ", min(sublim)," to ",max(sublim), " (n=",nbsub,")", sortlab, sep="")

		ssamp <- subdata[sublim,]
		nbseq <- seqdim(ssamp)[1]
	
		seqbar <- apply(ssamp,1,seqgbar, statl=statl, seql=seql)

		if (nplot>1) {
			subtitle <- levels(group)[np]
			if (!is.null(title)) subtitle <- paste(title,"-",subtitle)
		} 
		else subtitle <- title

		## The PLot
		barplot(seqbar,col=cpal,
			ylab=ylab,
			main=subtitle,
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

	}

	## Plotting the legend
	if (!is.null(legpos)) TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)
}
