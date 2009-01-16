## ================================
## PLot of the sequences frequency
## ================================

seqfplot <- function(seqdata, group=NULL, tlim=10, title=NULL, 
	pbarw=FALSE, 
	cpal=NULL, missing.color=NA, 
	ylab, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1, 
	use.layout=(!is.null(group) | withlegend!=FALSE), legend.prop=NA, rows=NA, cols=NA, ...) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	## Extracting some sequence characteristics
	seql <- seqdim(seqdata)[2]
	statl <- attr(seqdata,"alphabet")
	nr <- attr(seqdata,"nr")
	nbstat <- length(statl)

	## ================================
	## Setting color palette and labels
	## ================================
	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	if (is.null(missing.color)) missing.color <- attr(seqdata,"missing.color") 

	if (is.null(cpal)) cpal <- attr(seqdata,"cpal")

	if (is.null(xtlab)) xtlab <- colnames(seqdata)

	if (missing(ylab)) ylab.auto=TRUE 
	else ylab.auto=FALSE

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

	## =======
	## Ploting
	## =======
	for (np in 1:nplot) {
		subdata <- seqdata[gindex[[np]],]
		
		n <- seqdim(subdata)[1]

		tab <- seqtab(subdata)
		ndseq <- nrow(tab)
		if (ndseq<tlim) tlim <- ndseq
		tab <- tab[1:tlim,]

		seqlist <- suppressMessages(seqformat(row.names(tab),from="SPS",to="STS"))

		seqbar <- apply(seqlist,1, seqgbar, seql=seql, statl=statl)

		if (pbarw==TRUE) barw=tab$Percent 
		else barw=1

		if (nplot>1) {
			subtitle <- levels(group)[np]
			if (!is.null(title)) subtitle <- paste(title,"-",subtitle)
		} 
		else subtitle <- title

		if (ylab.auto)
			ylab <- paste("% freq. (n=",n,")",sep="")

		barplot(seqbar,col=cpal, width=barw,
			ylab=ylab,
			main=subtitle,
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
			y.lab <- tab$Percent

			y.lab.pos <- 0.7
			for (p in 2:length(y.lab))
				y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*0.2)+0.7)
			} 
		else { 
			y.lab <- tab$Percent[tab$Percent>=1]

			y.lab.pos <- (tab$Percent[1]/2)+1
			for (p in 2:length(y.lab))
				y.lab.pos <- c(y.lab.pos, sum(y.lab[1:(p-1)]+0.4)+(y.lab[p]/2+0.5))
			} 
				
		axis(2, at=y.lab.pos, 
			labels=paste(round(y.lab,1),sep=""), 
			tick=FALSE,
			mgp=c(1.5,0,0), 
			las=1, cex.axis=cex.plot)
	}	

	## Plotting the legend
	if (!is.null(legpos)) TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)
}
