## ============================================
## PLOT OF THE STATES DISTRIBUTION BY TIME UNIT
## ============================================

seqcplot <- function(seqdata, group=NULL, dist.matrix=NA, pbarw=TRUE, method="dist", 
	entropy=FALSE, mline=FALSE, fline=FALSE, 
	dist.central=TRUE, dist.center=TRUE,
	cpal=NULL, title=NULL, 
	axes="all",
	withlegend="auto", cex.legend=1, legend.prop=NA, ltext=NULL, xtlab=NULL,
	fontsize=1,
	use.layout=(!is.null(group) | withlegend!=FALSE), rows=NA, cols=NA, ...) {

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
#	if (!missing(withborder)) 
#		warning(" [!] option 'withborder' is obsolete, use 'border=NA' to plot without border")

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
		axisp <- 1
		legpos <- NULL
	}

	## ===============
	## Distance matrix
	## ===============
	if (method %in% c("dist", "mscore") && !length(dist.matrix)>1)
		dist.matrix <- seqdist(seqdata, method="LCS")

	dmax <- max(dist.matrix)
	
	## ========
	## Plotting
	## ========
	for (np in 1:nplot) {
	
		## Selecting sequences according to group
		subdata <- seqdata[gindex[[np]],]
		seql <- dim(subdata)[2]

		## Selecting distances according to group
		if (method %in% c("dist", "mscore")) 
			subdist <- dist.matrix[gindex[[np]],gindex[[np]]]
		else subdist <- NA

		prof <- seqctype(subdata, dist.matrix=subdist, method=method)
		prof.freq <- prof$Frequencies
		statd <- seqstatd(subdata)

		mod <- apply(prof.freq,2, max)

		if (!pbarw) {
			prof.freq[prof.freq>0] <- 1
		}

		## Rescaling for the graphic
		prof.freq <- prof.freq/2
		mod <- mod/2

		if (nplot>1) {
			subtitle <- levels(group)[np]
			if (!is.null(title)) subtitle <- paste(title,"-",subtitle)
		} 
		else subtitle <- title

		n <- seqdim(subdata)[1]

		barplot(prof.freq,
			space=0,
			## mgp=c(2.5,0.6,0),
			cex.names=fontsize,
			ylim=c(0,2.0),
			col=cpal,
			main=subtitle,
			axisnames=FALSE,
			ylab=paste("(N=",n,")",sep=""),
			axes=FALSE,
			...)

		if (fline) {
			freq <- statd$Frequencies
			for (j in 1:nbstat) 
				lines(1:seql-0.5, freq[j,], type="l", lwd=3.5, col=cpal[j])
		}
	
		if (mline)
			lines(1:seql-0.5, mod, type="l", lwd=3.5, col="red")

		if (entropy) lines(1:seql-0.5, statd$Entropy, type="l", lwd=3.5, col="blue")
	
		## Time axis for the sequence
		axis(1, at=1:seql-0.5, labels=xtlab, 
			pos=-0.04, 
			## mgp=c(.5,.5,0), 
			cex.axis=fontsize)

		## Axis for the state frequencies
		if (pbarw) 
			axis(2, at=seq(0,0.5,0.125), labels=paste(c(0,25,50,75,100),"%", sep=""), 
				las=2, cex.axis=fontsize)

		if (method %in% c("dist", "mscore")) {
			ctfreq <- round((length(prof$Index)/n)*100,1)
			text(seql/2, 0.6, paste("Centrotype (freq=", ctfreq ,"%)", sep=""), cex=fontsize)

			legend.text <- NULL

			if (dist.central) {
				boxplot(prof$Distances/(dmax/seql), at=1.3, add=TRUE, boxwex=0.2,
					col="cyan", 
					horizontal=TRUE, axes=F,axisnames=F)

				mean.central.dist <- round(mean(prof$Distances),1)

				legend.text <- c(legend.text, paste("(A) dist. to centrotype, mean=",mean.central.dist,sep=""))

			}
			

			if (dist.center) {
				dc <- disscenter(subdist)
				boxplot(dc/(dmax/seql), at=1.1, add=TRUE, boxwex=0.2,
					col="cyan", 
					horizontal=TRUE, axes=F,axisnames=F)

				mean.center.dist <- round(mean(dc),1)
				legend.text <- c(legend.text, paste("(B) dist. to center, mean=",mean.center.dist,sep=""))

				ctoc <- dc[prof$Index[1]]
				lines(ctoc/2, 1.1, type="p", lwd=2, col="red")
			}

			## Axis for the boxplot of distances to the central seq
			if (dist.central | dist.center) {
				axis(1, at=seq(0,seql,seql/4), labels=seq(0,dmax,dmax/4), 
					pos=1.0, mgp=c(.5,.5,0), cex.axis=fontsize)
				axis(2, at=c(1.1,1.3), labels=c("B","A"), las=2, cex.axis=fontsize)
				legend("top", legend=legend.text, 
					## title="Distance", 
					box.lty=0,
					cex=fontsize)
			}
		}

	}

	if (!is.null(legpos)) TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)
}
