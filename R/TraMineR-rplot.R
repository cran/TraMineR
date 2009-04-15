## ============================================
## PLOT A REPRESENTATIVE SEQUENCE
## ============================================

TraMineR.rplot <- function(data, np, title, method, dist.matrix, 
	## dist.method='LCS', norm=FALSE, indel=1, sm, with.miss = FALSE,
	pbarw=TRUE, entropy=FALSE, mline=FALSE, fline=FALSE, 
	dist.rep=TRUE, dist.center=TRUE, dmax=NULL,
	cpal, ylab, axisp, xtlab, cex.plot, ...) {

	seql <- ncol(data)
	statl <- attr(data,"alphabet")
	nbstat <- length(statl)

	n <- nrow(data)

	## Name of the representative sequence
	if (method=="dist") 
		ctname <- "centrotype"
	else if (method=="modseq")
		ctname <- "modal seq."
	else if (method=="mscore")
		ctname <- "max score seq."
	else if (method=="smode")
		ctname <- "modal state seq."
	else if (method=="prob")
		ctname <- "max prob seq."

	if (!is.null(dist.matrix) && (ncol(dist.matrix) != n || nrow(dist.matrix) != n))
		stop(paste("You must provide a ",n,"x",n," distance matrix",sep=""), call.=FALSE)

	## ============================
	## Max distance for axis limits
	## ============================
	prof <- seqrep(data, method=method, dist.matrix=dist.matrix)

	if (is.null(dmax)) 
		if (!is.null(dist.matrix))
			dmax <- max(dist.matrix)
		else 
			dmax <- max(prof$Distances)

	prof.freq <- prof$Frequencies
	prof.freq[is.na(prof.freq)] <- 0
	statd <- seqstatd(data)
		
	mod <- apply(prof.freq,2, max)

	if (!pbarw)
		prof.freq[prof.freq>0] <- 1
	## Rescaling for the graphic
	prof.freq <- prof.freq/2
	mod <- mod/2

	barplot(prof.freq,
		space=0,
		## mgp=c(2.5,0.6,0),
		cex.names=cex.plot,
		ylim=c(0,2.0),
		col=cpal,
		main=title,
		axisnames=FALSE,
		ylab=paste("(N=",n,")",sep=""),
		axes=FALSE,
		...)

		## Plotting a line with the state frequencies
		if (fline) {
			freq <- statd$Frequencies
			for (j in 1:nbstat) 
				lines(1:seql-0.5, freq[j,], type="l", lwd=3.5, col=cpal[j])
		}
	
		if (mline)
			lines(1:seql-0.5, mod, type="l", lwd=3.5, col="red")

		if (entropy) 
			lines(1:seql-0.5, statd$Entropy/2, type="l", lwd=3.5, col="blue")
	
		## Time axis for the sequence
		axis(1, at=1:seql-0.5, labels=xtlab, 
			pos=-0.04, 
			## mgp=c(.5,.5,0), 
			cex.axis=cex.plot)

		## Axis for the state frequencies
		if (pbarw) {
			axis(2, at=seq(0,0.5,0.125), labels=c("0","",".5","",1), 
				las=2, cex.axis=cex.plot)
			text(-12, 0.25,"State freq.", cex=cex.plot, srt=90)
		}

		## Frequency of the representative sequence
		nbrep <- length(prof$Index)
		ctfreq <- round((nbrep/n)*100,1)
		text(seql/2, 0.6, paste(ctname, " (n=",nbrep,", freq=", ctfreq ,"%)", sep=""), cex=cex.plot)

		## 
		legend.text <- NULL

		## Distance to representative sequence
		if (dist.rep) {
			boxplot(prof$Distances/(dmax/seql), at=1.3, add=TRUE, boxwex=0.2,
				col="cyan", 
				horizontal=TRUE, axes=F,axisnames=F)

			mean.central.dist <- round(mean(prof$Distances),1)

			legend.text <- c(legend.text, paste("(A) dist. to ", ctname, ", mean=", mean.central.dist, sep=""))
		}
			
		## Distance to center
		if (!is.null(dist.matrix) && dist.center) {
			dc <- disscenter(dist.matrix)
			boxplot(dc/(dmax/seql), at=1.1, add=TRUE, boxwex=0.2,
				col="cyan", 
				horizontal=TRUE, axes=F,axisnames=F)

			mean.center.dist <- round(mean(dc),1)
			legend.text <- c(legend.text, paste("(B) dist. to center, mean=",mean.center.dist,sep=""))

			ctoc <- dc[prof$Index[1]]
			lines(ctoc/2, 1.1, type="p", lwd=2, col="red")
		}

		## Axis for the boxplot of distances to the central seq
		if (dist.rep | dist.center) {
			axis(1, at=seq(0,seql,seql/4), 
				labels=round(seq(0,dmax,dmax/4),1), 
				pos=1.0, mgp=c(.5,.5,0), 
				cex.axis=cex.plot)
			axis(2, at=c(1.1,1.3), 
				labels=c("B","A"), 
				las=2, 
				cex.axis=cex.plot)
			legend("top", legend=legend.text, 
				## title="Distance", 
				box.lty=0,
				cex=cex.plot)
		}
}

