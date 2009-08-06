## ============================================
## PLOT A REPRESENTATIVE SEQUENCE
## ============================================

plot.stslist.rep <- function(x, cpal=NULL,
	pbarw=TRUE, dmax=NULL,
	ylab=NULL, xaxis=TRUE, xtlab=NULL, cex.plot=1, ...) {

	## Extracting attributes
	n <- attr(x,"nbseq")
	
	if (is.null(xtlab))
		xtlab <- colnames(x)
	seql <- length(xtlab)
	statl <- attr(x,"alphabet")
	nbstat <- length(statl)
	criterion <- attr(x,"criterion")
	Quality <- attr(x,"Quality")
	Distances <- attr(x,"Distances")
	freq <- attr(x,"Frequencies")

	if (is.null(dmax)) dmax <- attr(x,"dmax")

	if (is.null(cpal)) cpal <- attr(x,"cpal")

	## Name of the representative sequence
	if (criterion=="dist") ctname <- "centrality"
	else if (criterion=="freq") ctname <- "frequency"
	else if (criterion=="mscore") ctname <- "rep. score"
	else if (criterion=="smode") ctname <- "modal state seq."
	else if (criterion=="prob") ctname <- "probability"
	else if (criterion=="density") ctname <- "density"

	## ============================
	## Max distance for axis limits
	## ============================
	nbrep <- nrow(x)

	## Scaling factor for the distance summaries
	dist.scaling <- dmax/seql

	## space between text lines
	vspace <- 0.1	

	## 
	if (pbarw) 
		barw <- Quality$pctcapt[1:nbrep]/100
	else barw=1

	seqbar <- apply(x, 1, seqgbar, seql=seql, statl=statl)

	## The plot
	ymax <- 2.5

	barplot(seqbar,col=cpal, width=barw,
		ylab=paste("(N=",n,")",sep=""),
		xlim=c(-2,seql),
		ylim=c(0,ymax),
		horiz=TRUE,
		axes=FALSE,
		axisnames=FALSE,
		...)

	## Time axis for the sequence
	axis(1, at=1:seql-0.5, labels=xtlab, 
		pos=-0.04, 
		## mgp=c(.5,.5,0), 
		cex.axis=cex.plot)

	## Frequency of the representative sequence
	nbprox <- sum(Quality$nbprox[1:nbrep])
	ctfreq <- round((nbprox/n)*100,1)
	text(seql/2, 1.3, 
		paste("Criterion=",ctname,", ",nbprox," neighbours (", ctfreq ,"%)", sep=""), 
		cex=cex.plot)

	## 
	legend.text <- NULL

	## ==========
	## Statistics
	## ==========

	## Start position
	y.sym.pos <- 1.8
	
	## Distance to representative sequence
	repcol <- brewer.pal(10,"Paired")
	repsymb <- c(21:25,15:19)

	dist.rep.pos <- y.sym.pos		

	## Distance to representative seq.
	for (i in 1:nbrep) {
		lines(mean(Distances[,i],na.rm=TRUE)/dist.scaling, dist.rep.pos, 
			type="b", pch=repsymb[i], lwd=3, col=repcol[i], cex=1+barw[i])
	}

	mean.central.dist <- round(mean(Distances,na.rm=TRUE),1)

	legend.B <- paste("(B) Mean dist. to representative seq.", sep="")

	## Symbols for the representative sequences
	spval <- 0.2/nbrep
	y.lab.pos <- (barw[1]/2)+spval
	lines(-1, y.lab.pos, type="b", pch=repsymb[1], lwd=3, col=repcol[1], cex=cex.plot+barw[1])
	for (p in 2:nbrep) {
		y.lab.pos <- sum(barw[1:p-1])+(p*spval)+(barw[p]/2)
		lines(-1, y.lab.pos, type="b", pch=repsymb[p], lwd=3, col=repcol[p], cex=cex.plot+barw[p])
	}	

	## Distance to center
	y.sym.pos <- y.sym.pos + 2*vspace
	dist.center.pos <- y.sym.pos

	for (i in 1:nbrep) {
		lines(Quality$inertia[i]/dist.scaling, 
			y.sym.pos, 
			type="b", pch=repsymb[i], lwd=3, col=repcol[i], cex=cex.plot+barw[i])
	}

	legend.A <- paste("(A) Discrepancy (mean dist. to center)",sep="")

	## Axis for the boxplot of distances to the representative sequences
	dypos <- 1.7

	nbdec <- if (dmax>=4) 0 else 1  
		
	axis(1, at=seq(0,seql,seql/4), 
		labels=round(seq(0,dmax,dmax/4),nbdec), 
		pos=dypos, mgp=c(.5,.5,0), 
		cex.axis=cex.plot)
		
	axis(2, at=c(dist.rep.pos,dist.center.pos), 
		labels=c("B","A"), 
		las=2, 
		cex.axis=cex.plot)

	legend(seql/2, ymax, legend=c(legend.A, legend.B),
		xjust=0.5,
		yjust=0.7,
		## title="Distance", 
		box.lty=0,
		cex=cex.plot)
}

