## ============================================
## PLOT OF THE STATES DISTRIBUTION BY TIME UNIT
## ============================================

seqdplot <- function(seqdata, cpal=NULL, title=NULL, 
	withlegend=TRUE, withborder=TRUE, ltext=NULL, xtlab=NULL, bmar=1) {

	if (!inherits(seqdata,"stslist")) {
		stop(" [!] data is not a sequence object, see seqdef function to create one\n")
		return()
		}

	n <- seqdim(seqdata)[1]
	lseq <- seqdim(seqdata)[2]

	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)
	
	tseq <- seqstatd(seqdata, digits=NULL)$Frequencies
		
	if (is.null(xtlab)) xtlab <- colnames(seqdata)
	colnames(tseq) <- xtlab
	
	## ELIMINATION DES ETATS
	nseq <- tseq
	## for (i in 1:nrow(tseq)) if (sum(tseq[i,])==0) nseq <- nseq[-i,]

	if (is.null(cpal)) cpal <- attr(seqdata,"cpal")
	else cpal <- cpal

	## Adding some space for the title and legend
	if (!is.null(title)) mt <- 2
	else mt <- 0
	
	if (withlegend==TRUE) ml <- 6
	else ml <- 0

	par(mar = c(bmar+3, 4, bmar+1+mt+ml, 2) + 0.1, xpd=TRUE)

	if (withborder==TRUE) bordcol <- "black" else bordcol <- NA

	barplot(nseq,
		space=0,
		mgp=c(2.5,0.6,0),
		cex.names=0.9,
		col=cpal,
		border=bordcol,
		main=title,
		ylab=paste("Freq. (n=",n,")",sep=""))

	## axis(1,at=1:lseq-0.5,mgp=c(3,0.5,0))

	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	## Computing some parameters for the legend's plotting
	leg.ncol <- if (round(nbstat/3,0)>1)  round(nbstat/3,0) else 2
	leg.inset <- -0.22 + ((2-leg.ncol)*0.025)

	if (withlegend==TRUE) 
		legend("top",inset=c(0,leg.inset),
			legend=ltext,
			fill=cpal,
			ncol=leg.ncol,
			bty="n")

	}
