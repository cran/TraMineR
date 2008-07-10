## =============================
## Plotting individual sequences
## =============================

seqiplot <- function(seqdata, tlim=1:10, sortv=NULL, statl=NULL, title=NULL, cpal=NULL, 
	withlegend=TRUE, withborder=TRUE, space=NULL, ltext=NULL, slab=FALSE, xtlab=NULL,
	bmar=1) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use 'seqdef' function to create one")
	
	if (is.null(statl)) statl <- attr(seqdata,"alphabet")

	if (tlim[1]==0) tlim <- 1:seqdim(seqdata)[1]

	if (!is.null(sortv)) {
		seqdata <- seqdata[order(sortv),]
		sortlab <- paste(", sorted")
	}
	else sortlab <- NULL

	nbstat <- length(statl)
	seql <- seqdim(seqdata)[2]
	ssamp <- seqdata[tlim,]
	nbseq <- seqdim(ssamp)[1]
	
	seqbar <- apply(ssamp,1,seqgbar,seql, statl,nbstat)

	if (is.null(cpal)) cpal <- c(attr(seqdata,"cpal"),"white")
	else cpal <- c(cpal,"white")

	if (withborder==TRUE) bordcol <- "black" else bordcol <- NA

	## Adding some space for the title and legend
	if (!is.null(title)) mt <- 2
	else mt <- 0
	
	if (withlegend==TRUE) ml <- 6
	else ml <- 0

	par(mar = c(bmar+3, 3, bmar+mt+ml, 2) + 0.1, xpd=TRUE)

	sppar <- space

	## The PLot
	barplot(seqbar,col=cpal,
		names.arg=tlim,
		ylab=paste("Sequences ",min(tlim),"-",max(tlim),sortlab,sep=""),
		mgp=c(1.5,0,0),
		main=title,
		horiz=TRUE,
		border=bordcol,
		yaxt="n",
		axes=FALSE,
		las=1, 
		space=sppar
		)

	tstart <-  attr(seqdata,"start")

	if (is.null(xtlab)) xtlab <- colnames(seqdata)
	axis(1,at=1:seql-0.5,labels=xtlab,mgp=c(3,0.5,0))

	if (slab==TRUE) {
		for (i in 1:nbseq) 
			for (j in 1:seql)  
				text(j,i+(i*0.2),ssamp[i,j],cex=1.5,adj=c(1.5,1.3))
	}
	
	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	## Computing some parameters for the legend's plotting
	leg.ncol <- if (round(nbstat/3,0)>1)  round(nbstat/3,0) else 2
	leg.inset <- -0.2 + ((2-leg.ncol)*0.025)

	if (withlegend==TRUE) 
		legend("top",inset=c(0,leg.inset),
			legend=ltext,
			fill=cpal,
			ncol=leg.ncol,
			bty="n")
}
		
