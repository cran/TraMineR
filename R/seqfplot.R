## ================================
## PLot of the sequences frequency
## ================================

seqfplot <- function(seqdata, tlim=10, title=NULL, cpal=NULL, pbarw=FALSE, 
	withlegend=TRUE, ltext=NULL, xtlab=NULL,bmar=1) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	n <- seqdim(seqdata)[1]

	tab <- seqtab(seqdata,tlim)
	seqlist <- suppressMessages(seqformat(row.names(tab),from="SPS2",to="STS"))
	seqlist <- seqdecomp(seqlist)

	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)
	seql <- seqdim(seqlist)[2]

	seqbar <- matrix(0,nrow=nrow(seqlist),ncol=seql*(nbstat+1))
	
	for (i in 1:tlim) {
		for (j in 1:seql) { 
			if (!is.na(seqlist[i,j])) seqbar[i,((j-1)*(nbstat+1))+which(statl==seqlist[i,j])] <- 1
			else seqbar[i,((j-1)*(nbstat+1))+nbstat+1] <- 1
		}
	}

	if (is.null(cpal)) cpal <- c(attr(seqdata,"cpal"),"white")
	else cpal <- c(cpal,"white")

	if (pbarw==TRUE) barw=tab$Percent else barw=1

	## allowing the legend to be plotted outside the plot region
	if (!is.null(title)) mt <- 2
	else mt <- 0
	
	if (withlegend==TRUE) ml <- 6
	else ml <- 0

	par(mar = c(bmar+3, 4, bmar+mt+ml, 2) + 0.1, xpd=TRUE)

	barplot(t(seqbar),col=cpal,width=barw,
		names.arg=paste(round(tab$Percent,1),"%",sep=""),
		cex.names=0.9,
		ylab=paste("Freq. (n=",n,")",sep=""),
		mgp=c(2.5,0,0),
		main=title,
		horiz=TRUE,
		axes=FALSE,
		las=1)
	
	if (is.null(xtlab)) xtlab <- colnames(seqdata)[1:seql]
	axis(1,at=1:seql-0.5,labels=xtlab,mgp=c(2.5,0.5,0))

	## Computing some parameters for the legend's plotting
	leg.ncol <- if (round(nbstat/3,0)>1)  round(nbstat/3,0) else 2
	leg.inset <- -0.2 + ((2-leg.ncol)*0.025)

	## Plotting the legend	
	if (is.null(ltext)) ltext <- attr(seqdata,"labels")

	if (withlegend==TRUE) 
		legend("top",inset=c(0,leg.inset),
			legend=ltext,
			fill=cpal,
			ncol=leg.ncol,
			bty="n")
	}
