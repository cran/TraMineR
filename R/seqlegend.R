## =========================================
## Plotting the legend for a sequence object
## =========================================

seqlegend <- function(seqdata, cpal, ltext, 
	position="topleft", fontsize=1,...) {
	
	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	if (missing(cpal)) cpal <- attr(seqdata,"cpal")

	if (missing(ltext)) ltext <- attr(seqdata,"labels")

	plot(0, type= "n", axes=FALSE, xlab="", ylab="")
	legend(position, fill=cpal, legend=ltext, cex=fontsize,...)
}
