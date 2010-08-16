## ==============
## Mean durations
## ==============

seqmeant <- function(seqdata, weighted=TRUE, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	istatd <- suppressMessages(seqistatd(seqdata, with.missing=with.missing))

	weights <- attr(seqdata, "weights")

	if (!weighted || is.null(weights)) 
		weights <- rep(1, nrow(seqdata))
	## Also takes into account that in unweighted sequence objects created with 
	## older TraMineR versions the weights attribute is a vector of 1
	## instead of NULL  
	if (all(weights==1)) 
		weighted <- FALSE

	mtime <- apply(istatd*weights,2,sum)

	res <- mtime/sum(weights)

	res <- as.matrix(res)
	colnames(res) <- "Mean"

	class(res) <- c("stslist.meant", "matrix")

	col <- cpal(seqdata)
	if (with.missing) {
		col <- c(col, attr(seqdata,"missing.color"))
	}

	attr(res,"nbseq") <- sum(weights)
	attr(res,"cpal") <- col
	attr(res,"xtlab") <- colnames(seqdata)
	attr(res,"weighted") <- weighted

	return(res)
}

