## ==============
## Mean durations
## ==============

seqmeant <- function(seqdata, weighted=TRUE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	istatd <- suppressMessages(seqistatd(seqdata))

	weights <- attr(seqdata, "weights")

	if (!weighted || is.null(weights)) 
		weights <- rep(1, nrow(seqdata))

	mtime <- apply(istatd*weights,2,sum)

	res <- mtime/sum(weights)

	res <- as.matrix(res)
	colnames(res) <- "Mean"

	class(res) <- c("stslist.meant", "matrix")

	attr(res,"nbseq") <- nrow(seqdata)
	attr(res,"cpal") <- cpal(seqdata)
	attr(res,"xtlab") <- colnames(seqdata)
	attr(res,"weighted") <- weighted

	return(res)
}
