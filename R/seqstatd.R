## States frequency by time unit

seqstatd <- function(seqdata, weighted=TRUE, with.missing=FALSE, norm=TRUE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")
	
	## Retrieving the alphabet
	statl <- attr(seqdata,"alphabet")

	if (with.missing)
		statl <- c(statl, attr(seqdata,"nr"))

	nbstat <- length(statl)

	seql <- ncol(seqdata)

	sd <- matrix(nrow=nbstat,ncol=seql)
	row.names(sd) <- statl
	colnames(sd) <- colnames(seqdata)

	## Weights
	weights <- attr(seqdata, "weights")

	if (!weighted || is.null(weights)) 
		weights <- rep(1, nrow(seqdata))

	for (i in 1:nbstat)
		for (j in 1:seql)
			sd[i,j] <- sum((seqdata[,j]==statl[i])*weights)

	## sd <-	apply(seqdata,2,table)
	N <- apply(sd,2,sum)
	for (i in 1:seql) sd[,i] <- sd[,i]/N[i]

	E <- apply(sd,2,entropy)
	## Maximum entropy is the entropy of the alphabet
	if (norm==TRUE) {
		E.max <- entropy(rep(1/nbstat,nbstat)) 
		E <- E/E.max
	}

	res <- list(sd,N,E)
	names(res) <- c("Frequencies", "ValidStates", "Entropy")
	
	class(res) <- c("stslist.statd","list")

	attr(res,"nbseq") <- nrow(seqdata)
	attr(res,"cpal") <- cpal(seqdata)
	attr(res,"xtlab") <- colnames(seqdata)
	attr(res,"weighted") <- weighted

	return(res)
}
