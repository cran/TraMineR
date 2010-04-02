## ====================
## Modal state sequence
## ====================

seqmodst <- function(seqdata, weighted=TRUE, with.missing=FALSE, dist=FALSE, ...) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one", call.=FALSE)

	slength <- ncol(seqdata)
	statelist <- alphabet(seqdata)
	cnames <- colnames(seqdata)

	## State distribution
	freq <- seqstatd(seqdata, weighted, with.missing)$Frequencies

	ctype <- matrix(nrow=1, ncol=slength)
	stfreq <- matrix(nrow=1, ncol=slength)
	colnames(stfreq) <- cnames
	rownames(stfreq) <- "Freq."
	colnames(ctype) <- cnames

	## Constructing the transversal modal sequence
	for (i in 1:slength) {
		smax <- which(freq[,i]==max(freq[,i]))[1]
		stfreq[,i] <- freq[smax,i]
		ctype[,i] <- statelist[smax]
	}
	
	res <- suppressMessages(seqdef(ctype, alphabet=alphabet(seqdata), 
		labels=attr(seqdata,"labels"), cpal=cpal(seqdata)))

	nbocc <- length(seqfind(res, seqdata))

	## Distance to modal state sequence
	if (dist)
		dist.modst <- seqdist(seqdata, refseq=res, ...)
	else 
		dist.modst <- NULL

	class(res) <- c("stslist.modst", class(res))

	attr(res, "Frequencies") <- stfreq
	attr(res, "nbseq") <- nrow(seqdata)
	attr(res, "Distances") <- dist.modst
	attr(res, "Occurences") <- nbocc

	return(res)
 } 

