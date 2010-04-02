## =====================================
## Number of transitions in the sequence
## =====================================

seqtransn <- function(seqdata, with.missing=FALSE, norm=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	## Number of transitions
	dss <- seqdss(seqdata, with.missing=with.missing)
	trans <- seqlength(dss)-1

	if (norm) {
		seql <- seqlength(seqdata)
		trans <- trans/(seql-1)
	}

	colnames(trans) <- "Trans."

	return(trans)
}
