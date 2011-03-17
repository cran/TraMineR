## ======================================
## State distribution for each individual
## ======================================

seqistatd <- function(seqdata, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, see seqdef function to create one")
	return()
	}

	statl <- alphabet(seqdata)
	if (with.missing) 
		statl <- c(statl, attr(seqdata,"nr"))

	nbstat <- length(statl)
	nbseq <- nrow(seqdata)

	iseqtab <- matrix(nrow=nbseq, ncol=nbstat)

	colnames(iseqtab) <- statl
	rownames(iseqtab) <- rownames(seqdata)

	message(" [>] computing state distribution for ", nbseq," sequences ...")

	for (i in 1:nbstat) {
		iseqtab[,i] <- apply(seqdata,1,function(x) sum(x==statl[i],na.rm=TRUE))
	}			
	
	return(iseqtab)

	}	
