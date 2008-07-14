## ======================================
## State distribution for each individual
## ======================================

seqistatd <- function(seqdata) {

	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, see seqdef function to create one")
	return()
	}

	statl <- alphabet(seqdata)
	nbstat <- length(statl)
	iseqtab <- matrix(nrow=seqdim(seqdata)[1],ncol=nbstat)

	colnames(iseqtab) <- statl

	message(" [>] Computing state distribution for ",seqdim(seqdata)[1]," sequences ...")

	for (i in 1:nbstat) {
		iseqtab[,i] <- apply(seqdata,1,function(x) sum(x==statl[i],na.rm=TRUE))
	}			
	
	return(iseqtab)

	}	
