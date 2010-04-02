## =======================================
## Extracts distinct states from sequences
## =======================================

seqdss <- function(seqdata, with.missing=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	nbseq <- nrow(seqdata)

	sl <- seqlength(seqdata) 
	maxsl <- max(sl)

	trans <- matrix(nrow=nbseq, ncol=maxsl)
	statl <- attr(seqdata,"alphabet")

	seqdatanum <- TraMineR:::seqasnum(seqdata, with.missing=with.missing)

	if (!with.missing)
		seqdatanum[is.na(seqdatanum)] <- -99

	for (i in 1:nbseq) {
		idx <- 1
		j <- 1

		tmpseq <- seqdatanum[i,]
		
		while (idx <= sl[i]) {
			iseq <- tmpseq[idx]

			while (idx < sl[i] & (tmpseq[idx+1]==iseq || tmpseq[idx+1]==-99)) { 
				idx <- idx+1
			}

			## The range of the numeric alphabet 
			## obtained with seqasnum is 0..n
			if (iseq!=-99) {
				trans[i,j] <- statl[(iseq+1)]
				j <- j+1
			}
			idx <- idx+1
		}
	}
	
	trans <- suppressMessages(
		seqdef(trans, alphabet=statl, 
		cnames=paste("ST",seq(1:maxsl),sep=""), 
		cpal=cpal(seqdata),
		id=rownames(seqdata)))

	return(trans)
}

